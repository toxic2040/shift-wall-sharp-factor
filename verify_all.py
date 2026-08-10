#!/usr/bin/env python3
"""Run the complete Paper 4 replay from an isolated, resumable work tree."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import subprocess
import sys
import threading
import time


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "SHA256SUMS"
CLASSIFICATION = "REPRODUCIBLE_MANUSCRIPT_REPLAY"
EXPECTED_SYMPY = "1.14.0"
WORKERS = 6
MANIFEST_LINE = re.compile(r"([0-9a-f]{64})  ([!-~]+)")
TEXT_SUFFIXES = {"", ".cff", ".gitignore", ".json", ".jsonl", ".md", ".py", ".tex", ".txt"}

BRIDGE_REL = Path("verification/lower_stripes")
ALLT_REL = Path("verification/all_stripes")
H2_REL = ALLT_REL / "columns"
H3_REL = ALLT_REL / "core"
MAIN_REL = ALLT_REL / "join"
VENDOR_REL = Path("vendor/odd-cycles-square-tail-replay")

BRIDGE_OUTPUTS = (
    BRIDGE_REL / "BRIDGE_T2_MINIMUM.json",
    BRIDGE_REL / "BRIDGE_T3_BARRIER.json",
    BRIDGE_REL / "SQ_T4_BRIDGE_TRANSFER.json",
    BRIDGE_REL / "A1_PI_ENDPOINT_AUDIT.json",
)
CORE_OUTPUTS = (
    H3_REL / "A_core_sweep.jsonl",
    H3_REL / "A_CORE_SWEEP_SUMMARY.json",
)
PREFIX_OUTPUTS = tuple(H2_REL / f"H2_PREFIX_T{t}.json" for t in range(5, 9))
TAIL_SPECS = (
    ("t10_m200000", 10, 200_000, "PASS", 0),
    ("t10_m100000", 10, 100_000, "FAIL", 1),
    ("t9_m68133", 9, 68_133, "PASS", 0),
    ("t9_m55000", 9, 55_000, "FAIL", 1),
    ("t8_m19999", 8, 19_999, "PASS", 0),
    ("t7_m12398", 7, 12_398, "PASS", 0),
    ("t7_m6199", 7, 6_199, "FAIL", 1),
    ("t6_m4006", 6, 4_006, "PASS", 0),
    ("t6_m2003", 6, 2_003, "FAIL", 1),
    ("t5_m1000", 5, 1_000, "PASS", 0),
    ("t4_m1000_crosscheck", 4, 1_000, "PASS", 0),
)
TAIL_OUTPUTS = tuple(
    H2_REL / (
        "H2_TAIL_T4_M1000_crosscheck.json"
        if key == "t4_m1000_crosscheck"
        else f"H2_TAIL_T{t}_M{m0}.json"
    )
    for key, t, m0, _status, _code in TAIL_SPECS
)
H2_JOIN_OUTPUTS = tuple(H2_REL / f"H2_COL_T{t}.json" for t in range(5, 11))
OTHER_OUTPUTS = (
    H2_REL / "H2_A_GATE.json",
    MAIN_REL / "verify_joint_core_t11.json",
    MAIN_REL / "verify_joint_core_bidegree.json",
    MAIN_REL / "verify_joint_core_sharp.json",
    MAIN_REL / "dlaw_leg1_prefix_anchor.json",
    MAIN_REL / "DLAW_LEG2.json",
    MAIN_REL / "DLAW_LEG2_CONTROL.json",
    MAIN_REL / "ALLT_JOIN.json",
)
PROGRESS_OUTPUTS = (
    H2_REL / "h2_prefix_progress.jsonl",
    MAIN_REL / "dlaw_leg1_prefix_progress.jsonl",
)
MUTABLE_OUTPUTS = frozenset(
    path.as_posix()
    for path in (
        *BRIDGE_OUTPUTS,
        *CORE_OUTPUTS,
        *PREFIX_OUTPUTS,
        *TAIL_OUTPUTS,
        *H2_JOIN_OUTPUTS,
        *OTHER_OUTPUTS,
    )
)
PURGED_OUTPUTS = (
    *BRIDGE_OUTPUTS,
    *CORE_OUTPUTS,
    *PREFIX_OUTPUTS,
    *TAIL_OUTPUTS,
    *H2_JOIN_OUTPUTS,
    *OTHER_OUTPUTS,
    *PROGRESS_OUTPUTS,
)
BYTE_COMPARE_OUTPUTS = (
    *BRIDGE_OUTPUTS,
    CORE_OUTPUTS[1],
    *PREFIX_OUTPUTS,
    *TAIL_OUTPUTS,
    *H2_JOIN_OUTPUTS,
    *OTHER_OUTPUTS,
)

# Encoded so the scanner does not whitelist itself by containing the text.
FORBIDDEN_PATH_MARKERS = tuple(
    bytes.fromhex(item)
    for item in ("2f686f6d652f", "2f55736572732f", "2e2e2f", "7e2f")
)
FORBIDDEN_AUTHORSHIP_MARKERS = tuple(
    bytes.fromhex(item)
    for item in (
        "6f70656e6169",
        "63686174677074",
        "636f646578",
        "616e7468726f706963",
        "636c61756465",
        "67656d696e69",
        "67726f6b",
        "67656e657261746564206279",
        "616920617373697374616e74",
    )
)

ACTIVE_PROCESSES: dict[int, subprocess.Popen[str]] = {}
ACTIVE_LOCK = threading.Lock()


class ReplayError(RuntimeError):
    """A fail-closed replay gate did not pass."""


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def parse_manifest(text: str) -> dict[str, str]:
    if not text or not text.endswith("\n"):
        raise ReplayError("SHA256SUMS is empty or lacks a final newline")
    entries: dict[str, str] = {}
    for number, line in enumerate(text.splitlines(), 1):
        match = MANIFEST_LINE.fullmatch(line)
        if match is None:
            raise ReplayError(f"malformed SHA256SUMS line {number}")
        digest, name = match.groups()
        pure = PurePosixPath(name)
        if pure.is_absolute() or ".." in pure.parts or name != pure.as_posix():
            raise ReplayError(f"unsafe SHA256SUMS path on line {number}: {name}")
        if name in entries:
            raise ReplayError(f"duplicate SHA256SUMS path: {name}")
        entries[name] = digest
    if list(entries) != sorted(entries):
        raise ReplayError("SHA256SUMS paths are not sorted")
    return entries


def checkout_inventory() -> set[str]:
    inventory: set[str] = set()
    for path in HERE.rglob("*"):
        relative = path.relative_to(HERE)
        if relative.parts and relative.parts[0] == ".git":
            continue
        if path.is_symlink() or path.is_file():
            inventory.add(relative.as_posix())
    return inventory


def validate_release_tree() -> tuple[dict[str, str], str]:
    if not MANIFEST.is_file():
        raise ReplayError("SHA256SUMS is missing")
    manifest_text = MANIFEST.read_text(encoding="ascii")
    entries = parse_manifest(manifest_text)
    if MANIFEST.name in entries:
        raise ReplayError("SHA256SUMS must not contain its own digest")
    allowed = set(entries) | {MANIFEST.name}
    actual = checkout_inventory()
    missing = sorted(set(entries) - actual)
    extra = sorted(actual - allowed)
    if missing:
        raise ReplayError(f"manifested files are missing: {missing[:8]}")
    if extra:
        raise ReplayError(f"unmanifested files are present: {extra[:8]}")
    failures: list[str] = []
    for name, expected in entries.items():
        path = HERE / name
        if path.is_symlink() or not path.is_file():
            failures.append(f"unsafe or missing path: {name}")
        elif sha256_file(path) != expected:
            failures.append(f"hash mismatch: {name}")
    if failures:
        raise ReplayError("; ".join(failures[:8]))
    scan_portability(entries)
    return entries, sha256_bytes(manifest_text.encode("ascii"))


def scan_portability(entries: dict[str, str]) -> None:
    failures: list[str] = []
    for name in entries:
        path = HERE / name
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        payload = path.read_bytes().lower()
        if any(marker in payload for marker in FORBIDDEN_PATH_MARKERS):
            failures.append(f"nonportable path marker: {name}")
        if any(marker in payload for marker in FORBIDDEN_AUTHORSHIP_MARKERS):
            failures.append(f"external-artifact authorship marker: {name}")
    if failures:
        raise ReplayError("; ".join(failures[:8]))


def validate_environment() -> None:
    failures: list[str] = []
    if sys.flags.optimize != 0:
        failures.append("Python optimization must be disabled")
    if not sys.flags.isolated:
        failures.append("invoke Python with -I")
    if not sys.flags.dont_write_bytecode:
        failures.append("invoke Python with -B")
    if sys.version_info[:2] != (3, 14):
        failures.append("Python 3.14 is required")
    is_gil_enabled = getattr(sys, "_is_gil_enabled", None)
    if is_gil_enabled is None or is_gil_enabled():
        failures.append("the free-threading Python build is required")
    try:
        installed_sympy = version("sympy")
    except PackageNotFoundError:
        installed_sympy = "missing"
    if installed_sympy != EXPECTED_SYMPY:
        failures.append(
            f"SymPy {EXPECTED_SYMPY} is required, found {installed_sympy}"
        )
    if failures:
        raise ReplayError("; ".join(failures))


def ignored_copy_path(_directory: str, names: list[str]) -> set[str]:
    ignored = {name for name in names if name in {".git", "tmp", "__pycache__"}}
    ignored.update(name for name in names if name.endswith(".pyc"))
    return ignored


def load_state(path: Path) -> dict:
    if not path.exists():
        return {"schema": "shift-wall-replay-state-v1", "completed_stages": {}}
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReplayError(f"invalid replay state: {error}") from error
    if state.get("schema") != "shift-wall-replay-state-v1":
        raise ReplayError("unrecognized replay state schema")
    return state


def setup_work_tree(
    work_dir: Path, manifest_entries: dict[str, str], manifest_sha: str
) -> tuple[Path, Path, dict]:
    work_dir = work_dir.resolve()
    if work_dir == HERE or work_dir.is_relative_to(HERE):
        raise ReplayError("--work-dir must be outside the release checkout")
    work_dir.mkdir(parents=True, exist_ok=True)
    tree = work_dir / "tree"
    state_path = work_dir / "state.json"
    state = load_state(state_path)
    prior_sha = state.get("release_manifest_sha256")
    if prior_sha is not None and prior_sha != manifest_sha:
        raise ReplayError("work directory belongs to a different release manifest")

    if not tree.exists():
        if any(work_dir.iterdir()):
            allowed = {state_path.name}
            unexpected = sorted(path.name for path in work_dir.iterdir()
                                if path.name not in allowed)
            if unexpected:
                raise ReplayError(
                    f"work directory is not empty and has no replay tree: {unexpected}"
                )
        shutil.copytree(HERE, tree, ignore=ignored_copy_path)
        for relative in PURGED_OUTPUTS:
            (tree / relative).unlink(missing_ok=True)
        state.update({
            "release_manifest_sha256": manifest_sha,
            "workers": WORKERS,
            "initialized": True,
            "completed_stages": {},
        })
        atomic_json(state_path, state)
        print(f"[INIT] isolated replay tree: {tree}")
    elif prior_sha is None:
        raise ReplayError("existing replay tree has no matching state record")

    validate_static_work_tree(tree, manifest_entries)
    (work_dir / "logs").mkdir(exist_ok=True)
    return tree, state_path, state


def validate_static_work_tree(tree: Path, manifest_entries: dict[str, str]) -> None:
    failures: list[str] = []
    for name, expected in manifest_entries.items():
        if name in MUTABLE_OUTPUTS:
            continue
        path = tree / name
        if path.is_symlink() or not path.is_file():
            failures.append(f"missing or unsafe static replay file: {name}")
        elif sha256_file(path) != expected:
            failures.append(f"static replay file changed: {name}")
    if failures:
        raise ReplayError("; ".join(failures[:8]))


def terminate_active(_signum: int | None = None, _frame=None) -> None:
    with ACTIVE_LOCK:
        processes = list(ACTIVE_PROCESSES.values())
    for process in processes:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


def command_for(script: Path, *arguments: str) -> list[str]:
    return [sys.executable, "-I", "-B", str(script), *arguments]


def stage_is_complete(state: dict, name: str, tree: Path) -> bool:
    record = state.get("completed_stages", {}).get(name)
    if not isinstance(record, dict) or record.get("status") != "PASS":
        return False
    for relative, expected in record.get("outputs", {}).items():
        path = tree / relative
        if not path.is_file() or sha256_file(path) != expected:
            return False
    return True


def save_stage(
    state: dict, state_path: Path, name: str, tree: Path, outputs: tuple[Path, ...]
) -> None:
    digests = {}
    for relative in outputs:
        path = tree / relative
        if not path.is_file():
            raise ReplayError(f"stage {name} did not create {relative}")
        digests[relative.as_posix()] = sha256_file(path)
    state.setdefault("completed_stages", {})[name] = {
        "status": "PASS",
        "outputs": digests,
    }
    atomic_json(state_path, state)


def run_logged(command: list[str], cwd: Path, log_path: Path, quiet: bool = False) -> int:
    environment = os.environ.copy()
    environment.update({"PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1"})
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            bufsize=1,
            start_new_session=True,
        )
        with ACTIVE_LOCK:
            ACTIVE_PROCESSES[process.pid] = process
        try:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                log.flush()
                if not quiet:
                    print(line, end="", flush=True)
            return process.wait()
        finally:
            with ACTIVE_LOCK:
                ACTIVE_PROCESSES.pop(process.pid, None)


def run_stage(
    *,
    name: str,
    command: list[str],
    cwd: Path,
    tree: Path,
    work_dir: Path,
    state: dict,
    state_path: Path,
    outputs: tuple[Path, ...] = (),
    expected_codes: tuple[int, ...] = (0,),
) -> None:
    if stage_is_complete(state, name, tree):
        print(f"[RESUME] {name}: already complete")
        return
    print(f"\n=== {name} ===", flush=True)
    code = run_logged(command, cwd, work_dir / "logs" / f"{name}.log")
    if code not in expected_codes:
        raise ReplayError(f"stage {name} exited {code}, expected {expected_codes}")
    save_stage(state, state_path, name, tree, outputs)


def append_jsonl(path: Path, payload: dict) -> None:
    with path.open("a", encoding="utf-8") as sink:
        sink.write(json.dumps(payload, sort_keys=True) + "\n")
        sink.flush()
        os.fsync(sink.fileno())


def load_tail_progress(path: Path) -> dict[str, dict]:
    records: dict[str, dict] = {}
    if not path.exists():
        return records
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, 1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ReplayError(
                    f"malformed tail progress line {line_number}: {error}"
                ) from error
            if isinstance(record, dict) and isinstance(record.get("key"), str):
                records[record["key"]] = record
    return records


def valid_tail_record(record: dict | None, tree: Path) -> bool:
    if not record or record.get("success") is not True:
        return False
    relative = record.get("output")
    expected = record.get("output_sha256")
    if not isinstance(relative, str) or not isinstance(expected, str):
        return False
    path = tree / relative
    return path.is_file() and sha256_file(path) == expected


def run_tail_task(
    spec: tuple[str, int, int, str, int], tree: Path, work_dir: Path
) -> dict:
    key, t, m0, expected_status, expected_code = spec
    h2 = tree / H2_REL
    output_relative = (
        H2_REL / "H2_TAIL_T4_M1000_crosscheck.json"
        if key == "t4_m1000_crosscheck"
        else H2_REL / f"H2_TAIL_T{t}_M{m0}.json"
    )
    output = tree / output_relative
    output.unlink(missing_ok=True)
    command = command_for(h2 / "h2_tail_cert.py", str(t), str(m0),
                          "--output", str(output))
    started = time.monotonic()
    try:
        code = run_logged(
            command,
            h2,
            work_dir / "logs" / f"tail_{key}.log",
            quiet=True,
        )
        payload = (
            json.loads(output.read_text(encoding="utf-8"))
            if output.exists()
            else {}
        )
        actual_status = payload.get("status")
        success = code == expected_code and actual_status == expected_status
        return {
            "key": key,
            "t": t,
            "m0": m0,
            "expected_exit": expected_code,
            "actual_exit": code,
            "expected_status": expected_status,
            "actual_status": actual_status,
            "output": output_relative.as_posix(),
            "output_sha256": sha256_file(output) if output.is_file() else None,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "success": success,
        }
    except Exception as error:  # one failed attempt must not kill the batch
        return {
            "key": key,
            "t": t,
            "m0": m0,
            "expected_exit": expected_code,
            "actual_exit": None,
            "expected_status": expected_status,
            "actual_status": None,
            "output": output_relative.as_posix(),
            "output_sha256": None,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "success": False,
            "error": repr(error),
        }


def run_tail_batch(tree: Path, work_dir: Path) -> None:
    progress_path = work_dir / "tail_progress.jsonl"
    recorded = load_tail_progress(progress_path)
    pending = [spec for spec in TAIL_SPECS
               if not valid_tail_record(recorded.get(spec[0]), tree)]
    if not pending:
        print("[RESUME] fixed-column tails: all 11 attempts already complete")
        return
    print(
        f"\n=== fixed-column tails: {len(pending)} pending attempts, "
        f"{WORKERS} concurrent single-process workers ===",
        flush=True,
    )
    failures: list[dict] = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {
            pool.submit(run_tail_task, spec, tree, work_dir): spec
            for spec in pending
        }
        for future in as_completed(futures):
            spec = futures[future]
            try:
                record = future.result()
            except Exception as error:  # defensive boundary around each entry
                key, t, m0, expected_status, expected_code = spec
                record = {
                    "key": key,
                    "t": t,
                    "m0": m0,
                    "expected_exit": expected_code,
                    "actual_exit": None,
                    "expected_status": expected_status,
                    "actual_status": None,
                    "output": None,
                    "output_sha256": None,
                    "success": False,
                    "error": repr(error),
                }
            append_jsonl(progress_path, record)
            label = "PASS" if record["success"] else "FAIL"
            print(
                f"[{label}] tail {record['key']}: exit={record['actual_exit']} "
                f"status={record['actual_status']} "
                f"elapsed={record.get('elapsed_seconds')}s",
                flush=True,
            )
            if not record["success"]:
                failures.append(record)
    if failures:
        names = [record["key"] for record in failures]
        raise ReplayError(f"tail attempts did not match their contracts: {names}")


def load_core_records(path: Path) -> dict[int, dict]:
    records: dict[int, dict] = {}
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, 1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ReplayError(
                    f"bad core JSONL line {line_number}: {error}"
                ) from error
            if not isinstance(record, dict) or not isinstance(record.get("r"), int):
                raise ReplayError(f"bad core record on line {line_number}")
            records[record["r"]] = record
    return records


def compare_core_jsonl(tree: Path) -> None:
    banked = load_core_records(HERE / CORE_OUTPUTS[0])
    replayed = load_core_records(tree / CORE_OUTPUTS[0])
    expected_rows = set(range(7, 504))
    if set(banked) != expected_rows or set(replayed) != expected_rows:
        raise ReplayError("core JSONL does not contain exactly rows 7 through 503")
    for r in sorted(expected_rows):
        left = dict(banked[r])
        right = dict(replayed[r])
        left.pop("secs", None)
        right.pop("secs", None)
        if left != right:
            raise ReplayError(f"core JSONL semantic mismatch at row r={r}")
    print("[PASS] finite-core JSONL: 497 rows match modulo reporting-only time")


def compare_banked_outputs(tree: Path) -> dict[str, str]:
    compared: dict[str, str] = {}
    for relative in BYTE_COMPARE_OUTPUTS:
        banked = HERE / relative
        replayed = tree / relative
        if not banked.is_file() or not replayed.is_file():
            raise ReplayError(f"missing banked or replayed output: {relative}")
        banked_bytes = banked.read_bytes()
        replayed_bytes = replayed.read_bytes()
        if replayed_bytes != banked_bytes:
            raise ReplayError(f"byte mismatch in regenerated output: {relative}")
        compared[relative.as_posix()] = sha256_bytes(banked_bytes)
    compare_core_jsonl(tree)
    print(f"[PASS] {len(compared)} deterministic outputs reproduced byte-for-byte")
    return compared


def execute_replay(tree: Path, work_dir: Path, state: dict, state_path: Path) -> dict:
    vendor = tree / VENDOR_REL
    bridge = tree / BRIDGE_REL
    h2 = tree / H2_REL
    h3 = tree / H3_REL
    main_dir = tree / MAIN_REL

    run_stage(
        name="upstream_panel3_replay",
        command=command_for(vendor / "verify_all.py"),
        cwd=vendor,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
    )

    bridge_stages = (
        ("bridge_t2_minimum", "verify_bridge_t2_minimum.py", BRIDGE_OUTPUTS[0]),
        ("bridge_t3_barrier", "verify_bridge_t3_barrier.py", BRIDGE_OUTPUTS[1]),
        ("bridge_t4_transfer", "verify_sq_t4_transfer.py", BRIDGE_OUTPUTS[2]),
        ("a1_pi_endpoint_audit", "audit_a1_pi_endpoint.py", BRIDGE_OUTPUTS[3]),
    )
    for name, script, output in bridge_stages:
        run_stage(
            name=name,
            command=command_for(bridge / script),
            cwd=bridge,
            tree=tree,
            work_dir=work_dir,
            state=state,
            state_path=state_path,
            outputs=(output,),
        )

    run_stage(
        name="finite_core_123753",
        command=command_for(h3 / "A_core_sweep.py"),
        cwd=h3,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
        outputs=CORE_OUTPUTS,
    )
    run_stage(
        name="finite_core_independent_crosscheck",
        command=command_for(h3 / "A_crosscheck.py"),
        cwd=h3,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
    )

    prefix_runs = (
        ("prefix_t5_run", ("--plan", "5:1001:1", "--workers", "6")),
        ("prefix_t6_run_1", ("--plan", "6:2004:3", "--workers", "6")),
        ("prefix_t6_run_2", ("--r-start", "2004", "--plan", "6:4007:2", "--workers", "6")),
        ("prefix_t7_run_1", ("--plan", "7:6200:5", "--workers", "6")),
        ("prefix_t7_run_2", ("--r-start", "6200", "--plan", "7:12399:4", "--workers", "6")),
        ("prefix_t8_run", ("--plan", "8:20000:12", "--workers", "6")),
    )
    for name, arguments in prefix_runs:
        run_stage(
            name=name,
            command=command_for(h2 / "h2_prefix.py", *arguments),
            cwd=h2,
            tree=tree,
            work_dir=work_dir,
            state=state,
            state_path=state_path,
        )
    for t, m0 in ((5, 1000), (6, 4006), (7, 12398), (8, 19999)):
        output = H2_REL / f"H2_PREFIX_T{t}.json"
        run_stage(
            name=f"prefix_t{t}_assemble",
            command=command_for(h2 / "h2_prefix.py", "--assemble", f"{t}:{m0}"),
            cwd=h2,
            tree=tree,
            work_dir=work_dir,
            state=state,
            state_path=state_path,
            outputs=(output,),
        )

    run_stage(
        name="a_gate",
        command=command_for(h2 / "h2_a_gate.py"),
        cwd=h2,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
        outputs=(H2_REL / "H2_A_GATE.json",),
    )
    run_tail_batch(tree, work_dir)

    join_specs = (
        (5, 1000, "1000", False),
        (6, 4006, "2003,4006", False),
        (7, 12398, "6199,12398", False),
        (8, 19999, "19999", False),
        (9, 68133, "55000,68133", True),
        (10, 200000, "100000,200000", True),
    )
    for t, m0, attempts, prefix_free in join_specs:
        script = "h2_join_prefix_free.py" if prefix_free else "h2_join.py"
        run_stage(
            name=f"column_t{t}_join",
            command=command_for(h2 / script, str(t), str(m0),
                                "--attempts", attempts),
            cwd=h2,
            tree=tree,
            work_dir=work_dir,
            state=state,
            state_path=state_path,
            outputs=(H2_REL / f"H2_COL_T{t}.json",),
        )

    for name, script, output in (
        (
            "joint_t11_t16",
            "verify_joint_core_t11.py",
            MAIN_REL / "verify_joint_core_t11.json",
        ),
        (
            "joint_uniform_t17",
            "verify_joint_core_bidegree.py",
            MAIN_REL / "verify_joint_core_bidegree.json",
        ),
        (
            "dlaw_leg1",
            "dlaw_leg1_prefix.py",
            MAIN_REL / "dlaw_leg1_prefix_anchor.json",
        ),
        ("dlaw_leg2", "dlaw_leg2_bootstrap.py", MAIN_REL / "DLAW_LEG2.json"),
    ):
        run_stage(
            name=name,
            command=command_for(main_dir / script),
            cwd=main_dir,
            tree=tree,
            work_dir=work_dir,
            state=state,
            state_path=state_path,
            outputs=(output,),
        )
    run_stage(
        name="dlaw_leg2_negative_control",
        command=command_for(main_dir / "dlaw_leg2_bootstrap.py", "--negative-control"),
        cwd=main_dir,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
        outputs=(MAIN_REL / "DLAW_LEG2_CONTROL.json",),
    )
    run_stage(
        name="joint_sharpened_corridor",
        command=command_for(main_dir / "verify_joint_core_sharp.py"),
        cwd=main_dir,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
        outputs=(MAIN_REL / "verify_joint_core_sharp.json",),
    )
    run_stage(
        name="master_join_22",
        command=command_for(main_dir / "verify_allt_join.py"),
        cwd=main_dir,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
        outputs=(MAIN_REL / "ALLT_JOIN.json",),
    )

    compared = compare_banked_outputs(tree)
    run_stage(
        name="release_anchors_40",
        command=command_for(tree / "verify_release_anchors.py"),
        cwd=tree,
        tree=tree,
        work_dir=work_dir,
        state=state,
        state_path=state_path,
    )
    return compared


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="persistent directory outside the release checkout",
    )
    args = parser.parse_args()
    started = time.monotonic()
    work_dir = args.work_dir.resolve()
    report_path = work_dir / "REPLAY_REPORT.json"
    report: dict = {
        "schema": "shift-wall-full-replay-report-v1",
        "classification": CLASSIFICATION,
        "status": "FAIL",
        "workers": WORKERS,
        "interpreter": sys.version,
        "sympy": None,
    }
    try:
        validate_environment()
        report["sympy"] = version("sympy")
        manifest_entries, manifest_sha = validate_release_tree()
        report["release_manifest_sha256"] = manifest_sha
        report["manifested_files"] = len(manifest_entries)
        tree, state_path, state = setup_work_tree(
            work_dir, manifest_entries, manifest_sha
        )
        compared = execute_replay(tree, work_dir, state, state_path)
        report.update({
            "status": "PASS",
            "banked_outputs_reproduced_byte_for_byte": len(compared),
            "finite_core_rows_semantically_reproduced": 497,
            "tail_attempts": len(TAIL_SPECS),
            "tail_expected_failures_reproduced": 4,
            "tail_expected_passes_reproduced": 7,
            "master_join_gates": 22,
            "release_anchor_gates": 40,
            "elapsed_seconds": round(time.monotonic() - started, 3),
        })
        atomic_json(report_path, report)
        print(
            f"\n[PASS] complete Paper 4 replay; report={report_path}; "
            f"elapsed={report['elapsed_seconds']}s"
        )
        return 0
    except (ReplayError, OSError, subprocess.SubprocessError) as error:
        report.update({
            "error": str(error),
            "elapsed_seconds": round(time.monotonic() - started, 3),
        })
        try:
            atomic_json(report_path, report)
        except OSError:
            pass
        print(f"\n[FAIL] {error}", file=sys.stderr)
        return 1
    finally:
        terminate_active()


def handle_signal(signum: int, frame) -> None:
    terminate_active(signum, frame)
    raise KeyboardInterrupt


if __name__ == "__main__":
    for handled_signal in (signal.SIGINT, signal.SIGTERM):
        signal.signal(handled_signal, handle_signal)
    raise SystemExit(main())
