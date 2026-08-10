#!/usr/bin/env python3
"""Run the complete square-tail manuscript replay with fail-closed checks."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path, PurePosixPath
import re
import signal
import subprocess
import sys
import tempfile
import threading


HERE = Path(__file__).resolve().parent
VERIFY = HERE / "verification"
RESULTS = HERE / "results"
MANIFEST = HERE / "SHA256SUMS"
CLASSIFICATION = "REPRODUCIBLE_MANUSCRIPT_REPLAY"
STATIC_TIMEOUT_SECONDS = 1_800
GENERATOR_TIMEOUT_SECONDS = 7_200

STATIC_COMPONENTS = (
    ("wall/ladder reduction", "A1_reduce_check.py", ()),
    ("profile connection", "A8_conn2.py", ()),
    ("profile and first-lowering gates", "A8_l_curvature.py", ()),
    ("coefficient-ratio signs", "A12_ratio_sign.py", ()),
    ("uniform decrement bound", "A13_delta5.py", ()),
    ("linear extremality", "A10_linear_extremal.py", ()),
    ("joint window", "A14_uniform_sc.py", ()),
    ("shift-wall transcription", "A36_bw_bridge.py", ()),
    ("t=1 stripe characterization", "A45_t1_stripe.py", ()),
)

TAILS = (
    (2, 312, 22, True),
    (3, 115, 30, False),
    (4, 249, 38, False),
    (5, 699, 46, False),
    (6, 2_003, 54, False),
    (7, 6_199, 62, False),
    (8, 19_999, 70, False),
)

VERIFICATION_FILES = {
    "engine.py",
    "A1_verify_tail.py",
    "A1_verify_tail_independent.py",
    "A1_reduce_check.py",
    "A8_conn2.py",
    "A8_l_curvature.py",
    "A10_linear_extremal.py",
    "A12_ratio_sign.py",
    "A13_delta5.py",
    "A14_uniform_sc.py",
    "A31_symbolic_strip.py",
    "A34_top_edge.py",
    "A35_t2_decrement_sign.py",
    "A36_bw_bridge.py",
    "A37_targeted_prefix.py",
    "A45_t1_stripe.py",
    "A38_targeted_core.py",
    "A2_stable.py",
    "A15_top_coverage.py",
    "A15_TOP_COVERAGE.md",
    "A39_global_infimum_retarget_probe.py",
    "A40_global_infimum_reduction.py",
    "A41_global_infimum_windows.py",
    "A42_global_infimum_join.py",
    "A43_core_uniqueness.py",
    "A44_dyadic_screen.py",
}

# The sharp-constant strengthening.  Theorem A (rho>=1, the factor 2) is proved
# by the seven fixed tails above and is left untouched.  Theorem B (rho >=
# rho(129,2), attained only at (129,2)) is proved by the components below and
# implies Theorem A; both are replayed independently.
SHARP_RESULT_FILES = {
    "global_infimum_tails.json",
    "global_infimum_reduction.json",
    "global_infimum_windows.json",
    "core_uniqueness.json",
    "dyadic_screen.json",
    "global_infimum_join.json",
}

RESULT_FILES = {
    *(f"fixed_tail_t{t}.json" for t, _, _, _ in TAILS),
    "finite_core.json",
    "fixed_prefixes.json",
    "symbolic_join.json",
    "t2_orientation.json",
    "terminal_cells.json",
    *SHARP_RESULT_FILES,
}

EXPECTED_MANIFEST_PATHS = tuple(
    sorted(
        {
            ".gitignore",
            "CITATION.cff",
            "ENVIRONMENT.md",
            "LICENSE",
            "LICENSE-CC-BY-4.0",
            "OPEN.md",
            "README.md",
            "RESULTS.md",
            "VERIFICATION_SCOPE.md",
            "paper/PANEL3_odd_cycles_sech_square_tail.pdf",
            "paper/PANEL3_odd_cycles_sech_square_tail.tex",
            "requirements.txt",
            "verify_all.py",
            *(f"verification/{name}" for name in VERIFICATION_FILES),
            *(f"results/{name}" for name in RESULT_FILES),
        }
    )
)

CERTIFICATE_SCHEMAS = {
    **{
        f"results/fixed_tail_t{t}.json": "square-tail-fixed-column-v1"
        for t, _, _, _ in TAILS
    },
    "results/finite_core.json": "square-tail-finite-core-v1",
    "results/fixed_prefixes.json": "square-tail-fixed-prefix-v1",
    "results/symbolic_join.json": "square-tail-symbolic-join-v2",
    "results/t2_orientation.json": "square-tail-t2-orientation-v1",
    "results/terminal_cells.json": "square-tail-terminal-cells-v1",
    "results/global_infimum_tails.json": "square-tail-global-infimum-retarget-probe-v3",
    "results/global_infimum_reduction.json": "square-tail-global-infimum-reduction-v1",
    "results/global_infimum_windows.json": "square-tail-global-infimum-windows-v1",
    "results/core_uniqueness.json": "square-tail-core-uniqueness-v1",
    "results/dyadic_screen.json": "square-tail-dyadic-screen-v1",
    "results/global_infimum_join.json": "square-tail-global-infimum-join-v1",
}

FORBIDDEN_PATH_MARKERS = tuple(
    bytes.fromhex(encoded)
    for encoded in ("2f686f6d652f", "2f55736572732f", "2e2e2f", "7e2f")
)
TEXT_SUFFIXES = {"", ".cff", ".gitignore", ".json", ".md", ".py", ".tex", ".txt"}
MANIFEST_LINE = re.compile(r"([0-9a-f]{64})  ([!-~]+)")
ACTIVE_PROCESSES: dict[int, subprocess.Popen] = {}
ACTIVE_PROCESSES_LOCK = threading.Lock()


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def inside_checkout(path: Path) -> bool:
    try:
        path.resolve().relative_to(HERE)
    except ValueError:
        return False
    return True


def ignored_inventory_path(relative: Path) -> bool:
    if not relative.parts:
        return False
    return relative.parts[0] == ".git"


def checkout_inventory() -> set[str]:
    inventory: set[str] = set()
    for path in HERE.rglob("*"):
        relative = path.relative_to(HERE)
        if ignored_inventory_path(relative):
            continue
        if path.is_symlink() or path.is_file():
            inventory.add(relative.as_posix())
    return inventory


def parse_manifest(text: str) -> tuple[dict[str, str], list[str]]:
    entries: dict[str, str] = {}
    failures: list[str] = []
    lines = text.splitlines()
    if not lines:
        return entries, ["manifest is empty"]
    if not text.endswith("\n"):
        failures.append("manifest lacks its final newline")
    for number, line in enumerate(lines, 1):
        match = MANIFEST_LINE.fullmatch(line)
        if match is None:
            failures.append(f"malformed manifest line {number}")
            continue
        digest, name = match.groups()
        pure = PurePosixPath(name)
        if pure.is_absolute() or ".." in pure.parts or name != pure.as_posix():
            failures.append(f"unsafe manifest path on line {number}: {name}")
            continue
        if name in entries:
            failures.append(f"duplicate manifest path: {name}")
            continue
        entries[name] = digest
    if list(entries) != sorted(entries):
        failures.append("manifest paths are not sorted")
    return entries, failures


def validate_manifest(text: str, check_inventory: bool = True) -> list[str]:
    entries, failures = parse_manifest(text)
    expected = set(EXPECTED_MANIFEST_PATHS)
    actual = set(entries)
    for missing in sorted(expected - actual):
        failures.append(f"manifest omits {missing}")
    for extra in sorted(actual - expected):
        failures.append(f"manifest contains unexpected path {extra}")
    for name in sorted(expected & actual):
        path = HERE / name
        if path.is_symlink():
            failures.append(f"release path is a symlink: {name}")
            continue
        if not path.is_file() or not inside_checkout(path):
            failures.append(f"missing or escaping release file: {name}")
            continue
        if sha256(path.read_bytes()) != entries[name]:
            failures.append(f"hash mismatch: {name}")
    if check_inventory:
        allowed = expected | {MANIFEST.name}
        for extra in sorted(checkout_inventory() - allowed):
            failures.append(f"unmanifested release file: {extra}")
    return failures


def all_gates_pass(gates) -> bool:
    if isinstance(gates, dict):
        return bool(gates) and all(value is True for value in gates.values())
    if isinstance(gates, list):
        return bool(gates) and all(
            isinstance(gate, dict) and gate.get("passed") is True for gate in gates
        )
    return gates is None


def validate_certificate(name: str, payload: bytes) -> list[str]:
    failures: list[str] = []
    try:
        record = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        return [f"certificate is not valid UTF-8 JSON: {error}"]
    if not isinstance(record, dict):
        return ["certificate root is not an object"]
    if record.get("schema") != CERTIFICATE_SCHEMAS[name]:
        failures.append(f"schema is not {CERTIFICATE_SCHEMAS[name]}")
    if record.get("passed") is not True:
        failures.append("certificate does not record passed=true")
    if record.get("status") not in (None, "PASS"):
        failures.append("certificate status is not PASS")
    if not all_gates_pass(record.get("gates")):
        failures.append("certificate contains a failed or malformed gate set")

    short = PurePosixPath(name).name
    if short.startswith("fixed_tail_t"):
        t = int(short.removeprefix("fixed_tail_t").removesuffix(".json"))
        spec = next(item for item in TAILS if item[0] == t)
        _, anchor, degree, controls = spec
        config = record.get("configuration", {})
        polynomial = record.get("certificate_polynomial", {})
        tail_gates = record.get("gates")
        if not isinstance(tail_gates, dict) or not tail_gates or not all(
            value is True for value in tail_gates.values()
        ):
            failures.append("fixed-tail gate set is missing, empty, or failed")
        if config.get("t") != t or config.get("anchor_m0") != anchor:
            failures.append("fixed-tail column or anchor mismatch")
        if config.get("controls") is not controls:
            failures.append("fixed-tail control policy mismatch")
        if polynomial.get("degree") != degree:
            failures.append("fixed-tail polynomial degree mismatch")
        if polynomial.get("coefficient_count") != degree + 1:
            failures.append("fixed-tail coefficient count mismatch")
    elif short == "finite_core.json":
        if record.get("cells") != 125_250 or record.get("blocks") != 20:
            failures.append("finite-core count mismatch")
        if record.get("algorithm") != "square-tail-finite-core-block-v2":
            failures.append("finite-core algorithm mismatch")
        if record.get("identity_replay_cells") != 120:
            failures.append("finite-core identity replay count mismatch")
        if record.get("malformed_blocks") != [] or record.get("unexpected_blocks") != []:
            failures.append("finite-core block inventory is not clean")
        if record.get("finite_box_witness", {}).get("cell") != [129, 2]:
            failures.append("finite-core witness cell mismatch")
    elif short == "fixed_prefixes.json":
        if record.get("cells") != 26_888:
            failures.append("fixed-prefix count mismatch")
        if record.get("algorithm") != "square-tail-fixed-prefix-row-v2":
            failures.append("fixed-prefix algorithm mismatch")
        if record.get("identity_replay_cells") != 84:
            failures.append("fixed-prefix identity replay count mismatch")
        if any(record.get(key) != [] for key in ("missing_rows", "unexpected_rows", "malformed_rows")):
            failures.append("fixed-prefix row inventory is not clean")
    elif short == "terminal_cells.json":
        if record.get("version") != 1 or record.get("failures") != []:
            failures.append("terminal certificate version or failure list mismatch")
        if record.get("symbolic", {}).get("status") != "PASS":
            failures.append("terminal symbolic certificate failed")
        replay = record.get("definition_replay", {})
        if replay.get("cells") != 21 or replay.get("status") != "PASS":
            failures.append("terminal replay count mismatch")
    elif short == "t2_orientation.json":
        if record.get("version") != 1:
            failures.append("t=2 orientation version mismatch")
        if record.get("cells") != 309 or record.get("sign_failures") != []:
            failures.append("t=2 orientation count or sign mismatch")
        if record.get("normalization_replay", {}).get("failures") != []:
            failures.append("t=2 normalization replay failed")
        tail = record.get("tail_orientation", {})
        if (
            tail.get("anchor_m0") != 312
            or tail.get("range") != "integer r>=313"
            or tail.get("passed") is not True
            or tail.get("gate")
            != "V*G-corr has nonnegative coefficients and positive constant term"
        ):
            failures.append("t=2 oriented tail is missing or failed")
    elif short == "symbolic_join.json":
        if record.get("version") != 2 or record.get("failures") != []:
            failures.append("symbolic-join version or failure list mismatch")
        if record.get("finite_sq_cells") != 152_138:
            failures.append("symbolic-join finite count mismatch")
        if record.get("tail_coefficients") != 329:
            failures.append("symbolic-join tail coefficient mismatch")
        if record.get("columns") != list(range(2, 9)):
            failures.append("symbolic join does not bind all seven columns")
        if not isinstance(record.get("gates"), dict) or not record["gates"]:
            failures.append("symbolic-join gate set is missing or empty")
        if record.get("mutation_controls") != {"boundary_t8_rejected": True}:
            failures.append("symbolic-join boundary control mismatch")
    elif short == "global_infimum_tails.json":
        columns = record.get("columns", {})
        if sorted(columns) != ["2", "3", "4", "5"]:
            failures.append("retargeted tails do not cover exactly t=2..5")
        expected_anchor = {"2": 503, "3": 503, "4": 503, "5": 1_000}
        for key, anchor in expected_anchor.items():
            entry = columns.get(key, {})
            if entry.get("anchor_m0") != anchor:
                failures.append(f"retargeted tail t={key} anchor mismatch")
            if entry.get("target") != "203/200":
                failures.append(f"retargeted tail t={key} target mismatch")
            if entry.get("negative_coefficient_indices") != []:
                failures.append(f"retargeted tail t={key} has a negative coefficient")
            if entry.get("coefficient_capacity_vs_target") != "above":
                failures.append(f"retargeted tail t={key} capacity is not above target")
        if record.get("controls", {}).get("original_t2_anchor_rejected") is not True:
            failures.append("the superseded t=2 anchor was not rejected")
    elif short == "global_infimum_reduction.json":
        if record.get("joint_row_free_core", {}).get("constant") != "513/500":
            failures.append("row-free core constant mismatch")
        crossings = record.get("level_crossings", {}).get("crossings", {})
        expected_crossings = {"25": 1_468, "30": 3_634, "35": 8_414, "36": 9_884, "40": 18_457}
        for level, row in expected_crossings.items():
            entry = crossings.get(level, {})
            if entry.get("first_r_with_a_ge_threshold") != row:
                failures.append(f"a_r crossing at {level} mismatch")
        if record.get("joint_floor_with_row_corrections", {}).get("row") != 504:
            failures.append("joint floor is not anchored at the worst row 504")
    elif short == "global_infimum_windows.json":
        if record.get("cells") != 19_146 or record.get("rows") != 16_483:
            failures.append("finite-window count mismatch")
        if record.get("identity_replay_cells") != 84:
            failures.append("finite-window identity replay count mismatch")
        if any(record.get(key) != [] for key in ("missing_rows", "malformed_rows")):
            failures.append("finite-window row inventory is not clean")
        ranges = {(item["t"], item["r_lo"], item["r_hi"]) for item in record.get("ranges", [])}
        if ranges != {(5, 504, 1_000), (6, 504, 3_633), (7, 1_468, 8_413), (8, 9_884, 18_456)}:
            failures.append("finite-window ranges are not the four residual windows")
    elif short == "core_uniqueness.json":
        if record.get("cells") != 125_250:
            failures.append("core-uniqueness cell count mismatch")
        if record.get("cells_at_or_below_witness") != [[129, 2]]:
            failures.append("core uniqueness: more than one cell at or below the witness")
        if record.get("cells_equal_to_witness") != [[129, 2]]:
            failures.append("core uniqueness: equality set is not exactly the witness")
        if record.get("unique_minimiser_on_core") is not True:
            failures.append("core uniqueness verdict missing")
        if "unique_global_minimiser" in record:
            failures.append("core certificate must not assert global uniqueness")
        if record.get("runner_up", {}).get("cell") != [128, 2]:
            failures.append("core-uniqueness runner-up cell mismatch")
        if record.get("runner_up", {}).get("separation_positive") is not True:
            failures.append("core-uniqueness runner-up is not strictly above the witness")
    elif short == "dyadic_screen.json":
        if record.get("cells") != 125_250:
            failures.append("dyadic-screen cell count mismatch")
        if record.get("exact_fallbacks") != 1:
            failures.append("dyadic screen must fall back exactly once")
        if record.get("cells_equal_to_witness") != [[129, 2]]:
            failures.append("dyadic screen: sole fallback is not the witness equality")
        if record.get("agreement_gate", {}).get("disagreements") != []:
            failures.append("dyadic screen disagreed with the exact comparison")
    elif short == "global_infimum_join.json":
        if not record.get("uniqueness", "").startswith("global"):
            failures.append("join does not assert global uniqueness")
        counts = record.get("counts", {})
        if counts.get("finite_core_cells") != 125_250:
            failures.append("join finite-core count mismatch")
        if counts.get("residual_window_cells") != 19_146:
            failures.append("join residual-window count mismatch")
        if counts.get("unbounded_fixed_column_tails") != 4:
            failures.append("join tail count mismatch")
        if sorted(record.get("inputs", {})) != sorted(
            ["core_uniqueness", "finite_core", "finite_windows",
             "retargeted_tails", "terminal_cells", "top_joint_reduction"]
        ):
            failures.append("join does not bind all six inputs")

    mutations = record.get("mutation_controls")
    if mutations is not None:
        if not isinstance(mutations, dict) or not mutations or not all(
            isinstance(value, bool) and value for value in mutations.values()
        ):
            failures.append("certificate contains a failed mutation control")
    return failures


def portable_release_sources() -> list[str]:
    failures: list[str] = []
    for name in EXPECTED_MANIFEST_PATHS:
        path = HERE / name
        if path.suffix not in TEXT_SUFFIXES or not path.is_file():
            continue
        content = path.read_bytes()
        for marker in FORBIDDEN_PATH_MARKERS:
            if marker in content:
                failures.append(
                    f"nonportable path marker {marker.decode('ascii')!r} in {name}"
                )
    return failures


def preflight() -> tuple[str, list[str]]:
    failures: list[str] = []
    if len(sys.argv) != 1:
        failures.append("verify_all.py has no quick, smoke, or skip mode")
    if sys.flags.optimize:
        failures.append("optimized Python is forbidden because assertions carry gates")
    try:
        sympy_version = version("sympy")
    except PackageNotFoundError:
        failures.append("required dependency is missing: sympy==1.14.0")
    else:
        if sympy_version != "1.14.0":
            failures.append(f"SymPy version is {sympy_version}, expected 1.14.0")
    if (
        MANIFEST.is_symlink()
        or not MANIFEST.is_file()
        or not inside_checkout(MANIFEST)
    ):
        return "", [
            *failures,
            "SHA256SUMS is missing, a symlink, or escapes the release tree",
        ]
    try:
        manifest_text = MANIFEST.read_text(encoding="ascii")
    except (OSError, UnicodeError) as error:
        return "", [*failures, f"cannot read SHA256SUMS: {error}"]
    failures.extend(validate_manifest(manifest_text))
    failures.extend(portable_release_sources())
    return manifest_text, failures


def mutated_manifest_is_rejected(manifest_text: str) -> bool:
    entries, failures = parse_manifest(manifest_text)
    if failures or not entries:
        return False
    first_name = next(iter(entries))
    old = entries[first_name]
    replacement = ("0" if old[0] != "0" else "1") + old[1:]
    mutated = manifest_text.replace(
        f"{old}  {first_name}", f"{replacement}  {first_name}", 1
    )
    return bool(validate_manifest(mutated, check_inventory=False))


def run_command(
    label: str,
    script_name: str,
    arguments: tuple[str, ...],
    timeout: int,
) -> tuple[bytes, str | None]:
    command = [sys.executable, "-I", "-B", str(VERIFY / script_name), *arguments]
    environment = {"LC_ALL": "C"}
    try:
        process = subprocess.Popen(
            command,
            cwd=VERIFY,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
    except OSError as error:
        return b"", f"could not launch {script_name}: {error}"
    with ACTIVE_PROCESSES_LOCK:
        ACTIVE_PROCESSES[process.pid] = process
    timed_out = False
    try:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = process.communicate()
    except BaseException:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.communicate()
        raise
    finally:
        with ACTIVE_PROCESSES_LOCK:
            ACTIVE_PROCESSES.pop(process.pid, None)
    transcript = stdout + stderr
    if timed_out:
        return transcript, f"{label} exceeded {timeout} seconds"
    if process.returncode != 0:
        return transcript, f"{script_name} exited {process.returncode}"
    return transcript, None


def terminate_active_processes() -> None:
    """Stop every component process group after a concurrent sibling fails."""
    with ACTIVE_PROCESSES_LOCK:
        processes = list(ACTIVE_PROCESSES.values())
    for process in processes:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def display_transcript(label: str, transcript: bytes) -> None:
    print(f"\n=== {label} ===", flush=True)
    if transcript:
        sys.stdout.buffer.write(transcript)
        if not transcript.endswith(b"\n"):
            print()
        sys.stdout.flush()


def compare_certificate(relative_name: str, generated_path: Path) -> list[str]:
    failures: list[str] = []
    try:
        generated = generated_path.read_bytes()
    except OSError as error:
        return [f"generated certificate is missing: {error}"]
    failures.extend(validate_certificate(relative_name, generated))
    try:
        banked = (HERE / relative_name).read_bytes()
    except OSError as error:
        failures.append(f"banked certificate is missing: {error}")
        return failures
    if generated != banked:
        failures.append("rebuilt certificate differs from the banked record")
    if not failures:
        record = json.loads(generated)
        record["passed"] = False
        mutated = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
        if not validate_certificate(relative_name, mutated):
            failures.append("passed=false certificate mutation was not rejected")
    return failures


def main() -> int:
    print(f"classification={CLASSIFICATION}")
    print(f"runtime_cpus={os.cpu_count() or 1}")
    manifest_text, failures = preflight()
    if failures:
        print("VERIFY_ALL: PREFLIGHT FAILED")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    if not mutated_manifest_is_rejected(manifest_text):
        print("VERIFY_ALL: PREFLIGHT FAILED")
        print("  - manifest mutation was not rejected")
        return 1

    for label, script_name, arguments in STATIC_COMPONENTS:
        transcript, failure = run_command(
            label, script_name, arguments, STATIC_TIMEOUT_SECONDS
        )
        display_transcript(f"{label}: {script_name}", transcript)
        if failure:
            print(f"VERIFY_ALL: {failure}")
            return 1

    rejected_certificate_mutations = 0
    with tempfile.TemporaryDirectory(prefix="square-tail-replay-") as temporary:
        temporary_root = Path(temporary)
        generated_results = temporary_root / "results"
        generated_results.mkdir()

        tail_jobs = {}
        tail_workers = min(len(TAILS), os.cpu_count() or 1)
        print(f"\n=== fixed tails: {len(TAILS)} single-worker jobs on {tail_workers} slots ===")
        pool = ThreadPoolExecutor(max_workers=tail_workers)
        aborting = False
        try:
            for t, anchor, _, controls in TAILS:
                output = generated_results / f"fixed_tail_t{t}.json"
                arguments = [str(t), str(anchor), "--json", str(output)]
                if controls:
                    arguments.append("--controls")
                future = pool.submit(
                    run_command,
                    f"fixed tail t={t}",
                    "A1_verify_tail_independent.py",
                    tuple(arguments),
                    GENERATOR_TIMEOUT_SECONDS,
                )
                tail_jobs[future] = (t, output)
            for future in as_completed(tail_jobs):
                t, output = tail_jobs[future]
                transcript, failure = future.result()
                display_transcript(f"fixed tail t={t}", transcript)
                if failure:
                    terminate_active_processes()
                    for pending in tail_jobs:
                        pending.cancel()
                    print(f"VERIFY_ALL: {failure}")
                    return 1
                name = f"results/fixed_tail_t{t}.json"
                certificate_failures = compare_certificate(name, output)
                if certificate_failures:
                    terminate_active_processes()
                    for pending in tail_jobs:
                        pending.cancel()
                    print(f"VERIFY_ALL: invalid {name}")
                    for item in certificate_failures:
                        print(f"  - {item}")
                    return 1
                rejected_certificate_mutations += 1
        except BaseException:
            aborting = True
            terminate_active_processes()
            for pending in tail_jobs:
                pending.cancel()
            raise
        finally:
            pool.shutdown(wait=not aborting, cancel_futures=True)

        generated_specs = (
            (
                "terminal cells",
                "A34_top_edge.py",
                "terminal_cells.json",
                ("--output", str(generated_results / "terminal_cells.json")),
            ),
            (
                "t=2 decrement orientation",
                "A35_t2_decrement_sign.py",
                "t2_orientation.json",
                ("--output", str(generated_results / "t2_orientation.json")),
            ),
            (
                "fixed-column prefixes",
                "A37_targeted_prefix.py",
                "fixed_prefixes.json",
                (
                    "--jsonl",
                    str(temporary_root / "fixed_prefixes.work.jsonl"),
                    "--output",
                    str(generated_results / "fixed_prefixes.json"),
                ),
            ),
            (
                "finite core",
                "A38_targeted_core.py",
                "finite_core.json",
                (
                    "--jsonl",
                    str(temporary_root / "finite_core.work.jsonl"),
                    "--output",
                    str(generated_results / "finite_core.json"),
                    "--workers",
                    str(os.cpu_count() or 1),
                    "--require-fresh",
                ),
            ),
        )
        for label, script_name, result_name, arguments in generated_specs:
            transcript, failure = run_command(
                label, script_name, arguments, GENERATOR_TIMEOUT_SECONDS
            )
            display_transcript(f"{label}: {script_name}", transcript)
            if failure:
                print(f"VERIFY_ALL: {failure}")
                return 1
            name = f"results/{result_name}"
            certificate_failures = compare_certificate(
                name, generated_results / result_name
            )
            if certificate_failures:
                print(f"VERIFY_ALL: invalid {name}")
                for item in certificate_failures:
                    print(f"  - {item}")
                return 1
            rejected_certificate_mutations += 1

        join_output = generated_results / "symbolic_join.json"
        join_arguments = (
            "--tails-dir",
            str(generated_results),
            "--prefix",
            str(generated_results / "fixed_prefixes.json"),
            "--core",
            str(generated_results / "finite_core.json"),
            "--terminal",
            str(generated_results / "terminal_cells.json"),
            "--orientation",
            str(generated_results / "t2_orientation.json"),
            "--output",
            str(join_output),
        )
        transcript, failure = run_command(
            "symbolic coverage join",
            "A31_symbolic_strip.py",
            join_arguments,
            GENERATOR_TIMEOUT_SECONDS,
        )
        display_transcript("symbolic coverage join: A31_symbolic_strip.py", transcript)
        if failure:
            print(f"VERIFY_ALL: {failure}")
            return 1
        certificate_failures = compare_certificate(
            "results/symbolic_join.json", join_output
        )
        if certificate_failures:
            print("VERIFY_ALL: invalid results/symbolic_join.json")
            for item in certificate_failures:
                print(f"  - {item}")
            return 1
        rejected_certificate_mutations += 1

        # ---- the sharp-constant strengthening (Theorem B) -----------------
        sharp_specs = (
            (
                "retargeted fixed tails t=2..5",
                "A39_global_infimum_retarget_probe.py",
                "global_infimum_tails.json",
                (
                    "--jsonl", str(temporary_root / "global_infimum_tails.work.jsonl"),
                    "--output", str(generated_results / "global_infimum_tails.json"),
                    "--workers", str(os.cpu_count() or 1),
                ),
            ),
            (
                "top/joint reduction and level crossings",
                "A40_global_infimum_reduction.py",
                "global_infimum_reduction.json",
                (
                    "--core", str(generated_results / "finite_core.json"),
                    "--output", str(generated_results / "global_infimum_reduction.json"),
                ),
            ),
            (
                "four bounded windows",
                "A41_global_infimum_windows.py",
                "global_infimum_windows.json",
                (
                    "--jsonl", str(temporary_root / "global_infimum_windows.work.jsonl"),
                    "--output", str(generated_results / "global_infimum_windows.json"),
                ),
            ),
            (
                "core uniqueness of (129,2)",
                "A43_core_uniqueness.py",
                "core_uniqueness.json",
                (
                    "--jsonl", str(temporary_root / "core_uniqueness.work.jsonl"),
                    "--output", str(generated_results / "core_uniqueness.json"),
                    "--workers", str(os.cpu_count() or 1),
                    "--require-fresh",
                ),
            ),
            (
                "dyadic screen accelerator",
                "A44_dyadic_screen.py",
                "dyadic_screen.json",
                (
                    "--jsonl", str(temporary_root / "dyadic_screen.work.jsonl"),
                    "--output", str(generated_results / "dyadic_screen.json"),
                    "--uniqueness", str(generated_results / "core_uniqueness.json"),
                    "--workers", str(os.cpu_count() or 1),
                    "--require-fresh",
                ),
            ),
            (
                "global infimum join",
                "A42_global_infimum_join.py",
                "global_infimum_join.json",
                (
                    "--retargeted-tails", str(generated_results / "global_infimum_tails.json"),
                    "--top-joint-reduction", str(generated_results / "global_infimum_reduction.json"),
                    "--finite-windows", str(generated_results / "global_infimum_windows.json"),
                    "--finite-core", str(generated_results / "finite_core.json"),
                    "--terminal-cells", str(generated_results / "terminal_cells.json"),
                    "--core-uniqueness", str(generated_results / "core_uniqueness.json"),
                    "--output", str(generated_results / "global_infimum_join.json"),
                ),
            ),
        )
        for label, script_name, result_name, arguments in sharp_specs:
            transcript, failure = run_command(
                label, script_name, arguments, GENERATOR_TIMEOUT_SECONDS
            )
            display_transcript(f"{label}: {script_name}", transcript)
            if failure:
                print(f"VERIFY_ALL: {failure}")
                return 1
            name = f"results/{result_name}"
            certificate_failures = compare_certificate(
                name, generated_results / result_name
            )
            if certificate_failures:
                print(f"VERIFY_ALL: invalid {name}")
                for item in certificate_failures:
                    print(f"  - {item}")
                return 1
            rejected_certificate_mutations += 1

    final_failures: list[str] = []
    if (
        MANIFEST.is_symlink()
        or not MANIFEST.is_file()
        or not inside_checkout(MANIFEST)
    ):
        final_failures.append(
            "SHA256SUMS became missing, a symlink, or external during replay"
        )
    else:
        try:
            final_manifest_text = MANIFEST.read_text(encoding="ascii")
        except (OSError, UnicodeError) as error:
            final_failures.append(f"cannot reread SHA256SUMS: {error}")
        else:
            if final_manifest_text != manifest_text:
                final_failures.append("SHA256SUMS changed during replay")
            final_failures.extend(validate_manifest(final_manifest_text))
    final_failures.extend(portable_release_sources())
    if final_failures:
        print("VERIFY_ALL: FINAL INVENTORY FAILED")
        for failure in final_failures:
            print(f"  - {failure}")
        return 1

    print(
        "\nVERIFY_ALL: PASS "
        f"({len(STATIC_COMPONENTS)} symbolic/transcription components, "
        f"{len(CERTIFICATE_SCHEMAS)} byte-stable certificates, "
        f"{1 + rejected_certificate_mutations} rejected release mutations; "
        f"{CLASSIFICATION})"
    )
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except KeyboardInterrupt:
        print("\nVERIFY_ALL: interrupted", file=sys.stderr)
        exit_code = 130
    finally:
        terminate_active_processes()
    raise SystemExit(exit_code)
