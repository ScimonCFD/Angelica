"""
Helpers for the tutorial regression tests (``tests/test_tutorials.py``) and the
reference updater (``tests/update_tutorial_references.py``).

Two kinds of reference are stored under ``tests/tutorial_reference/``:

``scripts/<solver>__<tutorial>.txt``
    The full stdout of ``run.py``.  It is compared line by line: the text
    skeleton must match exactly, and every printed number must agree to within
    one unit of its last printed digit.  Any change to a result that is large
    enough to show up in the tutorial's own output is therefore caught.

``gui/<solver>__<tutorial>.json``
    Node pressures, node temperatures and component mass flows obtained by
    loading the ``.gui.json`` file and solving it exactly as the GUI does.

Run as a script (``python tests/_tutorial_golden.py <file.gui.json>``) it solves
one GUI file and prints the result as JSON; the tests call it in a subprocess
so the GUI cases can run in parallel.
"""
from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
TUTORIAL_ROOT = SRC_ROOT / "angelica" / "tutorials"
REF_ROOT = Path(__file__).resolve().parent / "tutorial_reference"

SCRIPTS = sorted(TUTORIAL_ROOT.glob("*/*/run.py"))
GUI_FILES = sorted(TUTORIAL_ROOT.glob("*/*/*.gui.json"))

TIMEOUT_S = 600
_WORKERS = max(2, min(8, os.cpu_count() or 2))


def tutorial_id(path: Path) -> str:
    """``steady_black_oil/07_hot_oil_pipeline/run.py`` → ``steady_black_oil__07_hot_oil_pipeline``."""
    return f"{path.parent.parent.name}__{path.parent.name}"


def script_ref_path(run_py: Path) -> Path:
    return REF_ROOT / "scripts" / f"{tutorial_id(run_py)}.txt"


def gui_ref_path(gui_json: Path) -> Path:
    return REF_ROOT / "gui" / f"{tutorial_id(gui_json)}.json"


def _env() -> dict:
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONHASHSEED"] = "0"
    env["MPLBACKEND"] = "Agg"
    return env


# ── running ──────────────────────────────────────────────────────────────────

def run_script(run_py: Path, cwd: Path) -> subprocess.CompletedProcess:
    """Run a tutorial script exactly as a user would, from an unrelated directory."""
    return subprocess.run(
        [sys.executable, "-W", "ignore", str(run_py)],
        cwd=cwd, capture_output=True, text=True, encoding="utf-8",
        timeout=TIMEOUT_S, env=_env(),
    )


def run_gui(gui_json: Path, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-W", "ignore", str(Path(__file__).resolve()), str(gui_json)],
        cwd=cwd, capture_output=True, text=True, encoding="utf-8",
        timeout=TIMEOUT_S, env=_env(),
    )


def run_all(paths, runner, cwd: Path) -> dict:
    """Run every tutorial in parallel; returns {path: CompletedProcess}."""
    with ThreadPoolExecutor(max_workers=_WORKERS) as pool:
        futures = {p: pool.submit(runner, p, cwd) for p in paths}
        return {p: f.result() for p, f in futures.items()}


def solve_gui_file(gui_json: Path) -> dict:
    """Load a GUI scene, solve it as the GUI does, and summarise the result."""
    import warnings

    sys.path.insert(0, str(SRC_ROOT))
    warnings.filterwarnings("ignore")
    from angelica.gui.io import (
        build_network_case_from_scene,
        build_solver_from_scene,
        load_scene_from_file,
    )

    scene = load_scene_from_file(gui_json)
    result = build_solver_from_scene(scene).solve(build_network_case_from_scene(scene))
    gb = getattr(result, "global_balance", None)
    geb = getattr(result, "global_energy_balance", None)
    return {
        "converged": bool(result.converged),
        "node_pressures_pa": {str(k): v for k, v in sorted(result.node_pressures_pa.items())},
        "node_temperatures_c": {str(k): v for k, v in sorted(result.node_temperatures_c.items())},
        "mass_flows_kg_per_s": {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows},
        "mass_error_pct": None if gb is None else gb.mass_error_pct,
        "energy_error_pct": None if geb is None else geb.energy_error_pct,
    }


# ── comparing script output ──────────────────────────────────────────────────

_NUM = re.compile(r"(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _tolerance(token: str) -> float:
    """One unit in the last printed digit of ``token``."""
    mantissa, _, exponent = token.lower().partition("e")
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return 10.0 ** (int(exponent or 0) - decimals)


def _parse(text: str):
    out = []
    for line in text.splitlines():
        line = line.rstrip()
        out.append((_NUM.sub("#", line), _NUM.findall(line), line))
    return out


# Convergence residuals (pressure corrections, nodal mass imbalances) sit at
# round-off level and differ between CPUs; for these lines only require that
# the residual has not grown by more than a factor of 10.
_RESIDUAL_LINE = re.compile(r"pressure correction|imbalance|abs_dp|rel_dp", re.IGNORECASE)
_RESIDUAL_GROWTH = 10.0


def compare_output(actual: str, reference: str) -> list[str]:
    """Return a list of human-readable differences (empty when they agree)."""
    act, ref = _parse(actual), _parse(reference)
    problems: list[str] = []
    if len(act) != len(ref):
        problems.append(f"line count {len(act)} != reference {len(ref)}")
    for i, ((sa, na, la), (sr, nr, lr)) in enumerate(zip(act, ref), start=1):
        if sa != sr or len(na) != len(nr):
            problems.append(f"line {i} text changed:\n    ref: {lr}\n    now: {la}")
        else:
            residual = bool(_RESIDUAL_LINE.search(lr))
            for ta, tr in zip(na, nr):
                a, r = float(ta), float(tr)
                tol = 1.01 * max(_tolerance(tr), _tolerance(ta))
                if residual:
                    if abs(a) > _RESIDUAL_GROWTH * abs(r) + tol:
                        problems.append(f"line {i}: residual grew {tr} -> {ta}\n    ref: {lr}\n    now: {la}")
                        break
                elif abs(a - r) > tol:
                    problems.append(f"line {i}: {tr} -> {ta}\n    ref: {lr}\n    now: {la}")
                    break
        if len(problems) >= 8:
            problems.append("...")
            break
    return problems


_BAD_NUMBER = re.compile(r"(?<![A-Za-z])(nan|inf)(?![A-Za-z])", re.IGNORECASE)
_NOT_CONVERGED = re.compile(r"converged\s*[:=]\s*false", re.IGNORECASE)


def sanity_problems(stdout: str) -> list[str]:
    problems = []
    if _BAD_NUMBER.search(stdout):
        problems.append("output contains nan/inf")
    if _NOT_CONVERGED.search(stdout):
        problems.append("output reports a non-converged solve")
    return problems


# ── comparing GUI results ────────────────────────────────────────────────────

GUI_RTOL = 1e-6
GUI_ATOL = {"node_pressures_pa": 1e-3, "node_temperatures_c": 1e-6, "mass_flows_kg_per_s": 1e-9}


def compare_gui(actual: dict, reference: dict) -> list[str]:
    problems: list[str] = []
    if actual["converged"] != reference["converged"]:
        problems.append(f"converged {reference['converged']} -> {actual['converged']}")
    for key, atol in GUI_ATOL.items():
        a, r = actual[key], reference[key]
        if set(a) != set(r):
            problems.append(f"{key}: keys changed {sorted(set(r) ^ set(a))}")
            continue
        for k in r:
            if not math.isclose(a[k], r[k], rel_tol=GUI_RTOL, abs_tol=atol):
                problems.append(f"{key}[{k}]: {r[k]!r} -> {a[k]!r}")
    return problems[:12]


if __name__ == "__main__":
    print(json.dumps(solve_gui_file(Path(sys.argv[1])), indent=1, sort_keys=True))
