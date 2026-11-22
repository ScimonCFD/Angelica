"""
Regression tests for every tutorial: scripts and GUI files.

For each ``run.py``:
  - it runs without error from an unrelated working directory;
  - its output contains no nan/inf and no non-converged solve;
  - every printed number matches the stored reference to within one unit of
    its last printed digit (``tests/tutorial_reference/scripts``).

For each ``.gui.json``:
  - it solves, as the GUI does, to a converged and finite result;
  - the global mass balance closes (< 0.01 %);
  - the global energy balance closes (< 0.01 %): boundary fluxes, heat
    sources, wall loss (computed independently of the assembly), shaft work
    and the Joule-Thomson term of valves;
  - node pressures, temperatures and mass flows match the stored reference
    (``tests/tutorial_reference/gui``).

After an intentional change of results, regenerate the references with

    python tests/update_tutorial_references.py

and review the diff of ``tests/tutorial_reference`` before committing.
"""
from __future__ import annotations

import json
import math

import pytest
from _tutorial_golden import (
    GUI_FILES,
    SCRIPTS,
    compare_gui,
    compare_output,
    gui_ref_path,
    run_all,
    run_gui,
    run_script,
    sanity_problems,
    script_ref_path,
    tutorial_id,
)

ENERGY_TOL_PCT = 0.01
UPDATE_HINT = "If the change is intended, run: python tests/update_tutorial_references.py"


@pytest.fixture(scope="session")
def script_results(tmp_path_factory):
    return run_all(SCRIPTS, run_script, tmp_path_factory.mktemp("tutorial_scripts"))


@pytest.fixture(scope="session")
def gui_results(tmp_path_factory):
    return run_all(GUI_FILES, run_gui, tmp_path_factory.mktemp("tutorial_gui"))


def test_tutorials_discovered():
    assert len(SCRIPTS) >= 43
    assert len(GUI_FILES) >= 43


@pytest.mark.parametrize("run_py", SCRIPTS, ids=tutorial_id)
def test_tutorial_script(run_py, script_results):
    proc = script_results[run_py]
    assert proc.returncode == 0, f"exit code {proc.returncode}\n{proc.stderr[-3000:]}"
    problems = sanity_problems(proc.stdout)
    assert not problems, "; ".join(problems)

    ref = script_ref_path(run_py)
    assert ref.exists(), f"no reference output {ref.name}. {UPDATE_HINT}"
    diffs = compare_output(proc.stdout, ref.read_text(encoding="utf-8"))
    assert not diffs, "output differs from reference:\n  " + "\n  ".join(diffs) + f"\n{UPDATE_HINT}"


@pytest.mark.parametrize("gui_json", GUI_FILES, ids=tutorial_id)
def test_tutorial_gui_file(gui_json, gui_results):
    proc = gui_results[gui_json]
    assert proc.returncode == 0, f"exit code {proc.returncode}\n{proc.stderr[-3000:]}"
    res = json.loads(proc.stdout)

    assert res["converged"], "GUI case did not converge"
    values = [*res["node_pressures_pa"].values(), *res["node_temperatures_c"].values(),
              *res["mass_flows_kg_per_s"].values()]
    assert all(math.isfinite(v) for v in values), "non-finite result"
    if res["mass_error_pct"] is not None:
        assert res["mass_error_pct"] < 0.01, f"mass balance error {res['mass_error_pct']:.4f} %"
    if res["energy_error_pct"] is not None:
        assert res["energy_error_pct"] < ENERGY_TOL_PCT, f"energy balance error {res['energy_error_pct']:.3f} %"

    ref = gui_ref_path(gui_json)
    assert ref.exists(), f"no reference result {ref.name}. {UPDATE_HINT}"
    diffs = compare_gui(res, json.loads(ref.read_text(encoding="utf-8")))
    assert not diffs, "result differs from reference:\n  " + "\n  ".join(diffs) + f"\n{UPDATE_HINT}"
