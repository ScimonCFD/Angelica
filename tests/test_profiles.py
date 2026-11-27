"""
Per-connection profile report (roadmap item 7.4, v1.34.0).

Consistency checks (no external data: the report rearranges the solution):
along a connection split into segments, each segment starts where the
previous one ends (distance, elevation, node), the distances add up to the
pipe length and the elevations to its height change, and the pressure drops
add up to the drop across the connection.  Tolerance 1e-9: sums of the same
floating-point numbers.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from angelica.io.reporting import profile_rows, profile_table  # noqa: E402


@pytest.fixture(scope="module")
def segmented():
    """Black-oil tutorial 07's 20 km line as a GUI scene, split into 5
    segments with a 100 m rise."""
    from angelica.gui.io import build_network_case_from_scene, build_solver_from_scene, load_scene_from_file
    scene = load_scene_from_file(ROOT / "src/angelica/tutorials/steady_black_oil/07_hot_oil_pipeline/"
                                 "07_hot_oil_pipeline.gui.json")
    comp = scene.links[0].components[0]
    scene.update_link_component_properties(scene.links[0].link_id, comp.component_id,
                                           {**comp.properties, "num_segments": "5", "height_change_m": "100.0"})
    case = build_network_case_from_scene(scene)
    return case, build_solver_from_scene(scene).solve(case)


def test_segmented_connection_profile_is_continuous(segmented):
    case, res = segmented
    rows = [r for r in profile_table(case, res) if r["connection"] == "link_1"]
    assert len(rows) == 5
    for a, b in zip(rows, rows[1:]):
        assert b["from_node"] == a["to_node"]
        assert b["distance_start_m"] == a["distance_end_m"] and b["elevation_start_m"] == a["elevation_end_m"]
        assert b["pressure_in_pa"] == a["pressure_out_pa"]
    assert rows[-1]["distance_end_m"] == pytest.approx(20_000.0, rel=1e-9)
    assert rows[-1]["elevation_end_m"] == pytest.approx(100.0, rel=1e-9)
    total = sum(r["pressure_drop_pa"] for r in rows)
    assert total == pytest.approx(rows[0]["pressure_in_pa"] - rows[-1]["pressure_out_pa"], rel=1e-9)
    assert all(r["mass_flow_kg_per_s"] == pytest.approx(rows[0]["mass_flow_kg_per_s"], rel=2e-3) for r in rows)
    assert rows[-1]["temperature_out_c"] < rows[0]["temperature_in_c"]       # the hot oil cools along the line


def test_code_built_case_lists_each_component():
    from _two_loop_networks import build_case, solver_for
    case = build_case("gas")
    res = solver_for("gas").solve(case)
    table = profile_table(case, res)
    assert [r["component"] for r in table] == [c.component_id for c in case.components]
    assert all(r["connection"] == r["component"] for r in table)


def test_exports_carry_the_profiles(tmp_path, segmented):
    from angelica.io.reporting import export_solve_result_csv, export_solve_result_workbook
    case, res = segmented
    export_solve_result_csv(res, tmp_path / "r.csv", case=case)
    text = (tmp_path / "r.csv").read_text(encoding="utf-8")
    assert "Profiles" in text and "distance_end_m" in text
    export_solve_result_csv(res, tmp_path / "plain.csv")                     # without the case: as before
    assert "Profiles" not in (tmp_path / "plain.csv").read_text(encoding="utf-8")
    pytest.importorskip("openpyxl")
    export_solve_result_workbook(res, tmp_path / "r.xlsx", case=case)
    from openpyxl import load_workbook
    ws = load_workbook(tmp_path / "r.xlsx")["Profiles"]
    assert ws.max_row == len(profile_rows(case, res))


def test_many_connections_each_restart_and_chain():
    """Hanoi water network (34+ connections, link ids 1…, including link_1
    and link_10…): each connection's rows start at distance 0 and chain."""
    from angelica.gui.io import build_network_case_from_scene, build_solver_from_scene, load_scene_from_file
    scene = load_scene_from_file(ROOT / "src/angelica/tutorials/steady_isothermal_incompressible/"
                                 "07_hanoi_epanet_darcy_benchmark/hanoi_epanet_darcy_benchmark.gui.json")
    case = build_network_case_from_scene(scene)
    table = profile_table(case, build_solver_from_scene(scene).solve(case))
    by_conn: dict = {}
    for r in table:
        by_conn.setdefault(r["connection"], []).append(r)
    assert len(by_conn) == len(scene.links)
    for rows in by_conn.values():
        assert rows[0]["distance_start_m"] == 0.0
        assert all(b["from_node"] == a["to_node"] for a, b in zip(rows, rows[1:]))


def test_code_built_case_with_empty_ids_and_a_mismatched_result():
    from dataclasses import replace

    from _two_loop_networks import build_case, solver_for
    case = build_case("liquid")
    case = replace(case, components=tuple(replace(c, component_id="") for c in case.components))
    res = solver_for("liquid").solve(case)
    table = profile_table(case, res)
    assert len({r["connection"] for r in table}) == len(case.components)      # one connection each
    assert all(r["distance_start_m"] == 0.0 for r in table)
    assert all(r["mass_flow_kg_per_s"] == f.mass_flow_kg_per_s for r, f in zip(table, res.component_flows))
    smaller = replace(case, components=case.components[:-1])
    with pytest.raises(ValueError, match="not the solve"):
        profile_table(smaller, res)
