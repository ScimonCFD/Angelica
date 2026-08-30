from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pytest

openpyxl = pytest.importorskip("openpyxl")
load_workbook = openpyxl.load_workbook

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.results import ComponentFlowResult, GlobalBalance, GlobalEnergyBalance, SolveResult
from angelica.io.reporting import export_solve_result_workbook


def _minimal_result(**kwargs) -> SolveResult:
    defaults = {
        "case_name": "Test case",
        "converged": True,
        "node_pressures_pa": {1: 101325.0, 2: 95000.0},
        "component_flows": [
            ComponentFlowResult(
                label="Pipe:1",
                mass_flow_kg_per_s=1.23,
                volumetric_flow_m3_per_h=4.56,
            )
        ],
        "laminar_history": [],
        "laminar_metrics": [],
        "turbulent_history": [],
        "turbulent_metrics": [],
    }
    defaults.update(kwargs)
    return SolveResult(**defaults)


class ReportingTests(unittest.TestCase):
    def test_export_solve_result_workbook_no_balance(self) -> None:
        result = _minimal_result()  # global_balance=None by default

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "report.xlsx"
            export_solve_result_workbook(result, str(output_path))

            self.assertTrue(output_path.exists())

            workbook = load_workbook(output_path)
            self.assertEqual(workbook.sheetnames, ["Pressures", "Flows"])

            pressures_sheet = workbook["Pressures"]
            flows_sheet = workbook["Flows"]

            self.assertEqual(
                [cell.value for cell in pressures_sheet[1]],
                ["Node", "Pressure (Pa)", "Pressure (kPa)"],
            )
            self.assertEqual(
                [cell.value for cell in pressures_sheet[2]],
                [1, 101325.0, 101.325],
            )
            self.assertEqual(
                [cell.value for cell in flows_sheet[1]],
                ["Component", "Mass flow (kg/s)", "Vol. flow (m^3/h)"],
            )
            self.assertEqual(
                [cell.value for cell in flows_sheet[2]],
                ["Pipe:1", 1.23, 4.56],
            )

    def test_export_solve_result_workbook_with_balance_sheet(self) -> None:
        gb = GlobalBalance(mass_inlet_kg_per_s=2.0, mass_outlet_kg_per_s=1.998)
        geb = GlobalEnergyBalance(
            enthalpy_in_kw=500.0,
            enthalpy_out_kw=480.0,
            heat_sources_kw=0.0,
            heat_wall_loss_kw=20.0,
        )
        result = _minimal_result(global_balance=gb, global_energy_balance=geb)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "report.xlsx"
            export_solve_result_workbook(result, str(output_path))

            workbook = load_workbook(output_path)
            self.assertEqual(workbook.sheetnames, ["Pressures", "Flows", "Balance"])

            balance_sheet = workbook["Balance"]
            rows = [[cell.value for cell in row] for row in balance_sheet.iter_rows()]
            labels = [row[0] for row in rows]
            self.assertIn("Mass flow in", labels)
            self.assertIn("Mass flow out", labels)
            self.assertIn("Mass balance error", labels)
            self.assertIn("Enthalpy in", labels)
            self.assertIn("Energy balance error (%)", labels)


if __name__ == "__main__":
    unittest.main()
