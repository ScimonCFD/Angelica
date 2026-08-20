"""Angelica: open-source platform for pipe network simulation."""

from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("angelica")
except PackageNotFoundError:
    __version__ = "unknown"

from .closures.convection_scheme import ConvectionScheme, HybridScheme, PowerLawScheme, UpwindScheme
from .core.case import FlowBoundary, InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from .core.components import Fitting, HeatSource, Pipe, PressureChanger, Pump
from .core.results import ComponentFlowResult, SolveResult
from .core.settings import SolverSettings
from .io.reporting import print_solve_result
from .properties.black_oil import BlackOilFluid
from .properties.compressible_fluid import CompressibleFluid
from .properties.gas_correlations import lee_gonzalez_eakin_viscosity
from .properties.dead_oil import build_thermal_dead_oil
from .properties.eos import EquationOfState, IdealGasEOS, PengRobinsonEOS
from .properties.single_component import SingleComponentFluid
from .properties.thermal_fluid import ThermalFluid
from .properties.water_liquid import build_water_thermal_fluid
from .solvers import (
    BaseSolver,
    BlackOilSolverSettings,
    CompressibleSolverSettings,
    NonIsothermalSolverSettings,
    SteadyBlackOilSolver,
    SteadyCompressibleSolver,
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
)

__all__ = [
    "BaseSolver",
    "BlackOilFluid",
    "BlackOilSolverSettings",
    "ComponentFlowResult",
    "CompressibleFluid",
    "CompressibleSolverSettings",
    "ConvectionScheme",
    "EquationOfState",
    "Fitting",
    "FlowBoundary",
    "HeatSource",
    "HybridScheme",
    "IdealGasEOS",
    "InletFluidBC",
    "NetworkCase",
    "NonIsothermalSolverSettings",
    "PengRobinsonEOS",
    "Pipe",
    "PowerLawScheme",
    "PressureBoundary",
    "PressureChanger",
    "Pump",
    "SingleComponentFluid",
    "SolveResult",
    "SolverSettings",
    "SteadyBlackOilSolver",
    "SteadyCompressibleSolver",
    "SteadyIsothermalIncompressibleSolver",
    "SteadyNonIsothermalIncompressibleSolver",
    "ThermalBoundary",
    "ThermalFluid",
    "UpwindScheme",
    "build_thermal_dead_oil",
    "build_water_thermal_fluid",
    "lee_gonzalez_eakin_viscosity",
    "print_solve_result",
]
