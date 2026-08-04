"""Angelica: open-source platform for pipe network simulation."""

__version__ = "1.3.6"

from .closures.convection_scheme import ConvectionScheme, HybridScheme, PowerLawScheme, UpwindScheme
from .core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from .core.components import Fitting, HeatSource, Pipe, PressureChanger, Pump
from .core.results import ComponentFlowResult, SolveResult
from .core.settings import SolverSettings
from .io.reporting import print_solve_result
from .properties.compressible_fluid import CompressibleFluid
from .properties.gas_correlations import lee_gonzalez_eakin_viscosity
from .properties.dead_oil import build_thermal_dead_oil
from .properties.eos import EquationOfState, IdealGasEOS, PengRobinsonEOS
from .properties.single_component import SingleComponentFluid
from .properties.thermal_fluid import ThermalFluid
from .solvers import (
    BaseSolver,
    CompressibleSolverSettings,
    NonIsothermalSolverSettings,
    SteadyCompressibleSolver,
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
)

__all__ = [
    "BaseSolver",
    "ComponentFlowResult",
    "CompressibleFluid",
    "CompressibleSolverSettings",
    "ConvectionScheme",
    "EquationOfState",
    "FlowBoundary",
    "Fitting",
    "HeatSource",
    "HybridScheme",
    "IdealGasEOS",
    "PengRobinsonEOS",
    "lee_gonzalez_eakin_viscosity",
    "NetworkCase",
    "PowerLawScheme",
    "NonIsothermalSolverSettings",
    "Pipe",
    "PressureBoundary",
    "Pump",
    "PressureChanger",
    "SolveResult",
    "SingleComponentFluid",
    "SolverSettings",
    "SteadyCompressibleSolver",
    "SteadyIsothermalIncompressibleSolver",
    "SteadyNonIsothermalIncompressibleSolver",
    "ThermalBoundary",
    "ThermalFluid",
    "UpwindScheme",
    "build_thermal_dead_oil",
    "print_solve_result",
]
