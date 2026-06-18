"""Angelica: open-source platform for pipe network simulation."""

from .closures.convection_scheme import ConvectionScheme, HybridScheme, PowerLawScheme, UpwindScheme
from .core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from .core.components import Fitting, HeatSource, Pipe, PressureChanger
from .core.results import ComponentFlowResult, SolveResult
from .core.settings import SolverSettings
from .io.reporting import print_solve_result
from .properties.dead_oil import build_thermal_dead_oil
from .properties.single_component import SingleComponentFluid
from .properties.thermal_fluid import ThermalFluid
from .solvers import (
    BaseSolver,
    NonIsothermalSolverSettings,
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
)

__all__ = [
    "BaseSolver",
    "ComponentFlowResult",
    "ConvectionScheme",
    "FlowBoundary",
    "Fitting",
    "HeatSource",
    "HybridScheme",
    "NetworkCase",
    "PowerLawScheme",
    "NonIsothermalSolverSettings",
    "Pipe",
    "PressureBoundary",
    "PressureChanger",
    "SolveResult",
    "SingleComponentFluid",
    "SolverSettings",
    "SteadyIsothermalIncompressibleSolver",
    "SteadyNonIsothermalIncompressibleSolver",
    "ThermalBoundary",
    "ThermalFluid",
    "UpwindScheme",
    "build_thermal_dead_oil",
    "print_solve_result",
]
