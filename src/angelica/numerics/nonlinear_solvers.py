from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class NonlinearProblem:
    residual_fn: Callable[[float], float]
    fixed_point_fn: Callable[[float], float] | None = None
    derivative_fn: Callable[[float], float] | None = None


class NonlinearSolver(ABC):
    def __init__(self, max_iterations: int = 50) -> None:
        self.max_iterations = max_iterations

    @abstractmethod
    def solve(
        self,
        problem: NonlinearProblem,
        initial_value: float,
        tolerance: float,
    ) -> float:
        raise NotImplementedError


class FixedPointSolver(NonlinearSolver):
    def solve(
        self,
        problem: NonlinearProblem,
        initial_value: float,
        tolerance: float,
    ) -> float:
        if problem.fixed_point_fn is None:
            raise ValueError("FixedPointSolver requires a fixed_point_fn.")

        value = initial_value
        for _ in range(self.max_iterations):
            updated_value = problem.fixed_point_fn(value)
            denominator = max(abs(updated_value), 1e-12)
            if abs(updated_value - value) / denominator <= tolerance:
                return updated_value
            value = updated_value
        return value


class SecantSolver(NonlinearSolver):
    def solve(
        self,
        problem: NonlinearProblem,
        initial_value: float,
        tolerance: float,
    ) -> float:
        value_previous = initial_value
        value_current = initial_value * 1.05 if abs(initial_value) >= 1e-12 else 1e-6
        residual_previous = problem.residual_fn(value_previous)
        residual_current = problem.residual_fn(value_current)

        for _ in range(self.max_iterations):
            denominator = residual_current - residual_previous
            if abs(denominator) < 1e-14:
                return value_current

            value_next = value_current - residual_current * (
                value_current - value_previous
            ) / denominator
            residual_next = problem.residual_fn(value_next)

            value_scale = max(abs(value_next), 1e-12)
            if (
                abs(residual_next) <= tolerance
                or abs(value_next - value_current) / value_scale <= tolerance
            ):
                return value_next

            value_previous, residual_previous = value_current, residual_current
            value_current, residual_current = value_next, residual_next

        return value_current


class NewtonSolver(NonlinearSolver):
    def solve(
        self,
        problem: NonlinearProblem,
        initial_value: float,
        tolerance: float,
    ) -> float:
        if problem.derivative_fn is None:
            raise ValueError("NewtonSolver requires a derivative_fn.")

        value = initial_value
        for _ in range(self.max_iterations):
            residual = problem.residual_fn(value)
            if abs(residual) <= tolerance:
                return value
            derivative = problem.derivative_fn(value)
            if abs(derivative) < 1e-14:
                return value
            value = value - residual / derivative
        return value


def build_nonlinear_solver(method_name: str, max_iterations: int = 50) -> NonlinearSolver:
    if method_name == "fixed_point":
        return FixedPointSolver(max_iterations=max_iterations)
    if method_name == "secant":
        return SecantSolver(max_iterations=max_iterations)
    if method_name == "newton":
        return NewtonSolver(max_iterations=max_iterations)
    raise ValueError(f"Unsupported nonlinear method: {method_name}")
