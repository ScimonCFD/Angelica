from __future__ import annotations

import math

from angelica.numerics import NonlinearProblem, build_nonlinear_solver

from .gravity import GRAVITY_M_PER_S2, elevation_pressure_term
from .pressure_drop import PressureDropCorrelation


class ColebrookFrictionFactorStrategy:
    def solve(
        self,
        correlation: ColebrookPipeCorrelation,
        pipe_state,
        initial_guess: float,
        tolerance: float,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> float:
        raise NotImplementedError


class DirectColebrookFrictionFactorSolver(ColebrookFrictionFactorStrategy):
    def solve(
        self,
        correlation: ColebrookPipeCorrelation,
        pipe_state,
        initial_guess: float,
        tolerance: float,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> float:
        problem = NonlinearProblem(
            residual_fn=lambda friction_factor: correlation.evaluate_colebrook(
                pipe_state,
                friction_factor,
            ),
            derivative_fn=lambda friction_factor: correlation._colebrook_derivative(
                pipe_state,
                friction_factor,
            ),
        )
        solver = build_nonlinear_solver(
            friction_factor_method,
            max_iterations=friction_factor_max_iterations,
        )
        friction_factor = solver.solve(problem, max(initial_guess, 1e-8), tolerance)
        return max(friction_factor, 1e-8)


class TransformedColebrookFrictionFactorSolver(ColebrookFrictionFactorStrategy):
    def solve(
        self,
        correlation: ColebrookPipeCorrelation,
        pipe_state,
        initial_guess: float,
        tolerance: float,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> float:
        problem = NonlinearProblem(
            residual_fn=lambda transformed_friction: correlation.evaluate_colebrook(
                pipe_state,
                math.exp(transformed_friction),
            ),
            derivative_fn=lambda transformed_friction: correlation._colebrook_transformed_derivative(
                pipe_state,
                transformed_friction,
            ),
        )
        solver = build_nonlinear_solver(
            friction_factor_method,
            max_iterations=friction_factor_max_iterations,
        )
        transformed_friction = solver.solve(problem, math.log(max(initial_guess, 1e-8)), tolerance)
        return max(math.exp(transformed_friction), 1e-8)


def build_colebrook_friction_factor_strategy(
    strategy_name: str,
) -> ColebrookFrictionFactorStrategy:
    if strategy_name == "direct":
        return DirectColebrookFrictionFactorSolver()
    if strategy_name == "transformed":
        return TransformedColebrookFrictionFactorSolver()
    raise ValueError(f"Unsupported Colebrook friction-factor strategy: {strategy_name}")


class LaminarPipeCorrelation(PressureDropCorrelation):
    def calculate_velocity(
        self,
        pipe_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        driving_term = delta_p - elevation_pressure_term(density, pipe_state.component.height_change_m)
        return (
            pipe_state.component.diameter_m**2
            * driving_term
            / (32.0 * viscosity * pipe_state.component.length_m)
        )

    def calculate_coupling(self, pipe_state, density: float, viscosity: float) -> float:
        return (-density / (32.0 * viscosity)) * (
            pipe_state.area_m2 * pipe_state.component.diameter_m**2 / pipe_state.component.length_m
        )


class ColebrookPipeCorrelation(PressureDropCorrelation):
    def calculate_velocity(
        self,
        pipe_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
        colebrook_friction_strategy: str = "transformed",
        friction_factor_method: str = "newton",
        friction_factor_max_iterations: int = 50,
        velocity_loop_method: str = "fixed_point",
        velocity_loop_max_iterations: int = 50,
        velocity_loop_tolerance: float | None = None,
    ) -> float:
        if tolerance is None:
            raise ValueError("ColebrookPipeCorrelation requires a tolerance value.")

        v_tol = velocity_loop_tolerance if velocity_loop_tolerance is not None else tolerance

        initial_velocity = self._safe_velocity_guess(pipe_state.velocity_m_per_s)
        problem = NonlinearProblem(
            fixed_point_fn=lambda velocity: self._fixed_point_mapping(
                pipe_state,
                velocity,
                delta_p,
                density,
                viscosity,
                tolerance,
                colebrook_friction_strategy,
                friction_factor_method,
                friction_factor_max_iterations,
            ),
            residual_fn=lambda velocity: self._velocity_residual(
                pipe_state,
                velocity,
                delta_p,
                density,
                viscosity,
                tolerance,
                colebrook_friction_strategy,
                friction_factor_method,
                friction_factor_max_iterations,
            ),
        )
        solver = build_nonlinear_solver(
            velocity_loop_method,
            max_iterations=velocity_loop_max_iterations,
        )
        solved_velocity = solver.solve(problem, initial_velocity, v_tol)
        self._assign_pipe_state(
            pipe_state,
            solved_velocity,
            density,
            viscosity,
            tolerance,
            colebrook_friction_strategy,
            friction_factor_method,
            friction_factor_max_iterations,
        )
        return pipe_state.velocity_m_per_s

    def _fixed_point_mapping(
        self,
        pipe_state,
        velocity: float,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float,
        colebrook_friction_strategy: str,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> float:
        safe_velocity = self._safe_velocity_guess(velocity)
        friction_factor, _ = self._friction_factor_for_velocity(
            pipe_state,
            safe_velocity,
            density,
            viscosity,
            tolerance,
            colebrook_friction_strategy,
            friction_factor_method,
            friction_factor_max_iterations,
        )
        driving_term = delta_p - elevation_pressure_term(
            density,
            pipe_state.component.height_change_m,
        )
        return (
            2.0
            * pipe_state.component.diameter_m
            * driving_term
            / (
                density
                * friction_factor
                * pipe_state.component.length_m
                * max(abs(safe_velocity), 1e-12)
            )
        )

    def _velocity_residual(
        self,
        pipe_state,
        velocity: float,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float,
        colebrook_friction_strategy: str,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> float:
        safe_velocity = self._safe_velocity_guess(velocity)
        friction_factor, _ = self._friction_factor_for_velocity(
            pipe_state,
            safe_velocity,
            density,
            viscosity,
            tolerance,
            colebrook_friction_strategy,
            friction_factor_method,
            friction_factor_max_iterations,
        )

        driving_term = delta_p - elevation_pressure_term(
            density,
            pipe_state.component.height_change_m,
        )
        friction_term = (
            density
            * friction_factor
            * pipe_state.component.length_m
            * safe_velocity
            * abs(safe_velocity)
            / (2.0 * pipe_state.component.diameter_m)
        )
        return driving_term - friction_term

    def _friction_factor_for_velocity(
        self,
        pipe_state,
        velocity: float,
        density: float,
        viscosity: float,
        tolerance: float,
        colebrook_friction_strategy: str,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> tuple[float, float]:
        reynolds = density * abs(velocity) * pipe_state.component.diameter_m / viscosity
        pipe_state.reynolds = reynolds
        if reynolds < 2300.0:
            return 64.0 / max(reynolds, 1e-12), reynolds
        initial_guess = max(64.0 / max(reynolds, 1e-12), 1e-6)
        friction_factor = self.solve_colebrook(
            pipe_state,
            initial_guess,
            tolerance,
            colebrook_friction_strategy,
            friction_factor_method,
            friction_factor_max_iterations,
        )
        return friction_factor, reynolds

    def _assign_pipe_state(
        self,
        pipe_state,
        velocity: float,
        density: float,
        viscosity: float,
        tolerance: float,
        colebrook_friction_strategy: str,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> None:
        safe_velocity = self._safe_velocity_guess(velocity)
        friction_factor, reynolds = self._friction_factor_for_velocity(
            pipe_state,
            safe_velocity,
            density,
            viscosity,
            tolerance,
            colebrook_friction_strategy,
            friction_factor_method,
            friction_factor_max_iterations,
        )
        pipe_state.velocity_m_per_s = safe_velocity
        pipe_state.reynolds = reynolds
        pipe_state.friction_factor = friction_factor

    @staticmethod
    def _safe_velocity_guess(velocity: float) -> float:
        if abs(velocity) < 1e-6:
            return 1e-6 if velocity >= 0.0 else -1e-6
        return velocity

    def calculate_coupling(self, pipe_state, density: float, viscosity: float) -> float:
        return (
            -2.0
            * pipe_state.area_m2
            * pipe_state.component.diameter_m
            / (
                pipe_state.friction_factor
                * max(abs(pipe_state.velocity_m_per_s), 1e-12)
                * pipe_state.component.length_m
            )
        )

    def solve_colebrook(
        self,
        pipe_state,
        initial_guess: float,
        tolerance: float,
        colebrook_friction_strategy: str,
        friction_factor_method: str,
        friction_factor_max_iterations: int,
    ) -> float:
        strategy = build_colebrook_friction_factor_strategy(colebrook_friction_strategy)
        return strategy.solve(
            self,
            pipe_state,
            max(abs(initial_guess), 1e-8),
            tolerance,
            friction_factor_method,
            friction_factor_max_iterations,
        )

    def _colebrook_transformed_derivative(
        self,
        pipe_state,
        transformed_friction: float,
    ) -> float:
        f = max(math.exp(transformed_friction), 1e-8)
        re = max(pipe_state.reynolds, 1e-12)
        g = (pipe_state.component.absolute_roughness_m / pipe_state.component.diameter_m) / 3.7 + 2.51 / (re * math.sqrt(f))
        # dR/dt = dR/df * f; simplifies to -(1/(2√f) + 2.51/(ln10·Re·√f·g))
        return -(1.0 / (2.0 * math.sqrt(f))) - 2.51 / (math.log(10) * re * math.sqrt(f) * g)

    def _colebrook_derivative(self, pipe_state, friction_factor: float) -> float:
        f = max(friction_factor, 1e-8)
        re = max(pipe_state.reynolds, 1e-12)
        g = (pipe_state.component.absolute_roughness_m / pipe_state.component.diameter_m) / 3.7 + 2.51 / (re * math.sqrt(f))
        return -(1.0 / (2.0 * f ** 1.5)) - 2.51 / (math.log(10) * re * f ** 1.5 * g)

    def evaluate_colebrook(self, pipe_state, friction_factor: float) -> float:
        friction_factor = max(friction_factor, 1e-8)
        log_term = (
            (pipe_state.component.absolute_roughness_m / pipe_state.component.diameter_m) / 3.7
            + 2.51 / (max(pipe_state.reynolds, 1e-12) * math.sqrt(friction_factor))
        )
        return 1.0 / math.sqrt(friction_factor) + 2.0 * math.log10(log_term)


class HazenWilliamsPipeCorrelation(PressureDropCorrelation):
    HAZEN_WILLIAMS_EXPONENT = 1.852
    HAZEN_WILLIAMS_COEFFICIENT = 10.67
    HAZEN_WILLIAMS_DIAMETER_EXPONENT = 4.871

    def calculate_velocity(
        self,
        pipe_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
        colebrook_friction_strategy: str = "transformed",
        friction_factor_method: str = "newton",
        friction_factor_max_iterations: int = 50,
        velocity_loop_method: str = "fixed_point",
        velocity_loop_max_iterations: int = 50,
        velocity_loop_tolerance: float | None = None,
    ) -> float:
        del (
            viscosity,
            tolerance,
            colebrook_friction_strategy,
            friction_factor_method,
            friction_factor_max_iterations,
            velocity_loop_method,
            velocity_loop_max_iterations,
            velocity_loop_tolerance,
        )
        driving_term_pa = delta_p - elevation_pressure_term(
            density,
            pipe_state.component.height_change_m,
        )
        if abs(driving_term_pa) < 1e-12:
            return 0.0

        resistance = self._hazen_williams_resistance(pipe_state)
        headloss_m = abs(driving_term_pa) / (density * GRAVITY_M_PER_S2)
        volumetric_flow_m3_per_s = (
            headloss_m / resistance
        ) ** (1.0 / self.HAZEN_WILLIAMS_EXPONENT)
        volumetric_flow_m3_per_s *= 1.0 if driving_term_pa >= 0.0 else -1.0
        return volumetric_flow_m3_per_s / pipe_state.area_m2

    def calculate_coupling(self, pipe_state, density: float, viscosity: float) -> float:
        del viscosity
        resistance = self._hazen_williams_resistance(pipe_state)
        volumetric_flow_m3_per_s = max(
            abs(pipe_state.area_m2 * pipe_state.velocity_m_per_s),
            1e-12,
        )
        derivative = 1.0 / (
            GRAVITY_M_PER_S2
            * resistance
            * self.HAZEN_WILLIAMS_EXPONENT
            * volumetric_flow_m3_per_s ** (self.HAZEN_WILLIAMS_EXPONENT - 1.0)
        )
        return -derivative

    def _hazen_williams_resistance(self, pipe_state) -> float:
        c_factor = max(pipe_state.component.hazen_williams_c, 1e-12)
        return (
            self.HAZEN_WILLIAMS_COEFFICIENT
            * pipe_state.component.length_m
            / (
                c_factor ** self.HAZEN_WILLIAMS_EXPONENT
                * pipe_state.component.diameter_m ** self.HAZEN_WILLIAMS_DIAMETER_EXPONENT
            )
        )


class DarcyWeisbachModel:
    def __init__(self) -> None:
        self.laminar_correlation = LaminarPipeCorrelation()
        self.turbulent_correlation = ColebrookPipeCorrelation()

    def laminar_velocity(self, pipe_state, delta_p: float, density: float, viscosity: float) -> float:
        return self.laminar_correlation.calculate_velocity(pipe_state, delta_p, density, viscosity)

    def turbulent_velocity(self, pipe_state, delta_p: float, density: float, viscosity: float, tolerance: float) -> float:
        return self.turbulent_correlation.calculate_velocity(
            pipe_state,
            delta_p,
            density,
            viscosity,
            tolerance,
        )

    def laminar_coupling(self, pipe_state, density: float, viscosity: float) -> float:
        return self.laminar_correlation.calculate_coupling(pipe_state, density, viscosity)

    def turbulent_coupling(self, pipe_state) -> float:
        return self.turbulent_correlation.calculate_coupling(pipe_state, 0.0, 0.0)
