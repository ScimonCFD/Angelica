GRAVITY_M_PER_S2 = 9.807


def elevation_pressure_term(density: float, height_change_m: float) -> float:
    return density * GRAVITY_M_PER_S2 * height_change_m
