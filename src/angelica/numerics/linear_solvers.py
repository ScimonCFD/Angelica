from __future__ import annotations

import numpy as np


def solve_linear_system(matrix, vector):
    return np.linalg.solve(matrix, vector)
