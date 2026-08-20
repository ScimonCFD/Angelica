from __future__ import annotations

import numpy as np
from scipy.sparse import issparse
from scipy.sparse.linalg import spsolve


def solve_linear_system(matrix, vector):
    if issparse(matrix):
        return spsolve(matrix, vector)
    return np.linalg.solve(matrix, vector)
