from __future__ import annotations

from abc import ABC, abstractmethod


class BaseSolver(ABC):
    @abstractmethod
    def solve(self, case):
        raise NotImplementedError
