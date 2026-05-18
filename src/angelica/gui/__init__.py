"""Minimal GUI tools for Angelica."""


def main() -> None:
    from .app import main as run_gui

    run_gui()


__all__ = ["main"]
