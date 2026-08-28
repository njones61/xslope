"""Scored scenarios for the Studio assistant.

``tools/assistant_suite.py`` is the runner; this package is what it runs and how
it scores. See :mod:`tools.assistant_scenarios.definitions` for the scenarios
themselves, :mod:`~tools.assistant_scenarios.scorers` for the criteria they are
graded on, and :mod:`~tools.assistant_scenarios.faults` for the models broken on
purpose that the diagnosis scenarios open.
"""

from .core import (Criterion, Fault, Scenario, ScoreCtx, Turn, digest,  # noqa: F401
                   load_model, lock, parse_body, repo, seep_flow, solve,
                   solve_variant, ssrm_fs, use_solve_cache)
from .definitions import SCENARIOS                                       # noqa: F401

__all__ = ["Criterion", "Fault", "Scenario", "ScoreCtx", "Turn", "SCENARIOS",
           "by_name", "families", "digest", "load_model", "lock", "parse_body",
           "repo", "seep_flow", "solve", "solve_variant", "ssrm_fs",
           "use_solve_cache"]


def by_name(name):
    """The scenario called ``name``, or ``None``."""
    return next((s for s in SCENARIOS if s.name == name), None)


def families():
    """The families, in the order they are defined."""
    seen = []
    for scenario in SCENARIOS:
        if scenario.family not in seen:
            seen.append(scenario.family)
    return seen
