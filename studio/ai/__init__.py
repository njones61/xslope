"""AI assistant for XSlope Studio (Phase A spike).

A dockable chat that drives the app and the ``xslope`` engine via a single
``run_python`` tool — an in-process Python kernel with the live project document
and the engine preloaded. Claude-only for now (the official ``anthropic`` SDK);
multi-provider and richer tooling are later phases. See plan_studio.md §13.
"""
