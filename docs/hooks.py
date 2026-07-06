import logging


def on_startup(command, dirty, **kwargs):
    """Silence griffe/mkdocstrings docstring-parsing warnings during build/serve."""
    logging.getLogger("mkdocs.plugins.griffe").setLevel(logging.ERROR)
