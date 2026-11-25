"""One module per tutorial page, walked by ``test/tutorial_walk_check.py``.

Each module names its page (``PAGE``) and defines ``walk(w)``, which does what the
page tells the reader to do, in the page's order, on copies of the page's files:
Studio steps through the real window offscreen, script steps through the engine.
``w`` is a :class:`tutorial_walks.harness.Walk`; the shared machinery lives there.
"""
