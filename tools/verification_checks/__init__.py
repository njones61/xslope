"""Standing checkers for the verification pages under docs/verification.

Three checks run against each page:

  deltas   every printed percentage / absolute FS difference is re-derived from
           two numbers the page prints in the same sentence or table row;
  tags     every value a test tag locks is printed in the section carrying the
           tag, and every value the page presents as locked has a tag behind it;
  figures  every caption matches the figure it labels.

``tutorials.py`` runs beside them on the tutorial pages under docs/tutorials,
which carry no configs and no manifest: every factor of safety a tutorial
attributes to a run of its own is scored against the tags in scope — the page's
own, and the tags anywhere under docs/ on a model file the page links.  It
reports a per-page tally and never fails on a finding.

Change-gating lives in ``certify.py``: ``certified.json`` records the content
hash of each page as last certified, so an unchanged page costs one hash
compare and a changed page runs the full battery.  See README.md.
"""
