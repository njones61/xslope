"""Standing checkers for the verification pages under docs/verification.

Three checks run against each page:

  deltas   every printed percentage / absolute FS difference is re-derived from
           two numbers the page prints in the same sentence or table row;
  tags     every value a test tag locks is printed in the section carrying the
           tag, and every value the page presents as locked has a tag behind it;
  figures  every caption matches the figure it labels.

Change-gating lives in ``certify.py``: ``certified.json`` records the content
hash of each page as last certified, so an unchanged page costs one hash
compare and a changed page runs the full battery.  See README.md.
"""
