#!/usr/bin/env bash
# Make a branch working copy that tests like the main checkout.
#
#   tools/new_worktree.sh <branch> [base]      (base defaults to main)
#
# Creates the git worktree in a folder beside the main checkout
# (<main>_worktrees/<branch>, or $XSLOPE_WORKTREES/<branch>) and builds the
# optional fast FEM kernel from the branch's own source, which a new worktree
# lacks (it is git-ignored), so its strength-reduction runs would otherwise use
# the slower reference path. The suite's row timings and the private test repo
# need nothing: run_tests.py finds both through the main checkout.
#
# Set PYTHON to choose the interpreter (default python3); it needs Cython and a
# C compiler for the kernel, and the worktree is still usable without them.
set -euo pipefail

branch="${1:?usage: tools/new_worktree.sh <branch> [base]}"
base="${2:-main}"
here="$(cd "$(dirname "$0")/.." && pwd)"
common="$(git -C "$here" rev-parse --path-format=absolute --git-common-dir)"
main="$(dirname "$common")"
dest="${XSLOPE_WORKTREES:-${main}_worktrees}/$branch"

if git -C "$main" show-ref --verify --quiet "refs/heads/$branch"; then
    git -C "$main" worktree add "$dest" "$branch"
else
    git -C "$main" worktree add "$dest" -b "$branch" "$base"
fi

if (cd "$dest" && "${PYTHON:-python3}" setup_kernel.py build_ext --inplace >/dev/null 2>&1); then
    echo "fast kernel built in $dest/xslope"
else
    echo "fast kernel NOT built (needs Cython and a C compiler); FEM runs use the reference path" >&2
fi
rm -rf "$dest/build"
echo "worktree ready: $dest"
