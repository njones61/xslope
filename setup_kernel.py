# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Local build helper for the optional compiled Mohr-Coulomb SSRM kernel.

The kernel (``xslope/_fem_kernel.pyx``) is an OPT-IN accelerator for the
viscoplastic Step-6 constitutive update in ``fem.solve_fem``. It is NOT built by
``pip install xslope`` and is not required to run xslope: the pure-NumPy path is
the permanent default and the oracle against which every locked factor of safety
is defined. Building this kernel only speeds up SSRM runs that pass
``fast_kernel=True``; results are bit-identical to the NumPy path (enforced by
the ``kernel_xcheck`` CI gate).

Requirements: Cython and a C compiler (plus NumPy headers, already a dependency).

    pip install Cython

Build in place (compiles ``xslope/_fem_kernel.<platform>.so`` next to the .pyx):

    python setup_kernel.py build_ext --inplace

Once built, ``fem.solve_fem(..., fast_kernel=True)`` (and ``solve_ssrm``) picks it
up automatically. If ``fast_kernel=True`` is requested without a built kernel,
solve_fem warns once and falls back to NumPy. The generated ``.c`` and ``.so``
are git-ignored build artifacts; only the ``.pyx`` source is tracked.

To remove the built artifacts:

    rm -f xslope/_fem_kernel.c xslope/_fem_kernel*.so && rm -rf build
"""

import numpy as np
from setuptools import setup, Extension

try:
    from Cython.Build import cythonize
except ImportError as exc:  # pragma: no cover - clear message when Cython absent
    raise SystemExit(
        "Cython is required to build the optional fast kernel.\n"
        "    pip install Cython\n"
        "then re-run:  python setup_kernel.py build_ext --inplace"
    ) from exc

extensions = [
    Extension(
        "xslope._fem_kernel",
        ["xslope/_fem_kernel.pyx"],
        include_dirs=[np.get_include()],
    )
]

setup(
    name="xslope-fem-kernel",
    description="Optional compiled Mohr-Coulomb SSRM kernel for xslope (opt-in).",
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
            "initializedcheck": False,
        },
    ),
    zip_safe=False,
)
