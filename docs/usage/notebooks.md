# Colab Notebooks

XSLOPE provides a set of Google Colab notebooks that allow you to run analyses directly in your browser without installing anything locally. Each notebook is pre-configured to install xslope and provides an interactive workflow for uploading input files and viewing results.

## Seepage Analysis

Solve 2D seepage problems using the finite element method. Supports importing mesh data, solving for pore pressures, and integrating seepage results with slope stability analysis.

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_seep.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

## Limit Equilibrium Method

Solve slope stability problems using the limit equilibrium method. Supports multiple solution methods (Bishop, Spencer, Janbu, etc.), circular and non-circular failure surfaces, seismic loading, and automated search for critical surfaces.

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_lem.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

## Slope Design (LEM)

Find the critical slope angle that produces a target factor of safety. Sweeps a range of slope angles, runs automated circular search for each, and identifies the optimal design geometry.

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_design.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

## Finite Element Method

Solve slope stability problems using the finite element method with the Shear Strength Reduction Method (SSRM). Uses viscoplastic stress redistribution to find the critical factor of safety.

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_fem.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
