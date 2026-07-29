# Claude Code Skill

## Overview

XSLOPE includes a custom [Claude Code](https://claude.ai/code) skill that allows you to build input templates and
run analyses using natural language and technical diagrams. With this skill you can:

- **Provide a diagram or sketch** of a slope problem and have Claude build a fully populated XSLOPE input template
- **Run analyses** by saying things like "find the factor of safety using Spencer's method" or "perform a seepage analysis"
- **Iterate on designs** by modifying materials, geometry, or boundary conditions through conversation

The skill understands the complete XSLOPE input template format and all three analysis types (LEM, seepage, and FEM).

The same knowledge is available inside the desktop app.
[XSlope Studio](../../studio/index.md) embeds a chat
**[AI assistant](../../studio/assistant.md)** that builds on this skill but works
*inside* the app: instead of writing an `.xlsx`, it populates the live project — which
renders immediately on the canvas — and can run analyses and script parameter studies
against the `xslope` API. It also supports models beyond Claude (OpenAI, local Ollama,
DeepSeek, Z.ai/GLM). This file-first skill and the document-first Studio assistant
coexist as two front ends to the same engine.

## Quick Start

Open a terminal in your xslope project directory and launch Claude Code:

```bash
claude
```

Use the `/xslope` command to get started. For example, paste or attach a diagram and type:

```text
/xslope Build an input template from this diagram
```

Once the input file is created, ask for an analysis in plain language — no need to repeat `/xslope`:

```text
Find the factor of safety using Spencer's method
```

```text
Now run a seepage analysis
```

The `/xslope` command loads the skill with detailed template layout and API knowledge. After the first
invocation, that context stays in the conversation, so follow-up requests can be typed as plain English.

## How It Works

### Building an Input Template

When you provide a technical diagram of a slope, the skill:

1. **Extracts data** from the image — slope geometry, material properties, boundary conditions, water tables, loads,
   reinforcement, and piles
2. **Checks for missing information** — if anything required is not visible in the diagram (unit weights, strength
   parameters, units, etc.), it will ask you before proceeding
3. **Copies the blank template** (`input_template.xlsx`) and populates each sheet using openpyxl
4. **Validates** the file by loading it with `load_slope_data()` and plotting the inputs
5. **Shows a summary** of everything that was populated (materials table, geometry, BCs, etc.) along with a download
   link and the validation plot

### Running an Analysis

Once an input file exists, you can ask for any of the following:

- **LEM analysis**: "Find the factor of safety using Spencer's method" — runs an automated circular search and
  plots the critical failure surface
- **Seepage analysis**: "Perform a seepage analysis" — builds a finite element mesh, solves the seepage equations,
  and plots head contours with flowlines and the phreatic surface
- **FEM analysis**: "Run an FEM analysis" — builds a mesh, runs the Shear Strength Reduction Method (SSRM), and
  plots deformation and shear strain

Each analysis produces multiple plots (inputs, mesh/search results, and final solution) so you can verify
every stage. The skill uses the same code patterns as `main_lem.py`, `main_seep.py`, and `main_fem.py`.

## Installation

The skill is a single markdown file (`.claude/commands/xslope.md`) that Claude Code automatically discovers
and makes available as a `/xslope` slash command.

### Option 1: Clone the repo (recommended)

If you clone the xslope repository, everything is ready to go — the skill file, the blank input template,
and the xslope package are all included:

```bash
git clone https://github.com/njones61/xslope.git
cd xslope
pip install -e .
```

The only additional requirement is [Claude Code](https://claude.ai/code).

### Option 2: Add the skill to an existing project

If you already have xslope installed (`pip install xslope`) and want to add the skill to a different
project directory:

1. Download the skill file: [xslope.md](xslope.md)

2. Place it in your project's `.claude/commands/` directory:

```bash
mkdir -p .claude/commands
cp ~/Downloads/xslope.md .claude/commands/xslope.md
```

3. Copy the blank input template from the xslope repository. The skill expects it at
`docs/inputs/input_template.xlsx`.

## Usage

### Typical Workflow

A typical session looks like this. Launch Claude Code from your xslope project directory:

```bash
claude
```

**Step 1 — Build the input file.** Use `/xslope` with a diagram or problem description. You can paste an
image directly into the terminal (drag-and-drop or clipboard paste) or reference a file path:

```text
/xslope Build an input template from this diagram. [paste/attach image]
```

Claude will examine the diagram, ask about any missing data (units, material properties, etc.), and produce
a populated `.xlsx` file. It shows a summary of what was populated and a validation plot of the geometry.

**Step 2 — Run an analysis.** Once the input file exists, just ask in plain English. No need to type `/xslope`
again — the skill context is already loaded:

```text
Find the factor of safety using Spencer's method
```

Claude will run the analysis and display all intermediate plots: the input geometry, the search results
(for auto-search), and the final solution with the critical failure surface.

**Step 3 — Iterate.** Continue the conversation to refine the problem:

```text
The FS is too low. Add a 500 psf berm at the toe and re-run.
```

```text
Change the cohesion of the clay layer to 300 psf and find the new FS.
```

```text
Now run a seepage analysis on that same problem.
```

### Example Output

When building an input file, Claude provides a structured summary:

```text
Input template created: inputs/earth_dam.xlsx

Geometry:
  - 3 profile lines defining 3 material layers (shell, core, foundation)
  - Origin at left edge of foundation, domain extends from x=0 to x=200
  - Max depth: 0

Materials:
  #  Name          gamma    c   phi   u
  1  shell           125    0    34   none
  2  core            122  100    26   none
  3  foundation      123    0    24   none

Failure Surfaces:
  - 2 starting circles: one at toe, one at base of core

Input file saved to: inputs/earth_dam.xlsx
```

## Supported Analysis Types

| Analysis | Example Prompt | Methods |
|----------|---------------|---------|
| LEM - Circular search | "Find FS using Spencer's method" | OMS, Bishop, Janbu, Corps of Engineers, Lowe & Karafiath, Spencer |
| LEM - Non-circular search | "Find FS with a non-circular surface using Janbu" | Janbu, Corps of Engineers, Lowe & Karafiath, Spencer |
| LEM - Single surface | "Analyze just the first circle using Bishop" | All methods |
| LEM - Reliability | "Run a reliability analysis" | All methods |
| Seepage | "Perform a seepage analysis" | FE steady-state, saturated/unsaturated |
| FEM - SSRM | "Run an FEM analysis" | Shear Strength Reduction Method |

## Tips

- **Provide clear diagrams.** The better the diagram, the more accurate the input file. Include coordinate labels,
  dimensions, material property tables, and boundary condition annotations when possible.
- **Indicate water tables clearly.** If your slope has a water table, mark it on the diagram with the
  standard inverted triangle symbol (&#x25BD;). A dashed line alone is not sufficient — it may be
  interpreted as a material boundary or construction line. Label the line "WT", "GWT", or "Water Table"
  for clarity.
- **Specify units.** Always state whether you are working in English (ft, pcf, psf) or Metric (m, kN/m3, kPa) units.
- **Check the validation plot.** Always review the `plot_inputs()` output before running an analysis to catch
  geometry errors early.
- **Start with LEM.** For most problems, start with a limit equilibrium analysis before moving to FEM. LEM is
  faster and gives a good first estimate of the factor of safety.
- **Multiple starting circles.** For multi-layer slopes, define one starting circle at the toe and one at the base
  of each layer to ensure the search finds the global minimum FS.
