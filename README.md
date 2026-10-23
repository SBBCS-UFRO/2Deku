# 2Deku

**Customizable 2D protein–ligand interaction diagrams from PLIP reports.**

2Deku turns a protein–ligand complex or an existing PLIP analysis into a readable
interaction diagram. It combines **PLIP** interaction data, **RDKit** ligand
coordinates and chemistry, and **Matplotlib** rendering to help prepare figures
for structural analysis, docking studies, and representative molecular dynamics
snapshots.

Use it from the command line, incorporate it into Python workflows, or access
its integrated viewer through **[CLIP! for UCSF ChimeraX](https://github.com/SBBCS-UFRO/CLIP)**.

The project is called **2Deku**; its Python package and command are named
`deku2d`. The standalone package currently declares version **0.1.0**. CLIP! has
its own version and installation process.

## Contents

- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Preparing your inputs](#preparing-your-inputs)
- [Command-line workflows](#command-line-workflows)
- [Styling and figure export](#styling-and-figure-export)
- [Comparing systems](#comparing-systems)
- [Python API](#python-api)
- [ChimeraX integration through CLIP!](#chimerax-integration-through-clip)
- [Interpretation and limitations](#interpretation-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Development and tests](#development-and-tests)
- [Scientific credit](#scientific-credit)

## Features

- Generate a single diagram, a panel with a shared legend, or a batch of figures.
- Read existing PLIP XML reports without running PLIP again.
- Run an external PLIP executable when starting from a PDB complex.
- Supply a SMILES template or an SDF/MOL file to define ligand bond orders and
  formal charges.
- Distinguish protein-donor and ligand-donor hydrogen bonds.
- Show residue identities, chain identifiers, interaction distances, and
  hydrophobic regions.
- Customize colors, fonts, line widths, residue boxes, legends, and layout search.
- Export **PNG**, **SVG**, and **PDF**, including multiple formats in one run.
- Compare interacting residue sets against a reference report.
- Use English or Spanish diagram legends.
- Render without a desktop display through Matplotlib's `Agg` backend.

### Interaction representation

| Interaction or element | Default appearance |
|---|---|
| Ligand | Black skeletal structure with labeled heteroatoms |
| Hydrophobic contacts | Yellow highlighting, orange arcs, and gray dotted connectors |
| Hydrogen bond: protein donor | Red dotted line with an arrow toward the ligand |
| Hydrogen bond: ligand donor | Green dotted line with an arrow toward the receptor residue |
| Salt bridge | Purple dotted line |
| π–π stacking | Cyan dotted line with a centroid marker |
| π–cation interaction | Orange dotted line with a centroid marker |
| Receptor residue | Labeled box, optionally tinted by residue category |
| Distance | Ångström label; enabled by default for hydrogen bonds and salt bridges |

**Five interaction classes are drawn.** Water bridges, halogen bonds, and metal
complexes are counted and reported as omitted when present in the PLIP report.
They are not converted into other interaction classes.

## Installation

### Requirements

| Component | Purpose |
|---|---|
| Python | Package metadata requires Python 3.9 or newer; Python 3.11 is used in the examples below |
| RDKit | Ligand chemistry, atom mapping, and 2D coordinates |
| NumPy | Geometry and layout calculations |
| Matplotlib | Diagram rendering and export |
| PLIP executable | Required only to calculate interactions from a PDB |

ChimeraX and PyMOL are not required for standalone use. Installing `deku2d` does
not install PLIP or configure its analysis dependencies.

### Install from the source folder

Open a terminal in the repository root: the folder containing `pyproject.toml`.
Create an isolated environment, then install the local package.

**macOS / Linux:**

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

**Windows PowerShell:**

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install .
```

The package installs NumPy, Matplotlib, and RDKit as dependencies. Their available
builds determine compatibility with your Python version and operating system.
For development, replace the final command with `python -m pip install -e .`.

Check the installation:

```bash
deku2d --version
deku2d --help
deku2d draw --help
```

You can also use `python -m deku2d` instead of `deku2d`. From the repository root,
this works without installing the local package if its dependencies are already
available. The optional `deku2d.sh` launcher uses `python3` by default; set
`DEKU2D_PYTHON` to select another interpreter.

### Configure PLIP only when needed

For PDB-based analysis, first prepare a working PLIP installation following the
[PLIP installation instructions](https://github.com/pharmai/plip#1-install-plip).
It may live in a separate environment.

Supply its executable explicitly:

```bash
deku2d draw --pdb complex.pdb --ligand-resname UNK \
  --sdf ligand.sdf --plip-exe /absolute/path/to/plip -o output/complex.png
```

Alternatively, set `DEKU2D_PLIP` before starting 2Deku:

```bash
export DEKU2D_PLIP="/absolute/path/to/plip"
```

On PowerShell, use `$env:DEKU2D_PLIP = "C:/path/to/plip.exe"`.
The current implementation checks for an executable file, so use its actual path
rather than assuming that a bare command name will be resolved from `PATH`.

## Quick start

The minimal example already includes a PLIP report and its matching corrected
PDB. **No PLIP installation is needed for this example.** Run it from the
repository root:

```bash
python -m deku2d draw \
  --plip-xml examples/mini/mol3_report.xml \
  --ligand-resname UNK \
  --smiles "OCc1[nH]c(O)cc1[C@@H]1C=CC[C@@H]2C[C@H](Sc3ncnc4[nH]cnc34)C[C@H]21" \
  --format png,svg,pdf \
  -o output/example.png
```

This creates `output/example.png`, `output/example.svg`, and
`output/example.pdf`. The output directory is created automatically.
The supplied example contains **13 supported interactions** involving **10
receptor residue boxes**.

On Windows, enter a multi-line Bash example as a single line, or replace the
continuation characters with the appropriate PowerShell syntax.

## Preparing your inputs

Each system needs a structure, a ligand identifier, and, preferably, a chemical
template for the ligand.

### Structure and interaction report

Choose either of these workflows:

| Input | Behavior |
|---|---|
| Complex PDB through `--pdb` | Runs PLIP and reads the resulting report |
| Existing XML through `--plip-xml` | Reads the report and looks for a neighboring `plipfixed*.pdb` |
| XML plus `--pdb` | Uses the supplied PDB if no neighboring corrected PDB exists |

Keep the report and its corresponding corrected PDB together. PLIP atom indices
refer to its processed structure and may differ from the original PDB serial
numbers. If multiple corrected PDB files share a report directory, the standalone
loader selects the first filename in sorted order; keep one matching result set
per directory.

The standalone pipeline can infer a constant serial-number offset using ligand
coordinates and reports when it does so. This is a recovery heuristic, not a
general remapping method for unrelated or reordered structures.

### Ligand selection

`--ligand-resname` is required. It identifies the ligand residue name, such as
`UNK`, `LIG`, or a PDB component code. Add `--chain` and `--resnum` when multiple
copies have the same name. Without those selectors, the first matching site or
residue may be used.

### Chemical template

| Input | What it provides | What to check |
|---|---|---|
| `--smiles` | Bond orders and formal charges assigned to the PDB ligand | Same heavy atoms and compatible connectivity, protonation, and tautomer |
| `--sdf` | Ligand chemistry from an SDF/MOL file | Same heavy atoms; preferably the same coordinates as the complex |
| Neither | Attempts bond-order inference from the PDB geometry | Requires suitable hydrogens and net charge; may fall back to single bonds |

Use one chemical template at a time. If both SMILES and SDF are supplied, the
standalone implementation gives SMILES priority. For automatic chemistry,
`--charge` sets the ligand's net charge and defaults to zero.

SMILES atom order does not need to match PDB atom order. Its heavy-atom count and
connectivity must match, and its formal charges should describe the analyzed
chemical state. SMILES stereochemistry is removed for template assignment; the
resulting flat drawing should not be used to communicate stereochemical identity.

SDF mapping first uses elements and coordinates. The standalone implementation
may fall back to atom order, with a warning, if coordinates do not match and the
element sequences agree. Review that warning before interpreting atom-specific
contacts. The CLIP! adapter rejects this atom-order fallback.

If bond orders cannot be inferred without a template, 2Deku warns that it is
showing single bonds. Treat that output as a **preview**, and provide a suitable
SMILES or SDF before preparing a scientific figure.

## Command-line workflows

### Draw one system

For your own existing PLIP result:

```bash
deku2d draw --plip-xml results/report.xml --ligand-resname UNK \
  --sdf ligand.sdf --chain A --resnum 900 \
  --title "Protein–ligand complex" --format png,svg,pdf -o output/complex.png
```

Replace the paths and ligand identifiers with those from your system. If a
corrected PDB is unavailable beside the report, add `--pdb matching_complex.pdb`
and verify that it corresponds to the report.

For a new PLIP calculation, replace `--plip-xml results/report.xml` with
`--pdb complex.pdb` and supply `--plip-exe`. The default command-line work directory
is based on the output name; override it with `--workdir`. A sufficiently recent
report and corrected PDB in that work directory may be reused. Choose a fresh
work directory when you need a new analysis after changing analysis settings.

### Create a panel

A panel combines systems into one figure with a shared legend:

```bash
deku2d panel --config examples/mini/panel.json --ncols 1 \
  --format png,svg,pdf -o output/panel.png
```

The supplied minimal configuration contains one system. Add entries to its
`systems` list to build a comparison panel. A two-system configuration has this
structure; replace the report paths and SMILES with your own data:

```json
{
  "ncols": 2,
  "label_style": "a)",
  "formats": ["png", "svg", "pdf"],
  "style": {"lang": "en", "dist": true, "tint": true, "dpi": 300},
  "systems": [
    {
      "name": "pose_a",
      "title": "Pose A",
      "plip_xml": "results/pose_a/report.xml",
      "ligand_resname": "UNK",
      "smiles": "CCO"
    },
    {
      "name": "pose_b",
      "title": "Pose B",
      "plip_xml": "results/pose_b/report.xml",
      "ligand_resname": "UNK",
      "smiles": "CCO"
    }
  ],
  "panel": {"output": "comparison", "ncols": 2}
}
```

`CCO` is an illustrative ethanol template, not the template for the bundled
example. Relative system input paths are resolved against the JSON file's
directory. Use distinct, filename-safe `name` values for batch outputs.

Useful system keys include `ligand_resname`, `pdb`, `plip_xml`, `smiles`, `sdf`,
`chain`, `resnum`, `charge`, `title`, `name`, `workdir`, `plip_exe`, `compare_xml`,
`compare_label`, `compare_resname`, and `as_salt`.

### Process a batch

With the configuration above saved as `systems.json` and its input files prepared:

```bash
deku2d batch --config systems.json --outdir output/batch
```

This saves a figure for each system. The optional top-level `panel` object also
requests a combined figure. `--format` overrides the configuration's `formats`.
For batch mode, a relative output directory is resolved against the configuration
file's directory, including when it is supplied through `--outdir`.

## Styling and figure export

### Frequently used options

| Option | Purpose | Default |
|---|---|---|
| `--lang en` / `--lang es` | Diagram legend language | `en` |
| `--fontsize` | Residue label font size | `17` |
| `--atom-fontsize` | Ligand atom font size | Residue font size + 2 in the default style |
| `--ppu` | Drawing scale in points per coordinate unit; lower values make text relatively larger | `34` |
| `--dpi` | Raster export resolution | `300` |
| `--min-width` | Minimum diagram width in drawing-coordinate units | `21` |
| `--dist` / `--no-dist` | Hydrogen-bond and salt-bridge distances | On |
| `--vdw-dist` / `--no-vdw-dist` | Hydrophobic contact distances | Off |
| `--tint` / `--no-tint` | Residue-category box colors | On |
| `--no-legend` | Hide the legend | Legend shown |
| `--seed` | Random seed for layout search | `3` |
| `--restarts` | Local-search restarts | `30` |
| `--iterations` | Moves per restart | `2500` |
| `--uniform` / `--no-uniform` | Uniform panel cells and centered ligands | On |
| `--mark-diff` / `--no-mark-diff` | Mark differences from a reference | On |
| `--legacy` | Disable the improved placement behavior | Off |
| `--style-json` | Load additional style parameters | None |
| `--format` | Comma-separated export formats | Output filename extension |
| `-q` / `--quiet` | Suppress the normal summary and scene warnings | Off |

Panel-specific options include `--ncols`, `--no-labels`, `--label-style`, and
`--no-titles`. See each subcommand's `--help` for the complete argument list.
The current CLI help and diagnostic messages are mainly in Spanish; `--lang`
changes the diagram legend, not the CLI language.

### Custom style file

Save this as `style.json` and pass `--style-json style.json` to a drawing command:

```json
{
  "fs": 18,
  "atom_fs": 20,
  "dpi": 300,
  "dist": true,
  "tint": true,
  "halo": "#FFE29A",
  "arc": "#E3A008",
  "red": "#D62728",
  "green": "#238B45",
  "bond_lw": 2.1,
  "box_lw": 1.4,
  "seed": 3
}
```

Style-file keys are the Python fields of `deku2d.config.Style`; for example,
`fs` corresponds to the CLI option `--fontsize`. Unknown keys produce an error.
The `texts` field can override individual legend strings.

Precedence is: default or legacy style → configuration `style` → `--style-json`
→ explicit CLI style options. PNG uses the requested DPI; SVG and PDF preserve
vector drawing elements. SVG text is exported as paths by default, and figure
backgrounds are white.

## Comparing systems

Supply `--compare-xml reference/report.xml` and optionally
`--compare-label "Reference pose"` to `draw`. Reference comparison marks residues
with supported interactions in the current diagram that are absent from the
reference interaction set, and can list reference-only residues. Residue name,
number, and chain identify the compared residues, so their identifiers must be
consistent across systems. This is a residue-set comparison, not an atom-by-atom
comparison of interaction geometries.

In a system configuration, use `compare_xml`, `compare_label`, and, if needed,
`compare_resname` when the reference ligand uses another residue name.

The advanced `--as-salt ASP12` option displays that residue's hydrogen bonds as
salt bridges. It can be repeated, or provided as an `as_salt` list in JSON.
This is a **manual visual reclassification**, explicitly reported as a warning;
it is not a new salt-bridge assignment made by PLIP. Disclose its use when
interpreting a figure.

## Python API

The same pipeline can be used without the command-line interface:

```python
from deku2d import Style, SystemSpec, scene_from_spec, render_single, save_figure

spec = SystemSpec(
    ligand_resname="UNK",
    plip_xml="examples/mini/mol3_report.xml",
    smiles="OCc1[nH]c(O)cc1[C@@H]1C=CC[C@@H]2C[C@H](Sc3ncnc4[nH]cnc34)C[C@H]21",
    name="example",
    title="Example complex",
)
style = Style(lang="en", dpi=300, seed=3)
scene = scene_from_spec(spec, style)
paths = save_figure(
    render_single(scene, style),
    "output/python_example.png",
    ["png", "svg", "pdf"],
    style.dpi,
)
print(paths)
print(scene.warnings)
```

Use `render_panel(scenes, style, ncols=2)` to combine scenes before saving.
`save_figure` creates the output directory and closes the figure after export.
Unlike JSON system paths, relative paths passed directly to `SystemSpec` are
relative to the Python process's working directory.

## ChimeraX integration through CLIP!

[CLIP!](https://github.com/SBBCS-UFRO/CLIP) integrates 2Deku into a native **2D view**
tab for UCSF ChimeraX. Its adapter connects diagrams to the selected PLIP result,
binding site, and interaction filters.

The integrated workflow provides automatic diagram generation after analysis
and on site/filter changes, a larger zoomable preview, collapsible manual ligand
inputs, display controls, and PNG/SVG/PDF export. It can attempt a compatible
SMILES from the report before falling back to the analyzed structure.

The standalone commands documented here do not install the ChimeraX plugin.
Use CLIP!'s installation guide for the bundle and its separate RDKit/Matplotlib
environment. CLIP! exposes a single-site viewer; standalone panels, batch jobs,
and reference-comparison options are not all exposed in that interface.
CLIP! offers seven interface languages, while 2Deku's built-in diagram legends
remain English and Spanish.

## Interpretation and limitations

2Deku arranges a diagram; PLIP supplies the interaction assignments and reported
measurements. The layout searches for residue-box positions around the ligand
while penalizing crossings, long connectors, and crowded labels. Some overlaps
or awkward connectors can remain. Try another seed or a larger layout search,
and review the exported figure.

For reproducible work, retain the structure, PLIP XML and corrected PDB,
chemical template, configuration, random seed, software versions, and warnings.
Matching seeds help repeat placement within the same software environment;
identical files are not guaranteed across dependency versions or platforms.

Important boundaries:

- Only the five supported interaction classes are drawn. Omitted classes may be
  important to the binding site.
- A PDB input should contain both receptor and ligand. An isolated ligand does
  not provide a protein–ligand interaction analysis.
- A docking pose or exported MD frame represents a **single static structure**.
  The program does not read trajectories or calculate contact occupancies.
- Diagram contacts and distances do not represent binding energies or affinity
  predictions.
- Bond orders, charges, protonation, and tautomer choice require review,
  particularly when no chemical template was supplied.
- Native PDB identifiers and serial-number mapping matter; arbitrarily renaming,
  renumbering, or changing coordinates can make an existing report incompatible.
- This README documents implemented workflows; it does not claim that every
  operating system or PLIP version has been validated.

## Troubleshooting

| Symptom | What to check |
|---|---|
| `deku2d` command not found | Activate the installation environment or use `python -m deku2d` |
| Missing RDKit, NumPy, or Matplotlib | Install into the same Python environment used to run the command |
| PLIP executable not found | Supply an actual executable path through `--plip-exe` or `DEKU2D_PLIP`, or use an existing XML report |
| Ligand or binding site not found | Check residue name, chain, and residue number in the processed structure and report |
| SMILES does not match | Check heavy-atom count, connectivity, protonation, and tautomer |
| SDF mapped by atom order | Prefer a coordinate-matched SDF or use a validated SMILES template |
| Report contacts cannot be mapped | Restore the matching corrected PDB and inspect serial-offset warnings |
| Single bonds shown unexpectedly | Provide a chemical template instead of relying on automatic inference |
| Crowded diagram | Try another `--seed`, more restarts/iterations, or revised font/scale settings |
| First render is slow or font-cache warning appears | Allow Matplotlib to build its cache; use a writable `MPLCONFIGDIR` if needed |

For issue reports, include the command, relevant configuration, package versions,
full error or warning, and a minimal input that you are able to share. Avoid
suppressing warnings with `--quiet` while diagnosing chemistry or mapping issues.

## Development and tests

Repository layout:

```text
deku2d/
  cli.py          Command-line interface
  config.py       Style fields and legend strings
  ligand.py       Chemical templates and ligand atom mapping
  plip_io.py      PLIP execution and XML parsing
  pipeline.py     Input-to-scene workflow
  layout.py       Diagram geometry and placement search
  render.py       Single-diagram rendering and export
  panel.py        Multi-system figures
tests/            Unit and CLI integration tests
examples/mini/    Minimal report, corrected PDB, complex, and panel configuration
examples/         Additional configurations and datasets
pyproject.toml    Package metadata, dependencies, and entry point
deku2d.sh         Optional Bash launcher
bin/deku2d        Optional Bash launcher
```

From the repository root, in an environment with the dependencies installed:

```bash
python -m unittest discover -s tests -v
```

The suite covers report parsing, ligand mapping, SMILES/SDF inputs, layout,
export formats, panels, CLI execution, and error handling. The PLIP execution
test is skipped when its configured executable is unavailable. Set
`DEKU2D_PLIP` before running the suite to configure that optional test.

## Scientific credit

2Deku relies on these projects:

- **[PLIP](https://github.com/pharmai/plip)** for interaction detection and reports.
  Use the citation recommended for the PLIP version used in your analysis.
- **[RDKit](https://www.rdkit.org/)** for ligand chemistry and 2D coordinates.
- **[Matplotlib](https://matplotlib.org/stable/users/project/citing.html)** for
  figure rendering and export.
- **NumPy** for geometry and numerical calculations.

Acknowledge 2Deku when using its diagrams and record the version or commit used.
The ChimeraX integration is provided through CLIP!; attribution to that integration
does not replace credit to the underlying analysis and rendering projects.
