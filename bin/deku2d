#!/usr/bin/env bash
# Lanzador de deku2d: se puede llamar desde cualquier carpeta (sigue enlaces simbólicos).
#   ./deku2d.sh draw --pdb complejo.pdb --ligand-resname UNK --smiles "..." -o fig.png
# Variable opcional: DEKU2D_PYTHON (intérprete con rdkit, matplotlib y numpy; por defecto python3).
SRC="${BASH_SOURCE[0]}"
while [ -L "$SRC" ]; do
  DIR="$(cd -P "$(dirname "$SRC")" && pwd)"; SRC="$(readlink "$SRC")"; [[ "$SRC" != /* ]] && SRC="$DIR/$SRC"
done
DIR="$(cd -P "$(dirname "$SRC")" && pwd)"
[ -f "$DIR/deku2d/__main__.py" ] || DIR="$(dirname "$DIR")"   # cuando se llama como bin/deku2d
PYTHONPATH="$DIR${PYTHONPATH:+:$PYTHONPATH}" exec "${DEKU2D_PYTHON:-python3}" -m deku2d "$@"
