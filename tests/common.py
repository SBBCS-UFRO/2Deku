"""Rutas y datos de las pruebas (ejemplo mínimo en ejemplos/mini)."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MINI = os.path.join(ROOT, "ejemplos", "mini")
XML = os.path.join(MINI, "mol3_report.xml")
FIXED = os.path.join(MINI, "plipfixed.complex_mol3_ehy6ig4f.pdb")
COMPLEX = os.path.join(MINI, "complex_mol3.pdb")
SMILES_MOL3 = "OCc1[nH]c(O)cc1[C@@H]1C=CC[C@@H]2C[C@H](Sc3ncnc4[nH]cnc34)C[C@H]21"
DATOS = os.path.join(ROOT, "ejemplos", "datos")
