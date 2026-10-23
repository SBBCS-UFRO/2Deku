"""Ligando: estructura química con coordenadas 2D y mapa «número de serie PDB (PLIP) ↔ átomo RDKit».

Rutas de entrada, de la más a la menos fiable:

1. PDB del complejo + SMILES plantilla (``AssignBondOrdersFromTemplate``): los órdenes de enlace y las cargas formales
   salen del SMILES; la geometría y los números de serie, del PDB. El SMILES debe tener los mismos átomos pesados.
2. PDB del complejo + SDF del ligando (con enlaces): el SDF se asocia a los átomos del PDB por coordenadas.
3. Solo PDB: conectividad por proximidad y, si hay hidrógenos, ``rdDetermineBonds``; si falla, enlaces simples (aviso).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDepictor

try:  # RDKit ≥ 2022.09
    from rdkit.Chem import rdDetermineBonds
except Exception:  # pragma: no cover
    rdDetermineBonds = None


class LigandError(ValueError):
    """Error comprensible al construir el ligando."""


@dataclass
class Ligand:
    mol: Chem.Mol                           # solo átomos pesados, con conformero 2D
    serial_of: Dict[int, int]               # índice RDKit -> número de serie PDB
    warnings: List[str] = field(default_factory=list)
    name: str = ""

    @property
    def serial_to_idx(self) -> Dict[int, int]:
        return {v: k for k, v in self.serial_of.items()}


def _element(line: str) -> str:
    el = line[76:78].strip()
    if el:
        return el.capitalize()
    nm = line[12:16].strip()
    return "".join(ch for ch in nm if ch.isalpha())[:1].upper()


def read_pdb_ligand(pdb: str, resname: str, chain: Optional[str] = None, resnum: Optional[str] = None,
                    include_h: bool = False) -> Tuple[List[str], List[str]]:
    """Líneas PDB del ligando (primer residuo que cumpla) y avisos."""
    sel: Dict[Tuple[str, str], List[str]] = {}
    with open(pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:6] not in ("HETATM", "ATOM  "):
                continue
            if line[17:20].strip() != resname:
                continue
            ch, rn = line[21].strip(), line[22:26].strip()
            if chain and ch != chain:
                continue
            if resnum and rn != str(resnum):
                continue
            if not include_h and _element(line) in ("H", "D"):
                continue
            sel.setdefault((ch, rn), []).append(line)
    if not sel:
        raise LigandError(f"No hay átomos del residuo '{resname}' en {pdb}.")
    warns = []
    keys = list(sel)
    if len(keys) > 1:
        warns.append(f"El PDB contiene {len(keys)} copias de '{resname}' {keys}; se usa la primera {keys[0]}.")
    return sel[keys[0]], warns


def _finish(mol: Chem.Mol, serial_of: Dict[int, int], warns: List[str], name: str) -> Ligand:
    """2D con CoordGen sobre una copia sin conformero 3D."""
    m = Chem.Mol(mol)
    m.RemoveAllConformers()
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    return Ligand(m, serial_of, warns, name)


def from_pdb_smiles(pdb: str, resname: str, smiles: str, chain: Optional[str] = None,
                    resnum: Optional[str] = None) -> Ligand:
    lines, warns = read_pdb_ligand(pdb, resname, chain, resnum)
    RDLogger.DisableLog("rdApp.*")
    try:
        mol = Chem.MolFromPDBBlock("".join(lines), removeHs=False, sanitize=False, proximityBonding=True)
        tpl = Chem.MolFromSmiles(smiles)
        if tpl is None:
            raise LigandError(f"SMILES no válido: {smiles}")
        flat = Chem.Mol(tpl)
        Chem.RemoveStereochemistry(flat)
        if mol.GetNumAtoms() != flat.GetNumAtoms():
            raise LigandError(f"El PDB tiene {mol.GetNumAtoms()} átomos pesados de '{resname}' y el SMILES {flat.GetNumAtoms()}: "
                              "deben coincidir (mismos átomos pesados; los hidrógenos del SMILES son implícitos).")
        try:
            new = AllChem.AssignBondOrdersFromTemplate(flat, mol)
        except Exception as e:
            raise LigandError(f"El SMILES no encaja con la conectividad del PDB ({e}). Revise tautómero, anillos y "
                              "heteroátomos.") from e
        Chem.SanitizeMol(new)
    finally:
        RDLogger.EnableLog("rdApp.*")
    serial = {a.GetIdx(): a.GetPDBResidueInfo().GetSerialNumber() for a in new.GetAtoms()}
    return _finish(new, serial, warns, resname)


def _sdf_mol(sdf: str) -> Chem.Mol:
    if sdf.lower().endswith(".mol"):
        m = Chem.MolFromMolFile(sdf, removeHs=False)
    else:
        sup = Chem.SDMolSupplier(sdf, removeHs=False)
        m = next((x for x in sup if x is not None), None)
    if m is None:
        raise LigandError(f"No se pudo leer el ligando de {sdf}.")
    return Chem.RemoveHs(m)


def from_pdb_sdf(pdb: str, resname: str, sdf: str, chain: Optional[str] = None, resnum: Optional[str] = None,
                 tol: float = 0.8) -> Ligand:
    """SDF con enlaces; se asocia a los números de serie del PDB por elemento y coordenadas (o por orden)."""
    lines, warns = read_pdb_ligand(pdb, resname, chain, resnum)
    m = _sdf_mol(sdf)
    xyz = m.GetConformer().GetPositions()
    pdb_atoms = [(int(l[6:11]), _element(l), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])) for l in lines]
    if len(pdb_atoms) != m.GetNumAtoms():
        raise LigandError(f"El PDB tiene {len(pdb_atoms)} átomos pesados de '{resname}' y el SDF {m.GetNumAtoms()}.")
    serial, used = {}, set()
    for s, el, p in pdb_atoms:
        d = np.linalg.norm(xyz - p, axis=1)
        for j in np.argsort(d):
            if m.GetAtomWithIdx(int(j)).GetSymbol() == el and d[j] <= tol and int(j) not in used:
                serial[int(j)] = s; used.add(int(j)); break
    if len(serial) != len(pdb_atoms):          # el SDF no está en el mismo sistema de coordenadas: se prueba el orden
        elems_pdb = [el for _, el, _ in pdb_atoms]
        elems_sdf = [a.GetSymbol() for a in m.GetAtoms()]
        if elems_pdb == elems_sdf:
            serial = {i: s for i, (s, _, _) in enumerate(pdb_atoms)}
            warns.append("El SDF no coincide en coordenadas con el PDB: se asoció por orden de átomos.")
        else:
            raise LigandError("No se pudo asociar el SDF al PDB ni por coordenadas ni por orden de átomos.")
    return _finish(m, serial, warns, resname)


def from_pdb_only(pdb: str, resname: str, chain: Optional[str] = None, resnum: Optional[str] = None,
                  charge: int = 0) -> Ligand:
    """Sin SMILES ni SDF: se intenta deducir enlaces a partir de la geometría 3D (necesita hidrógenos)."""
    lines, warns = read_pdb_ligand(pdb, resname, chain, resnum, include_h=True)
    mol = Chem.MolFromPDBBlock("".join(lines), removeHs=False, sanitize=False, proximityBonding=True)
    has_h = any(a.GetSymbol() == "H" for a in mol.GetAtoms())
    ok = False
    if rdDetermineBonds is not None and has_h:
        try:
            m2 = Chem.Mol(mol)
            rdDetermineBonds.DetermineBondOrders(m2, charge=charge, allowChargedFragments=True)
            Chem.SanitizeMol(m2)
            mol, ok = m2, True
        except Exception:
            ok = False
    if not ok:
        warns.append("Sin SMILES/SDF no se pudieron deducir órdenes de enlace ni cargas: se dibujan enlaces simples. "
                     "Aporte --smiles o --sdf para un esquema correcto.")
        try:
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_KEKULIZE
                             ^ Chem.SANITIZE_SETAROMATICITY)
        except Exception:
            pass
    for a in mol.GetAtoms():
        a.SetIntProp("_serial", a.GetPDBResidueInfo().GetSerialNumber())
    heavy = Chem.RemoveHs(mol, sanitize=False)
    serial = {a.GetIdx(): a.GetIntProp("_serial") for a in heavy.GetAtoms()}
    return _finish(heavy, serial, warns, resname)


def build_ligand(pdb: str, resname: str, smiles: Optional[str] = None, sdf: Optional[str] = None,
                 chain: Optional[str] = None, resnum: Optional[str] = None, charge: int = 0) -> Ligand:
    """Elige la ruta según los datos disponibles (SMILES > SDF > solo PDB)."""
    if smiles:
        return from_pdb_smiles(pdb, resname, smiles, chain, resnum)
    if sdf:
        return from_pdb_sdf(pdb, resname, sdf, chain, resnum)
    return from_pdb_only(pdb, resname, chain, resnum, charge)
