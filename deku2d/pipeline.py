"""De los archivos de entrada a un ``Scene``: PLIP (ejecutar o leer), ligando 2D y diseño."""
from __future__ import annotations

import os
from dataclasses import dataclass, fields
from typing import Any, Dict, List, Optional

from .config import Style
from .layout import Scene, build_scene
from .ligand import build_ligand
from .plip_io import find_fixed_pdb, parse_report, run_plip


@dataclass
class SystemSpec:
    """Un sistema ligando–receptor a dibujar (rutas absolutas o ya resueltas)."""
    ligand_resname: str
    pdb: Optional[str] = None          # complejo (PDB). Con plip_xml se prefiere el 'plipfixed' que está junto al XML
    smiles: Optional[str] = None       # SMILES plantilla (mismos átomos pesados que el PDB)
    sdf: Optional[str] = None          # alternativa al SMILES: SDF/MOL con enlaces
    plip_xml: Optional[str] = None     # report.xml ya calculado (evita ejecutar PLIP)
    chain: Optional[str] = None        # cadena del ligando si hay varias copias
    resnum: Optional[str] = None       # número de residuo del ligando si hay varias copias
    charge: int = 0                    # solo se usa sin SMILES/SDF
    title: str = ""
    name: str = ""
    workdir: Optional[str] = None      # carpeta de trabajo de PLIP (se reutiliza si el informe es más nuevo que el PDB)
    plip_exe: Optional[str] = None
    compare_xml: Optional[str] = None      # report.xml de la estructura de referencia (p. ej. la pose de docking)
    compare_label: str = ""                # nombre de la referencia en la leyenda (p. ej. "docking")
    compare_resname: Optional[str] = None  # código del ligando en el informe de referencia (por defecto el mismo)
    as_salt: Optional[List[str]] = None    # residuos (p. ej. "ASP12") cuyos enlaces de H se dibujan como puente salino (criterio del usuario)


def _pdb_ligand_coords(pdb: str, resname: str):
    """{serie PDB: (x, y, z)} de los átomos HETATM del ligando."""
    out = {}
    for ln in open(pdb, encoding="utf-8", errors="replace"):
        if ln.startswith(("HETATM", "ATOM")) and ln[17:20].strip() == resname:
            try:
                out[int(ln[6:11])] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                pass
    return out


def fix_serial_offset(interactions, pdb: str, resname: str, tol: float = 0.2):
    """Si los números de serie del XML no corresponden a las coordenadas de ``<ligcoo>`` en el PDB, estima el desplazamiento
    (p. ej. +1 cuando PLIP renumera) por coincidencia de coordenadas y lo aplica. Devuelve el desplazamiento aplicado (0 si no)."""
    coords = _pdb_ligand_coords(pdb, resname)
    if not coords:
        return 0

    def near(c, ser):
        q = coords.get(ser)
        return q is not None and sum((a - b) ** 2 for a, b in zip(q, c)) ** 0.5 <= tol
    probes = [(i.serials[0], i.coo) for i in interactions if i.coo and len(i.serials) == 1]
    if not probes or all(near(c, s) for s, c in probes):
        return 0
    votes = {}
    for s, c in probes:
        for ser, q in coords.items():
            if sum((a - b) ** 2 for a, b in zip(q, c)) ** 0.5 <= tol:
                votes[ser - s] = votes.get(ser - s, 0) + 1
    if not votes:
        return 0
    d = max(votes, key=votes.get)
    for i in interactions:
        i.serials = [x + d for x in i.serials]
    return d


def spec_from_dict(d: Dict[str, Any], base_dir: str = ".", default_workdir: Optional[str] = None,
                   plip_exe: Optional[str] = None) -> SystemSpec:
    """Crea un SystemSpec desde un diccionario JSON; las rutas relativas se resuelven respecto a ``base_dir``."""
    valid = {f.name for f in fields(SystemSpec)}
    alias = {"ligand": "ligand_resname", "resname": "ligand_resname", "report": "plip_xml", "xml": "plip_xml"}
    dd = {alias.get(k, k): v for k, v in d.items()}
    bad = sorted(set(dd) - valid)
    if bad:
        raise ValueError(f"Claves desconocidas en el sistema: {', '.join(bad)}. Válidas: {', '.join(sorted(valid))}")
    if "ligand_resname" not in dd:
        raise ValueError("Cada sistema necesita 'ligand_resname' (código de residuo del ligando en el PDB).")
    for k in ("pdb", "sdf", "plip_xml", "workdir", "compare_xml"):
        if dd.get(k) and not os.path.isabs(dd[k]):
            dd[k] = os.path.normpath(os.path.join(base_dir, dd[k]))
    if not dd.get("workdir") and default_workdir:
        dd["workdir"] = os.path.join(default_workdir, dd.get("name") or dd.get("title") or "sistema")
    if plip_exe and not dd.get("plip_exe"):
        dd["plip_exe"] = plip_exe
    return SystemSpec(**dd)


def scene_from_spec(spec: SystemSpec, style: Style) -> Scene:
    """Ejecuta o lee PLIP, construye el ligando y diseña el esquema."""
    notes: List[str] = []
    if spec.plip_xml:
        report = spec.plip_xml
        fixed = find_fixed_pdb(os.path.dirname(os.path.abspath(report)))
        pdb = fixed or spec.pdb
        if fixed is None and spec.pdb:
            notes.append("No hay 'plipfixed*.pdb' junto al report.xml: se usa el PDB indicado; sus números de serie "
                         "deben coincidir con los del informe.")
    else:
        if not spec.pdb:
            raise ValueError("Falta el PDB del complejo (--pdb) o un report.xml (--plip-xml).")
        workdir = spec.workdir or os.path.join(os.path.dirname(os.path.abspath(spec.pdb)),
                                               "deku2d_plip_" + os.path.splitext(os.path.basename(spec.pdb))[0])
        res = run_plip(spec.pdb, workdir, spec.plip_exe)
        report, pdb = res["report"], res["fixed"]
    if not pdb:
        raise ValueError("No se encontró un PDB con el ligando: indique --pdb.")
    site = parse_report(report, spec.ligand_resname, spec.chain, spec.resnum)
    if spec.as_salt:
        want = {(x[:3].upper(), int(x[3:].rstrip("ABCDEFGH") or 0)) for x in spec.as_salt}
        n_re = 0
        for i in site.interactions:
            if i.kind == "hbond" and (i.res, i.num) in want:
                i.kind, i.prot_donor = "salt", None
                n_re += 1
        if n_re:
            notes.append(f"{n_re} enlace(s) de H de {', '.join(sorted(spec.as_salt))} se dibujan como puente salino "
                         "(reclasificación manual; PLIP los informa como enlaces de H).")
    off = fix_serial_offset(site.interactions, pdb, site.hetid or spec.ligand_resname)
    if off:
        notes.append(f"Los números de serie del informe no coinciden con el PDB: se aplicó un desplazamiento de {off:+d} "
                     "deducido de las coordenadas (<ligcoo>).")
    lig = build_ligand(pdb, spec.ligand_resname, spec.smiles, spec.sdf, spec.chain, spec.resnum, spec.charge)
    sc = build_scene(lig, site.interactions, style, title=spec.title, name=spec.name or spec.ligand_resname)
    sc.warnings = notes + sc.warnings
    if spec.compare_xml:
        ref = parse_report(spec.compare_xml, spec.compare_resname or spec.ligand_resname)
        ref_res = {i.key for i in ref.interactions}
        cur = {k[:3] for k in sc.keys}
        sc.compare = {"label": spec.compare_label or "reference", "new": sorted(cur - ref_res),
                      "lost": sorted(ref_res - cur)}
        if sc.compare["lost"] and style.mark_diff:
            x0, y0, x1, y1 = sc.bounds
            sc.bounds = (x0, y0 - 0.9, x1, y1)              # espacio para la nota de residuos solo en la referencia
    if site.ignored:
        sc.warnings.append("PLIP informó interacciones de tipos que este estilo no dibuja: "
                           + ", ".join(f"{k}×{v}" for k, v in site.ignored.items()))
    return sc
