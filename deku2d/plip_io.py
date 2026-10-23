"""Entrada/salida de PLIP: ejecutar PLIP (subproceso) y leer ``report.xml``.

PLIP 2.x no se usa como módulo: se llama al ejecutable (por defecto el del entorno ``plip_env``).
Los índices de átomo del XML son números de serie del PDB **corregido** que escribe PLIP junto al informe
(``plipfixed*.pdb``), no necesariamente los del PDB original; por eso ``run_plip`` devuelve también esa ruta.
"""
from __future__ import annotations

import glob
import os
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .config import DEFAULT_PLIP_EXE

# etiquetas del XML que el estilo no dibuja (se cuentan y se avisa)
IGNORED_TAGS = ("water_bridge", "halogen_bond", "metal_complex")


@dataclass
class Interaction:
    kind: str                       # hydrophobic | hbond | salt | pipi | picat
    res: str                        # tipo de residuo en mayúsculas (ASP, TYR…)
    num: int                        # número de residuo
    chain: str                      # cadena ("" si no consta)
    serials: List[int]              # números de serie PDB (PLIP) de los átomos del ligando implicados
    dist: float                     # Å (D···A en enlaces de H; centroide en π)
    prot_donor: Optional[bool] = None   # solo enlaces de H: True si el dador es la proteína
    coo: Optional[tuple] = None         # coordenadas (x, y, z) del átomo del ligando según el XML (<ligcoo>), si es un solo átomo

    @property
    def key(self):
        return (self.res, self.num, self.chain)


@dataclass
class BindingSite:
    hetid: str
    chain: str
    position: str
    interactions: List[Interaction] = field(default_factory=list)
    ignored: Dict[str, int] = field(default_factory=dict)


def find_fixed_pdb(outdir: str) -> Optional[str]:
    """PDB corregido por PLIP en ``outdir`` (el que usan los índices del XML)."""
    c = sorted(glob.glob(os.path.join(outdir, "plipfixed*.pdb")))
    return c[0] if c else None


def run_plip(pdb: str, outdir: str, exe: Optional[str] = None, reuse: bool = True) -> Dict[str, str]:
    """Ejecuta ``plip -f PDB -x -t -o OUTDIR`` y devuelve ``{"report": ..., "fixed": ...}``.

    Con ``reuse=True`` y un ``report.xml`` más nuevo que el PDB no se vuelve a ejecutar PLIP.
    """
    exe = exe or DEFAULT_PLIP_EXE
    rep = os.path.join(outdir, "report.xml")
    if reuse and os.path.exists(rep) and os.path.getmtime(rep) >= os.path.getmtime(pdb) and find_fixed_pdb(outdir):
        return {"report": rep, "fixed": find_fixed_pdb(outdir)}
    if not os.path.exists(exe):
        raise FileNotFoundError(
            f"No se encuentra el ejecutable de PLIP: {exe}. Indique --plip-exe o la variable DEKU2D_PLIP, "
            "o aporte un report.xml ya calculado con --plip-xml.")
    os.makedirs(outdir, exist_ok=True)
    r = subprocess.run([exe, "-f", os.path.abspath(pdb), "-x", "-t", "-o", outdir], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(rep):
        raise RuntimeError(f"PLIP falló (código {r.returncode}):\n{(r.stderr or r.stdout)[-800:]}")
    return {"report": rep, "fixed": find_fixed_pdb(outdir) or pdb}


def list_sites(xml_path: str) -> List[BindingSite]:
    """Sitios de unión del informe (hetid, cadena, posición) sin leer las interacciones."""
    root = ET.parse(xml_path).getroot()
    return [BindingSite(bs.findtext("identifiers/hetid") or "", bs.findtext("identifiers/chain") or "",
                        bs.findtext("identifiers/position") or "") for bs in root.iter("bindingsite")]


def _ints(it: ET.Element, tag: str) -> List[int]:
    return [int(x.text) for x in it.iter(tag) if x.text and x.text.strip().lstrip("-").isdigit()]


def parse_report(xml_path: str, ligand_resname: Optional[str] = None, chain: Optional[str] = None,
                 resnum: Optional[str] = None) -> BindingSite:
    """Lee el sitio de unión del ligando y devuelve sus interacciones en el formato interno.

    Si el informe tiene varios sitios con el mismo ``ligand_resname`` se toma el primero que cumpla ``chain``/``resnum``.
    Sin ``ligand_resname`` solo se admite un único sitio.
    """
    root = ET.parse(xml_path).getroot()
    sites = list(root.iter("bindingsite"))
    sel = []
    for bs in sites:
        hid = bs.findtext("identifiers/hetid") or ""
        ch = bs.findtext("identifiers/chain") or ""
        pos = bs.findtext("identifiers/position") or ""
        if ligand_resname and hid != ligand_resname:
            continue
        if chain and ch != chain:
            continue
        if resnum and str(pos) != str(resnum):
            continue
        sel.append(bs)
    if not sel:
        avail = ", ".join(f"{(b.findtext('identifiers/hetid') or '?')}:{b.findtext('identifiers/chain') or '?'}"
                          f"{b.findtext('identifiers/position') or '?'}" for b in sites) or "ninguno"
        raise ValueError(f"No hay un sitio de unión que coincida con el ligando '{ligand_resname}' en {xml_path}. "
                         f"Disponibles: {avail}")
    if not ligand_resname and len(sel) > 1:
        names = ", ".join(sorted({b.findtext("identifiers/hetid") or "?" for b in sel}))
        raise ValueError(f"El informe tiene {len(sel)} sitios ({names}); indique --ligand-resname.")
    bs = sel[0]
    site = BindingSite(bs.findtext("identifiers/hetid") or "", bs.findtext("identifiers/chain") or "",
                       bs.findtext("identifiers/position") or "")
    ia = bs.find("interactions")
    if ia is None:
        return site

    def _coo(it):
        try:
            return tuple(float(it.findtext(f"ligcoo/{c}")) for c in "xyz")
        except (TypeError, ValueError):
            return None

    def add(kind, it, serials, dist, prot_donor=None, coo=None):
        site.interactions.append(Interaction(kind, (it.findtext("restype") or "").upper(), int(it.findtext("resnr")),
                                             it.findtext("reschain") or "", serials, dist, prot_donor, coo))
    for it in ia.iter("hydrophobic_interaction"):
        add("hydrophobic", it, [int(it.findtext("ligcarbonidx"))], float(it.findtext("dist")), coo=_coo(it))
    for it in ia.iter("hydrogen_bond"):
        pd = (it.findtext("protisdon") == "True")
        add("hbond", it, [int(it.findtext("acceptoridx" if pd else "donoridx"))], float(it.findtext("dist_d-a")), pd, coo=_coo(it))
    for tag, kind in (("salt_bridge", "salt"), ("pi_stack", "pipi"), ("pi_cation_interaction", "picat")):
        for it in ia.iter(tag):
            # solo los átomos del ligando (<lig_idx_list>); si el XML no la trae, todos los <idx> y se filtran luego con el mapa
            lig_idx = [int(x.text) for x in it.findall("lig_idx_list/idx")] or _ints(it, "idx")
            add(kind, it, lig_idx, float(it.findtext("dist") or it.findtext("centdist") or 0))
    for tag in IGNORED_TAGS:
        n = len(list(ia.iter(tag)))
        if n:
            site.ignored[tag] = n
    return site
