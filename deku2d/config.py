"""Parámetros por defecto del estilo de figura y textos por idioma.

El estilo reproduce el de ``estilo.png`` (referencia del usuario) tal como quedó en ``plip_2d_estilo.py``:
ligando negro con halo amarillo y arcos naranjas en los contactos hidrofóbicos, línea punteada gris clara caja→arcos,
cajas grises con nombre en mayúsculas y cadena, enlaces de H con línea punteada y flecha triangular abierta
(rojo: dador en la proteína; verde: dador en el ligando), puentes salinos morados, π–π cian y π–catión naranja.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, fields, replace
from typing import Any, Dict, Optional

DEFAULT_PLIP_EXE = os.environ.get("DEKU2D_PLIP", "/home/alecastro/anaconda3/envs/plip_env/bin/plip")

# Textos de la leyenda por idioma. Las claves son estables; los valores se pueden sobrescribir desde JSON.
LANG: Dict[str, Dict[str, str]] = {
    "en": {
        "hb_prot": "H-bond, protein donor",
        "hb_lig": "H-bond, ligand donor",
        "salt": "Salt bridge",
        "pipi": "π–π stacking",
        "picat": "π–cation",
        "vdw": "Hydrophobic / vdW contact",
        "region": "Hydrophobic region",
        "acid": "Acidic residue",
        "basic": "Basic residue",
        "polar": "Polar residue",
        "apolar": "Non-polar residue",
        "new": "Not in {ref}",
        "lost": "In {ref} only: {list}",
    },
    "es": {
        "hb_prot": "Enlace de H, dador en la proteína",
        "hb_lig": "Enlace de H, dador en el ligando",
        "salt": "Puente salino",
        "pipi": "Apilamiento π–π",
        "picat": "π–catión",
        "vdw": "Contacto hidrofóbico / vdW",
        "region": "Región hidrofóbica",
        "acid": "Residuo ácido",
        "basic": "Residuo básico",
        "polar": "Residuo polar",
        "apolar": "Residuo apolar",
        "new": "No está en {ref}",
        "lost": "Solo en {ref}: {list}",
    },
}


@dataclass
class Style:
    """Todos los parámetros que cambian el aspecto o el comportamiento del diseño."""

    # tamaño
    fs: float = 17.0            # letra de residuos (la leyenda usa fs-1; las distancias fs-5)
    atom_fs: Optional[float] = None   # letra de los átomos del ligando; por defecto fs+2 (fs en modo --legacy)
    ppu: float = 34.0           # puntos por Å: fija la escala figura/letra (menor = letra relativamente mayor)
    min_width: float = 21.0     # ancho mínimo de la figura (Å) para que quepa la leyenda
    dpi: int = 300
    # colores
    halo: str = "#FFE29A"
    arc: str = "#E3A008"
    grey: str = "#B5B5B5"
    red: str = "#FF4040"        # enlace de H con dador en la proteína
    green: str = "#2CA02C"      # enlace de H con dador en el ligando
    purple: str = "#8E44AD"     # puente salino
    cyan: str = "#1BA3B8"       # apilamiento π–π
    orange: str = "#E67E22"     # π–catión
    box_edge: str = "#8C8C8C"
    tint_acid: str = "#F8D9D4"      # relleno de las cajas según el tipo de residuo
    tint_basic: str = "#D6E6F8"
    tint_polar: str = "#DDF0E3"
    tint_apolar: str = "#F2F2F2"
    new_edge: str = "#1A1A1A"       # borde de las cajas de residuos que no están en la estructura de referencia
    # trazos (puntos tipográficos)
    bond_lw: float = 2.1
    link_lw: float = 3.0
    arc_lw: float = 2.6
    grey_lw: float = 1.9
    box_lw: float = 1.4
    new_lw: float = 3.4
    arrow_size: float = 0.15    # Å
    # contenido
    dist: bool = True           # distancia (Å) en enlaces de H y puentes salinos
    legend: bool = True
    tint: bool = True           # relleno de las cajas por tipo de residuo (ácido, básico, polar, apolar)
    vdw_dist: bool = False      # distancia (Å) junto a las líneas grises hidrofóbicas (útil en material suplementario)
    mark_diff: bool = True      # marca los residuos que no están en la estructura de referencia (ver compare_xml)
    uniform: bool = True        # en paneles, celdas del mismo tamaño y ligandos centrados en ellas
    lang: str = "en"
    # diseño
    seed: int = 3
    restarts: int = 30          # reinicios de la búsqueda local (el script original usaba 30)
    iterations: int = 2500      # movimientos por reinicio (el script original usaba 2500)
    improved: bool = True       # etiquetas de distancia sin solapes + líneas grises que esquivan el ligando
    texts: Optional[Dict[str, str]] = None   # sobrescribe textos de la leyenda

    def text(self, key: str) -> str:
        base = LANG.get(self.lang, LANG["en"])
        if self.texts and key in self.texts:
            return self.texts[key]
        return base[key]

    @property
    def atom_pt(self) -> float:
        """Tamaño de letra de los átomos del ligando."""
        if self.atom_fs is not None:
            return self.atom_fs
        return self.fs if self.legacy else self.fs + 2

    def tint_of(self, resname: str) -> str:
        """Color de relleno de una caja según el tipo de residuo."""
        r = resname.upper()
        if r in ("ASP", "GLU"):
            return self.tint_acid
        if r in ("LYS", "ARG", "HIS"):
            return self.tint_basic
        if r in ("SER", "THR", "ASN", "GLN", "TYR", "CYS"):
            return self.tint_polar
        return self.tint_apolar

    @property
    def legacy(self) -> bool:
        return not self.improved


def legacy_style(**kw: Any) -> Style:
    """Estilo que reproduce exactamente la figura del script original (sin las mejoras de etiquetas y líneas)."""
    return Style(improved=False, **kw)


def style_from_dict(d: Optional[Dict[str, Any]], base: Optional[Style] = None) -> Style:
    """Aplica un diccionario (p. ej. leído de JSON) sobre un Style; claves desconocidas dan un error claro."""
    st = base or Style()
    if not d:
        return st
    valid = {f.name for f in fields(Style)}
    bad = sorted(set(d) - valid)
    if bad:
        raise ValueError(f"Parámetros de estilo desconocidos: {', '.join(bad)}. Válidos: {', '.join(sorted(valid))}")
    return replace(st, **d)


def load_json(path: str) -> Dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)
