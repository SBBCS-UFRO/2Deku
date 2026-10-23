"""Línea de órdenes: ``python3 -m deku2d {draw,panel,batch}`` (o el lanzador ``deku2d``)."""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import replace
from typing import Any, Dict, List, Optional

from . import __version__
from .config import Style, legacy_style, load_json, style_from_dict
from .layout import Scene
from .panel import render_panel
from .pipeline import SystemSpec, scene_from_spec, spec_from_dict
from .render import render_single, save_figure


# ---------------------------------------------------------------- opciones comunes
def _common(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("estilo y salida")
    g.add_argument("--lang", choices=["en", "es"], help="idioma de la leyenda (por defecto en)")
    g.add_argument("--dist", action=argparse.BooleanOptionalAction, default=None,
                   help="distancias (Å) en enlaces de H y puentes salinos (por defecto sí)")
    g.add_argument("--fontsize", type=float, help="tamaño de letra de residuos y átomos (por defecto 17)")
    g.add_argument("--atom-fontsize", type=float, help="tamaño de letra de los átomos del ligando (por defecto fontsize + 2)")
    g.add_argument("--ppu", type=float, help="puntos por Å: menor = letra relativamente mayor (por defecto 34)")
    g.add_argument("--dpi", type=int, help="resolución de PNG (por defecto 300)")
    g.add_argument("--format", dest="formats", help="formatos de salida separados por coma: png,svg,pdf")
    g.add_argument("--min-width", type=float, help="ancho mínimo de la figura en Å (por defecto 21)")
    g.add_argument("--seed", type=int, help="semilla de la búsqueda de colocación (por defecto 3)")
    g.add_argument("--restarts", type=int, help="reinicios de la búsqueda local (por defecto 30)")
    g.add_argument("--iterations", type=int, help="movimientos por reinicio (por defecto 2500)")
    g.add_argument("--legacy", action="store_true",
                   help="reproduce el script original (sin las mejoras de etiquetas de distancia y líneas grises)")
    g.add_argument("--tint", action=argparse.BooleanOptionalAction, default=None,
                   help="relleno de las cajas por tipo de residuo: ácido, básico, polar, apolar (por defecto sí)")
    g.add_argument("--vdw-dist", action=argparse.BooleanOptionalAction, default=None,
                   help="distancias (Å) en las líneas grises hidrofóbicas (por defecto no)")
    g.add_argument("--mark-diff", action=argparse.BooleanOptionalAction, default=None,
                   help="marcar con borde grueso los residuos que no están en la referencia (por defecto sí)")
    g.add_argument("--uniform", action=argparse.BooleanOptionalAction, default=None,
                   help="paneles con celdas del mismo tamaño y ligandos centrados (por defecto sí)")
    g.add_argument("--no-legend", action="store_true", help="sin leyenda")
    g.add_argument("--style-json", help="JSON con parámetros de estilo (claves de deku2d.config.Style)")
    g.add_argument("--plip-exe", help="ruta del ejecutable de PLIP (o variable DEKU2D_PLIP)")
    g.add_argument("-q", "--quiet", action="store_true", help="no imprimir el resumen")


def _style(args: argparse.Namespace, cfg_style: Optional[Dict[str, Any]] = None) -> Style:
    """Prioridad: valores por defecto < 'style' del JSON de configuración < --style-json < opciones de la línea de órdenes."""
    st = legacy_style() if getattr(args, "legacy", False) else Style()
    st = style_from_dict(cfg_style, st)
    if getattr(args, "style_json", None):
        st = style_from_dict(load_json(args.style_json), st)
    ov: Dict[str, Any] = {}
    for opt, key in (("lang", "lang"), ("dist", "dist"), ("tint", "tint"), ("vdw_dist", "vdw_dist"), ("mark_diff", "mark_diff"), ("uniform", "uniform"), ("fontsize", "fs"), ("atom_fontsize", "atom_fs"), ("ppu", "ppu"), ("dpi", "dpi"),
                     ("min_width", "min_width"), ("seed", "seed"), ("restarts", "restarts"), ("iterations", "iterations")):
        v = getattr(args, opt, None)
        if v is not None:
            ov[key] = v
    if getattr(args, "no_legend", False):
        ov["legend"] = False
    return replace(st, **ov)


def _formats(args: argparse.Namespace, cfg: Optional[Dict[str, Any]] = None) -> Optional[List[str]]:
    f = getattr(args, "formats", None)
    if f:
        return [x.strip() for x in f.split(",") if x.strip()]
    if cfg and cfg.get("formats"):
        return list(cfg["formats"])
    return None


def _report(sc: Scene, paths: List[str], quiet: bool) -> None:
    if quiet:
        return
    print(f"{sc.name or 'sistema'}: {int(sc.stats.get('n_boxes', 0))} residuos, {len(sc.contacts)} interacciones, "
          f"coste {sc.cost:.1f}, {sc.stats.get('search_s', 0):.1f} s -> " + ", ".join(paths))
    for w in sc.warnings:
        print("  aviso:", w, file=sys.stderr)


# ---------------------------------------------------------------- órdenes
def cmd_draw(args: argparse.Namespace) -> int:
    st = _style(args)
    out = args.output
    base = os.path.splitext(out)[0]
    spec = SystemSpec(ligand_resname=args.ligand_resname, pdb=args.pdb, smiles=args.smiles, sdf=args.sdf,
                      plip_xml=args.plip_xml, chain=args.chain, resnum=args.resnum, charge=args.charge,
                      title=args.title or "", name=os.path.basename(base),
                      workdir=args.workdir or (base + "_plip"), plip_exe=args.plip_exe,
                      compare_xml=args.compare_xml, compare_label=args.compare_label, as_salt=args.as_salt)
    sc = scene_from_spec(spec, st)
    paths = save_figure(render_single(sc, st), out, _formats(args), st.dpi)
    _report(sc, paths, args.quiet)
    return 0


def cmd_panel(args: argparse.Namespace) -> int:
    cfg = load_json(args.config)
    st = _style(args, cfg.get("style"))
    base_dir = os.path.dirname(os.path.abspath(args.config))
    plip_exe = args.plip_exe or cfg.get("plip_exe")
    wd = cfg.get("workdir")
    if wd and not os.path.isabs(wd):
        wd = os.path.normpath(os.path.join(base_dir, wd))
    default_wd = wd or os.path.join(os.path.dirname(os.path.abspath(args.output)), os.path.splitext(os.path.basename(args.output))[0] + "_plip")
    scenes = []
    for d in cfg.get("systems", []):
        spec = spec_from_dict(d, base_dir, default_wd, plip_exe)
        sc = scene_from_spec(spec, st)
        scenes.append(sc)
        if not args.quiet:
            print(f"  {sc.name}: {int(sc.stats.get('n_boxes', 0))} residuos, {sc.stats.get('search_s', 0):.1f} s")
            for w in sc.warnings:
                print("    aviso:", w, file=sys.stderr)
    ncols = args.ncols or cfg.get("ncols", 2)
    fig = render_panel(scenes, st, ncols=ncols, labels=not args.no_labels and cfg.get("labels", True),
                       label_style=args.label_style or cfg.get("label_style", "a)"),
                       show_titles=not args.no_titles and cfg.get("titles", True))
    paths = save_figure(fig, args.output, _formats(args, cfg), st.dpi)
    if not args.quiet:
        print("panel ->", ", ".join(paths))
    return 0


def cmd_batch(args: argparse.Namespace) -> int:
    cfg = load_json(args.config)
    st = _style(args, cfg.get("style"))
    base_dir = os.path.dirname(os.path.abspath(args.config))
    outdir = args.outdir or cfg.get("outdir", "deku2d_salida")
    if not os.path.isabs(outdir):
        outdir = os.path.normpath(os.path.join(base_dir, outdir))
    plip_exe = args.plip_exe or cfg.get("plip_exe")
    wd = cfg.get("workdir") or os.path.join(outdir, "plip")
    if not os.path.isabs(wd):
        wd = os.path.normpath(os.path.join(base_dir, wd))
    fmts = _formats(args, cfg)
    scenes = []
    for d in cfg.get("systems", []):
        spec = spec_from_dict(d, base_dir, wd, plip_exe)
        sc = scene_from_spec(spec, st)
        scenes.append(sc)
        name = spec.name or spec.title or spec.ligand_resname
        paths = save_figure(render_single(sc, st), os.path.join(outdir, f"{name}.png"), fmts, st.dpi)
        _report(sc, paths, args.quiet)
    pc = cfg.get("panel")
    if pc and scenes:
        fig = render_panel(scenes, st, ncols=pc.get("ncols", 2), labels=pc.get("labels", True),
                           label_style=pc.get("label_style", "a)"), show_titles=pc.get("titles", True))
        paths = save_figure(fig, os.path.join(outdir, pc.get("output", "panel") + ".png"), fmts, st.dpi)
        if not args.quiet:
            print("panel ->", ", ".join(paths))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="deku2d", description="Esquemas 2D de interacciones ligando–receptor para artículos "
                                 "científicos (RDKit + matplotlib; interacciones de PLIP).")
    ap.add_argument("--version", action="version", version=f"deku2d {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("draw", help="un esquema desde un complejo PDB (o un report.xml de PLIP)")
    d.add_argument("--pdb", help="PDB del complejo (si hay --plip-xml se usa el 'plipfixed' contiguo al XML)")
    d.add_argument("--ligand-resname", "--resname", dest="ligand_resname", required=True,
                   help="código de residuo del ligando en el PDB (p. ej. UNK, 6IC)")
    d.add_argument("--smiles", help="SMILES plantilla con los mismos átomos pesados (define enlaces y cargas)")
    d.add_argument("--sdf", help="SDF/MOL del ligando con enlaces (alternativa al SMILES; se asocia por coordenadas)")
    d.add_argument("--plip-xml", help="report.xml de PLIP ya calculado")
    d.add_argument("--chain", help="cadena del ligando si hay varias copias")
    d.add_argument("--resnum", help="número de residuo del ligando si hay varias copias")
    d.add_argument("--charge", type=int, default=0, help="carga neta (solo sin SMILES/SDF)")
    d.add_argument("--title", help="título opcional sobre la figura")
    d.add_argument("--as-salt", action="append", metavar="RESNUM",
                   help="residuo (p. ej. ASP12) cuyo enlace de H se dibuja como puente salino; se puede repetir")
    d.add_argument("--compare-xml", help="report.xml de la estructura de referencia: marca los residuos que no están en ella")
    d.add_argument("--compare-label", default="", help="nombre de la referencia en la leyenda (p. ej. docking)")
    d.add_argument("--workdir", help="carpeta de trabajo de PLIP (por defecto <salida>_plip)")
    d.add_argument("-o", "--output", required=True, help="archivo de salida (.png, .svg o .pdf)")
    _common(d)
    d.set_defaults(func=cmd_draw)

    p = sub.add_parser("panel", help="panel a), b), c)… desde un JSON con la lista de sistemas")
    p.add_argument("--config", required=True, help="JSON: {systems:[{title,pdb,ligand_resname,smiles|sdf,plip_xml}], ncols, style}")
    p.add_argument("-o", "--output", required=True, help="archivo de salida (.png, .svg o .pdf)")
    p.add_argument("--ncols", type=int, help="paneles por fila (por defecto 2)")
    p.add_argument("--no-labels", action="store_true", help="sin rótulos a), b), …")
    p.add_argument("--label-style", help="a) | A) | (a) | a.")
    p.add_argument("--no-titles", action="store_true", help="no escribir el título de cada panel")
    _common(p)
    p.set_defaults(func=cmd_panel)

    b = sub.add_parser("batch", help="una figura por sistema (y opcionalmente un panel) desde un JSON")
    b.add_argument("--config", required=True, help="JSON: {outdir, formats, style, systems:[…], panel:{output,ncols}}")
    b.add_argument("--outdir", help="carpeta de salida (sobrescribe 'outdir' del JSON)")
    _common(b)
    b.set_defaults(func=cmd_batch)
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except (ValueError, FileNotFoundError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
