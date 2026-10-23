"""Dibujo con matplotlib de un ``Scene`` (estilo de ``estilo.png``) y figura individual con leyenda."""
from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Arc, FancyBboxPatch, Patch, Polygon  # noqa: E402
from rdkit import Chem  # noqa: E402

from .config import Style  # noqa: E402
from .layout import Scene, box_label  # noqa: E402

DOT = (0, (1, 1.8))
RC = {"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "path"}


def rc_context():
    return matplotlib.rc_context(RC)


def _col(st: Style, key: str) -> str:
    return getattr(st, key)


def _tri(ax, tip, frm, col, st: Style):
    v = tip - frm
    v = v / (np.linalg.norm(v) + 1e-9)
    nrm = np.array([-v[1], v[0]])
    size = st.arrow_size
    base = tip - v * size
    ax.add_patch(Polygon([tip, base + nrm * size * 0.5, base - nrm * size * 0.5], closed=True, fc="white", ec=col,
                         lw=2.0, zorder=6, joinstyle="miter"))
    return base


def draw_scene(ax, sc: Scene, st: Optional[Style] = None) -> None:
    """Pinta halo, ligando, arcos, líneas, distancias y cajas en ``ax`` (ya con límites y sin ejes)."""
    st = st or sc.style
    mol, xy = sc.mol, sc.xy
    # halo amarillo en átomos/anillos con contacto hidrofóbico
    hl = st.ppu * 0.55
    for b in mol.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in sc.halo_atoms and c in sc.halo_atoms:
            ax.plot(*zip(xy[a], xy[c]), color=st.halo, lw=hl, solid_capstyle="round", zorder=1)
    for a in sc.halo_atoms:
        ax.plot(*xy[a], "o", color=st.halo, ms=hl, mew=0, zorder=1)

    # ligando (Kekulé; doble enlace interior en anillos)
    kek = Chem.Mol(mol)
    try:
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        kek = Chem.Mol(mol)
    for b in kek.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        p, q = xy[i], xy[j]
        v = q - p
        L = np.linalg.norm(v) + 1e-9
        nrm = np.array([-v[1], v[0]]) / L
        bt = b.GetBondType()
        if bt == Chem.BondType.DOUBLE:
            inring = [r for r in range(len(sc.rings)) if i in sc.rings[r] and j in sc.rings[r]]
            if inring:
                side = 1 if np.dot(sc.ring_c[inring[0]] - (p + q) / 2, nrm) > 0 else -1
                ax.plot(*zip(p, q), color="k", lw=st.bond_lw, solid_capstyle="round", zorder=3)
                ax.plot(*zip(p + v * 0.13 + side * nrm * 0.17, q - v * 0.13 + side * nrm * 0.17), color="k",
                        lw=st.bond_lw, solid_capstyle="round", zorder=3)
            else:
                for s in (+1, -1):
                    ax.plot(*zip(p + s * nrm * 0.07, q + s * nrm * 0.07), color="k", lw=st.bond_lw,
                            solid_capstyle="round", zorder=3)
        elif bt == Chem.BondType.TRIPLE:
            ax.plot(*zip(p, q), color="k", lw=st.bond_lw, solid_capstyle="round", zorder=3)
            for s in (+1, -1):
                ax.plot(*zip(p + s * nrm * 0.1, q + s * nrm * 0.1), color="k", lw=st.bond_lw, solid_capstyle="round", zorder=3)
        else:
            ax.plot(*zip(p, q), color="k", lw=st.bond_lw, solid_capstyle="round", zorder=3)
    for a in mol.GetAtoms():
        e = a.GetSymbol()
        if e == "C":
            continue
        h = a.GetTotalNumHs()
        chg = a.GetFormalCharge()
        txt = e + ("$^{+}$" if chg > 0 else "$^{-}$" if chg < 0 else "") + ("H" if h == 1 else f"H$_{h}$" if h > 1 else "")
        ax.text(*xy[a.GetIdx()], txt, ha="center", va="center", fontsize=st.atom_pt, color="k", zorder=5,
                bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))

    # arcos hidrofóbicos + línea punteada gris caja → arcos
    for a in sc.arcs:
        for m in range(2):
            r = a.r0 + 0.2 * m
            ax.add_patch(Arc(a.center, 2 * r, 2 * r, theta1=a.ang - a.span / 2 + 5 * m, theta2=a.ang + a.span / 2 - 5 * m,
                             color=st.arc, lw=st.arc_lw, capstyle="round", zorder=2))
        ax.plot(*zip(a.start, a.tip), color=st.grey, lw=st.grey_lw, ls=DOT, dash_capstyle="round", zorder=2)
        if st.vdw_dist and a.dist:
            m = a.tip + (a.start - a.tip) * 0.5
            ax.text(*m, a.dist, fontsize=st.fs - 6, color="#7A7A7A", ha="center", va="center", zorder=3,
                    bbox=dict(boxstyle="square,pad=0.03", fc="white", ec="none", alpha=0.85))

    # interacciones polares: línea punteada (con flecha en enlaces de H) y distancia del color de la línea
    for p in sc.polar:
        col = _col(st, p.color_key)
        if p.kind == "hbond":
            if p.prot_donor:                                   # dador en la proteína → flecha en el ligando
                le = _tri(ax, p.end, p.start, col, st)
                ax.plot(*zip(p.start, le), color=col, lw=st.link_lw, ls=DOT, dash_capstyle="round", zorder=4)
            else:                                              # dador en el ligando → flecha en la caja
                tip = p.start - (p.start - p.end) / (np.linalg.norm(p.start - p.end) + 1e-9) * 0.04
                le = _tri(ax, tip, p.end, col, st)
                ax.plot(*zip(p.end, le), color=col, lw=st.link_lw, ls=DOT, dash_capstyle="round", zorder=4)
        else:
            ax.plot(*zip(p.start, p.end), color=col, lw=st.link_lw, ls=DOT, dash_capstyle="round", zorder=4)
            if p.kind in ("pipi", "picat"):
                ax.plot(*p.end, "o", color=col, ms=8, mec="white", mew=1.0, zorder=5)   # punto en el centro del anillo
        if st.dist and p.kind in ("hbond", "salt") and p.label_pos is not None:
            ax.text(*p.label_pos, f"{p.dist:.1f}", fontsize=st.fs - 5, color=col, ha="center", va="center", zorder=6,
                    bbox=dict(boxstyle="square,pad=0.04", fc="white", ec="none", alpha=0.85))

    # cajas de residuos: relleno por tipo de residuo y borde grueso si no están en la estructura de referencia
    cmp_ = sc.compare if (sc.compare and st.mark_diff) else None
    new = set(map(tuple, cmp_["new"])) if cmp_ else set()
    for k in sc.keys:
        is_new = k[:3] in new
        ax.add_patch(FancyBboxPatch((sc.pos[k][0] - sc.hw[k], sc.pos[k][1] - sc.hh), 2 * sc.hw[k], 2 * sc.hh,
                                    boxstyle="round,pad=0,rounding_size=0.08",
                                    fc=st.tint_of(k[0]) if st.tint else "white",
                                    ec=st.new_edge if is_new else st.box_edge, lw=st.new_lw if is_new else st.box_lw, zorder=7))
        ax.text(*sc.pos[k], box_label(k, sc.keys), ha="center", va="center", fontsize=st.fs, color="k", zorder=8)
    if cmp_ and cmp_["lost"]:
        lst = ", ".join(f"{r}{n}{c}" for r, n, c in cmp_["lost"])
        x0, y0 = sc.bounds[0], sc.bounds[1]
        ax.text(x0 + 0.15, y0 + 0.25, st.text("lost").format(ref=cmp_["label"], list=lst), ha="left", va="bottom",
                fontsize=st.fs - 3, color="#444444", zorder=9)


# ---------------------------------------------------------------- leyenda
def legend_entries(kinds: set, seen_dir: set, st: Style, classes: Optional[set] = None,
                   ref: str = "") -> List[Tuple[object, str]]:
    """Entradas de la leyenda para los tipos presentes (en el orden del estilo de referencia)."""
    items: List[Tuple[Line2D, str]] = []
    if True in seen_dir:
        items.append((Line2D([0], [0], color=st.red, lw=st.link_lw, ls=DOT), st.text("hb_prot")))
    if False in seen_dir:
        items.append((Line2D([0], [0], color=st.green, lw=st.link_lw, ls=DOT), st.text("hb_lig")))
    if "salt" in kinds:
        items.append((Line2D([0], [0], color=st.purple, lw=st.link_lw, ls=DOT), st.text("salt")))
    if "pipi" in kinds:
        items.append((Line2D([0], [0], color=st.cyan, lw=st.link_lw, ls=DOT, marker="o", ms=6), st.text("pipi")))
    if "picat" in kinds:
        items.append((Line2D([0], [0], color=st.orange, lw=st.link_lw, ls=DOT, marker="o", ms=6), st.text("picat")))
    if "hydrophobic" in kinds:
        items.append((Line2D([0], [0], color=st.grey, lw=2.2, ls=DOT), st.text("vdw")))
        items.append((Line2D([0], [0], color=st.arc, lw=3.0), st.text("region")))
    if st.tint and classes:
        for key, col in (("acid", st.tint_acid), ("basic", st.tint_basic), ("polar", st.tint_polar), ("apolar", st.tint_apolar)):
            if key in classes:
                items.append((Patch(fc=col, ec=st.box_edge, lw=1.2), st.text(key)))
    if ref:
        items.append((Patch(fc="white", ec=st.new_edge, lw=2.6), st.text("new").format(ref=ref)))
    return items


def scene_classes(scenes: Sequence[Scene], st: Style) -> set:
    """Tipos de residuo presentes (para la leyenda del relleno)."""
    out = set()
    names = {st.tint_acid: "acid", st.tint_basic: "basic", st.tint_polar: "polar", st.tint_apolar: "apolar"}
    for sc in scenes:
        for k in sc.keys:
            out.add(names[st.tint_of(k[0])])
    return out


def scene_ref(scenes: Sequence[Scene], st: Style) -> str:
    """Nombre de la referencia si algún panel marca residuos nuevos."""
    if not st.mark_diff:
        return ""
    for sc in scenes:
        if sc.compare and sc.compare["new"]:
            return str(sc.compare["label"])
    return ""


def legend_width_pt(items, ncol: int, fs: float) -> float:
    """Estimación del ancho de la leyenda en puntos (columnas llenas por columnas, como matplotlib)."""
    if not items:
        return 0.0
    nrow = math.ceil(len(items) / ncol)
    cols = [items[c * nrow:(c + 1) * nrow] for c in range(ncol)]
    f = fs - 1
    w = sum(max((len(t) for _, t in col), default=0) * 0.6 * f + 2.6 * f + 10 for col in cols if col)
    return w + (ncol - 1) * 1.8 * f


# ---------------------------------------------------------------- figura individual
def render_single(sc: Scene, st: Optional[Style] = None) -> plt.Figure:
    st = st or sc.style
    with rc_context():
        x0, y0, x1, y1 = sc.bounds
        items = legend_entries(sc.used_kinds, sc.seen_dir, st, scene_classes([sc], st), scene_ref([sc], st)) if st.legend else []
        ncol = min(2, len(items)) or 1
        need = max(st.min_width, legend_width_pt(items, ncol, st.fs) / st.ppu + 2.0)
        if x1 - x0 < need:
            cx = (x0 + x1) / 2
            x0, x1 = cx - need / 2, cx + need / 2
        y0 -= (0.7 + 1.15 * math.ceil(len(items) / 2)) if items else 0.3
        if sc.title:
            y1 += 1.2
        fig = plt.figure(figsize=((x1 - x0) * st.ppu / 72, (y1 - y0) * st.ppu / 72))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.axis("off")
        draw_scene(ax, sc, st)
        if sc.title:
            ax.text((x0 + x1) / 2, y1 - 0.25, sc.title, ha="center", va="top", fontsize=st.fs + 1, fontweight="bold")
        if items:
            ax.legend(*zip(*items), loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=ncol, fontsize=st.fs - 1,
                      frameon=False, handlelength=2.6, columnspacing=1.8)
        return fig


def save_figure(fig: plt.Figure, out: str, formats: Optional[Iterable[str]] = None, dpi: int = 300) -> List[str]:
    """Guarda ``fig`` en ``out`` (la extensión fija el formato si no se indica ``formats``) y devuelve las rutas."""
    base, ext = os.path.splitext(out)
    fmts = [f.lower().lstrip(".") for f in (formats or [ext.lstrip(".") or "png"])]
    paths = []
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with rc_context():
        for f in fmts:
            if f not in ("png", "svg", "pdf"):
                raise ValueError(f"Formato no soportado: {f} (use png, svg o pdf)")
            p = f"{base}.{f}"
            fig.savefig(p, dpi=dpi, facecolor="white")
            paths.append(p)
    plt.close(fig)
    return paths
