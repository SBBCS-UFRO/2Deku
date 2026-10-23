"""Figura de varios paneles (a), b), c)…) con una sola leyenda compartida y escala común.

Todos los paneles usan los mismos ``ppu`` y tamaños de letra, de modo que la letra mide lo mismo en cada uno;
cada panel ocupa exactamente la extensión de su esquema. Los paneles se reparten en filas de ``ncols`` y cada fila
se centra horizontalmente.
"""
from __future__ import annotations

import math
import string
from typing import List, Optional, Sequence

import matplotlib.pyplot as plt

from .config import Style
from .layout import Scene
from .render import draw_scene, legend_entries, legend_width_pt, rc_context, scene_classes, scene_ref


def panel_labels(n: int, style: str = "a)") -> List[str]:
    """Rótulos a), b), … (admite 'a)', 'A)', '(a)' o 'a.')."""
    letters = string.ascii_lowercase
    out = []
    for i in range(n):
        ch = letters[i % 26] * (1 + i // 26)
        if style.startswith("A"):
            ch = ch.upper()
        out.append(f"({ch})" if style.startswith("(") else f"{ch}." if style.endswith(".") else f"{ch})")
    return out


def render_panel(scenes: Sequence[Scene], st: Style, ncols: int = 2, labels: bool = True, label_style: str = "a)",
                 show_titles: bool = True, gap_in: float = 0.3, margin_in: float = 0.2) -> plt.Figure:
    if not scenes:
        raise ValueError("No hay esquemas para el panel.")
    with rc_context():
        kinds, seen = set(), set()
        for s in scenes:
            kinds |= s.used_kinds
            seen |= s.seen_dir
        items = legend_entries(kinds, seen, st, scene_classes(scenes, st), scene_ref(scenes, st)) if st.legend else []
        lims = []
        for s in scenes:
            lims.append(tuple(s.bounds))
        if st.uniform and len(scenes) > 1:          # misma celda para todos y ligando centrado: misma escala y altura visual
            cen = [s.xy.mean(0) for s in scenes]
            hx = max(max(c[0] - b[0], b[2] - c[0]) for c, b in zip(cen, lims))
            hy = max(max(c[1] - b[1], b[3] - c[1]) for c, b in zip(cen, lims))
            lims = [(c[0] - hx, c[1] - hy, c[0] + hx, c[1] + hy) for c in cen]
        sizes = [((b[2] - b[0]) * st.ppu / 72, (b[3] - b[1]) * st.ppu / 72) for b in lims]
        rows = [list(range(i, min(i + ncols, len(scenes)))) for i in range(0, len(scenes), ncols)]
        row_w = [sum(sizes[i][0] for i in r) + gap_in * (len(r) - 1) for r in rows]
        row_h = [max(sizes[i][1] for i in r) for r in rows]
        ncol_leg = min(3, len(items)) or 1
        leg_w = legend_width_pt(items, ncol_leg, st.fs) / 72
        while ncol_leg > 1 and leg_w + 2 * margin_in > max(row_w):     # reducir columnas hasta que quepa
            ncol_leg -= 1
            leg_w = legend_width_pt(items, ncol_leg, st.fs) / 72
        leg_h = (math.ceil(len(items) / ncol_leg) * st.fs * 1.7 + 12) / 72 if items else 0.0
        W = max(max(row_w), leg_w) + 2 * margin_in
        H = sum(row_h) + gap_in * (len(rows) - 1) + leg_h + 2 * margin_in
        fig = plt.figure(figsize=(W, H))
        lbls = panel_labels(len(scenes), label_style) if labels else [""] * len(scenes)
        y_top = H - margin_in
        for r, rh in zip(rows, row_h):
            x = (W - sum(sizes[i][0] for i in r) - gap_in * (len(r) - 1)) / 2
            for i in r:
                w, h = sizes[i]
                ax = fig.add_axes([x / W, (y_top - (rh + h) / 2) / H, w / W, h / H])
                x0, y0, x1, y1 = lims[i]
                ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.axis("off")
                draw_scene(ax, scenes[i], st)
                if lbls[i]:
                    ax.text(x0 + 0.1, y1 - 0.1, lbls[i], ha="left", va="top", fontsize=st.fs + 3, fontweight="bold", zorder=9)
                    if show_titles and scenes[i].title:
                        off = (len(lbls[i]) * 0.65 * (st.fs + 3) + 8) / st.ppu
                        ax.text(x0 + 0.1 + off, y1 - 0.18, scenes[i].title, ha="left", va="top", fontsize=st.fs, zorder=9)
                elif show_titles and scenes[i].title:
                    ax.text(x0 + 0.1, y1 - 0.18, scenes[i].title, ha="left", va="top", fontsize=st.fs, zorder=9)
                x += w + gap_in
            y_top -= rh + gap_in
        if items:
            fig.legend(*zip(*items), loc="lower center", bbox_to_anchor=(0.5, margin_in / H), ncol=ncol_leg,
                       fontsize=st.fs - 1, frameon=False, handlelength=2.6, columnspacing=1.8)
        return fig
