"""Diseño del esquema: geometría 2D del ligando, zonas hidrofóbicas, colocación de las cajas de residuos y de las
etiquetas de distancia. No dibuja nada; produce un ``Scene`` que ``render`` pinta.

Colocación de cajas (igual que el script original):
  * ranuras sobre una envolvente que sigue el contorno de la molécula (dos capas);
  * búsqueda local (reinicios aleatorios + movimientos de una caja) que minimiza una suma de penalizaciones:
    capa de la ranura, longitud de las líneas, cruces con enlaces del ligando (×12), distancia a los arcos hidrofóbicos,
    cruces entre líneas (×5), líneas que atraviesan cajas ajenas (×15) y cajas demasiado próximas;
  * restricción dura: las cajas no se solapan.
El coste se evalúa con numpy (tablas unarias precalculadas + términos de pares vectorizados), con el mismo orden de
sorteos aleatorios que el original: en modo ``legacy`` (30×2500, sin mejoras) el resultado es idéntico.

Mejoras (``Style.improved``): las líneas grises caja→arcos que cruzan el ligando se penalizan (×12), y las etiquetas de
distancia se colocan por búsqueda de la posición con menos solapes (otras etiquetas, cajas, líneas, enlaces).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem

from .config import Style
from .ligand import Ligand
from .plip_io import Interaction


# ---------------------------------------------------------------- geometría básica
def _ccw(a, b, c):
    return (c[..., 1] - a[..., 1]) * (b[..., 0] - a[..., 0]) > (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])


def seg_cross(p1, p2, q1, q2):
    """¿Se cruzan los segmentos p1p2 y q1q2? (admite difusión de numpy en los últimos ejes)."""
    return (_ccw(p1, q1, q2) != _ccw(p2, q1, q2)) & (_ccw(p1, p2, q1) != _ccw(p1, p2, q2))


def clip_box(c, hw, hh, toward):
    """Punto del borde de la caja (centro c, semiejes hw, hh) en la dirección de ``toward``."""
    v = np.asarray(toward, float) - c
    if abs(v[0]) < 1e-9 and abs(v[1]) < 1e-9:
        return np.asarray(c, float)
    t = min(hw / abs(v[0]) if abs(v[0]) > 1e-9 else 1e9, hh / abs(v[1]) if abs(v[1]) > 1e-9 else 1e9)
    return c + v * t


def seg_rect_hit(p, q, centre, half) -> bool:
    """¿El segmento pq toca el rectángulo de centro ``centre`` y semiejes ``half``?"""
    cx, cy = centre; w, h = half
    for a in (p, q):
        if abs(a[0] - cx) < w and abs(a[1] - cy) < h:
            return True
    cs = [np.array(v) for v in ((cx - w, cy - h), (cx + w, cy - h), (cx + w, cy + h), (cx - w, cy + h))]
    return any(bool(seg_cross(np.asarray(p), np.asarray(q), cs[i], cs[(i + 1) % 4])) for i in range(4))


def _seg_rect_arr(p, q, centres, half):
    """Segmentos p[s]→q[s] (S,2) frente a rectángulos de centros ``centres`` (L,2) y semiejes ``half``: matriz (S,L)."""
    P = p[:, None, :]
    Q = q[:, None, :]
    cx = centres[None, :, 0]; cy = centres[None, :, 1]
    w, h = half

    def inside(A):
        return (np.abs(A[..., 0] - cx) < w) & (np.abs(A[..., 1] - cy) < h)
    corners = np.stack([np.stack([cx - w, cy - h], -1), np.stack([cx + w, cy - h], -1),
                        np.stack([cx + w, cy + h], -1), np.stack([cx - w, cy + h], -1)], 2)     # (1,L,4,2)
    e2 = np.roll(corners, -1, axis=2)
    cr = seg_cross(P[:, :, None, :], Q[:, :, None, :], corners, e2).any(axis=2)
    return inside(P) | inside(Q) | cr


def rect_overlap(c1, h1, c2, h2) -> bool:
    return abs(c1[0] - c2[0]) < h1[0] + h2[0] and abs(c1[1] - c2[1]) < h1[1] + h2[1]


# ---------------------------------------------------------------- modelos
@dataclass
class Contact:
    kind: str
    res: str
    num: int
    chain: str
    idx: List[int]                  # átomos RDKit del ligando
    dist: float
    prot_donor: Optional[bool] = None

    @property
    def key(self):
        return (self.res, self.num, self.chain)


@dataclass
class Cluster:
    center: np.ndarray
    r0: float
    atoms: List[int]
    is_ring: bool
    ident: tuple


@dataclass
class PolarLine:
    key: tuple
    kind: str
    prot_donor: Optional[bool]
    dist: float
    j: int
    atom: int
    start: np.ndarray               # extremo en la caja
    end: np.ndarray                 # extremo en el ligando
    color_key: str                  # red|green|purple|cyan|orange
    label_pos: Optional[np.ndarray] = None


@dataclass
class ArcSet:
    key: tuple
    center: np.ndarray
    r0: float
    ang: float                      # grados
    span: float
    is_ring: bool
    tip: np.ndarray                 # extremo de la línea gris sobre los arcos
    start: np.ndarray               # extremo de la línea gris en la caja
    dist: str = ""                  # texto con las distancias de los contactos hidrofóbicos de este grupo


@dataclass
class Scene:
    name: str
    title: str
    mol: Chem.Mol
    xy: np.ndarray
    rings: List[List[int]]
    ring_c: List[np.ndarray]
    halo_atoms: set
    labelled: set
    keys: List[tuple]
    pos: Dict[tuple, np.ndarray]
    hw: Dict[tuple, float]
    hh: float
    polar: List[PolarLine]
    arcs: List[ArcSet]
    contacts: List[Contact]
    bounds: Tuple[float, float, float, float]       # x0, y0, x1, y1 (Å) sin leyenda
    used_kinds: set
    seen_dir: set
    cost: float
    stats: Dict[str, float] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    style: Optional[Style] = None
    compare: Optional[Dict[str, object]] = None     # {'label', 'new': [claves], 'lost': [residuos]} frente a la referencia


# ---------------------------------------------------------------- construcción
def box_label(k: tuple, keys: Sequence[tuple]) -> str:
    """Texto de una caja; si el residuo tiene varias cajas se numeran (THR58A·1, THR58A·2)."""
    base = f"{k[0]}{k[1]}{k[2]}"
    n = sum(1 for q in keys if q[:3] == k[:3])
    return f"{base}·{k[3] + 1}" if n > 1 else base


def contacts_from_interactions(inter: Sequence[Interaction], lig: Ligand) -> Tuple[List[Contact], List[str]]:
    """Traduce números de serie PLIP a átomos RDKit y descarta lo que no toca al ligando."""
    s2i = lig.serial_to_idx
    out, dropped = [], 0
    for it in inter:
        idx = [s2i[s] for s in it.serials if s in s2i]
        if not idx:
            dropped += 1
            continue
        out.append(Contact(it.kind, it.res, it.num, it.chain, idx, it.dist, it.prot_donor))
    warns = []
    if dropped:
        warns.append(f"{dropped} interacciones sin átomos del ligando identificables (serie PDB fuera del mapa); "
                     "compruebe que el PDB de entrada sea el 'plipfixed' del mismo informe.")
    return out, warns


class _Search:
    """Evaluador vectorizado del coste y búsqueda local con el mismo orden de sorteos que el script original."""

    def __init__(self, slots, layer, keys, hw, hh, targets, unary):
        self.slots, self.n = slots, len(keys)
        self.S = len(slots)
        self.hw = np.array([hw[k] for k in keys])
        self.hh = hh
        self.U = unary
        self.T = np.array([t for k in keys for t in targets[k]], float).reshape(-1, 2)
        self.Lo = np.array([i for i, k in enumerate(keys) for _ in targets[k]], int)
        self.m = len(self.Lo)
        self.ii, self.jj = np.triu_indices(self.m, 1)
        self.pi, self.pj = np.triu_indices(self.n, 1)
        self.evals = 0

    def cost(self, assign) -> float:
        self.evals += 1
        a = np.asarray(assign, int)
        n = self.n
        C = self.slots[a]
        c = float(self.U[np.arange(n), a].sum())
        if self.m:
            P = C[self.Lo]
            ii, jj = self.ii, self.jj
            if len(ii):
                far = np.linalg.norm(P[ii] - P[jj], axis=1) > 0.05
                c += 5.0 * np.count_nonzero(far & seg_cross(P[ii], self.T[ii], P[jj], self.T[jj]))
            w = (self.hw + 0.15)[None, :]; h = self.hh + 0.15
            ins_p = (np.abs(P[:, None, 0] - C[None, :, 0]) < w) & (np.abs(P[:, None, 1] - C[None, :, 1]) < h)
            ins_t = (np.abs(self.T[:, None, 0] - C[None, :, 0]) < w) & (np.abs(self.T[:, None, 1] - C[None, :, 1]) < h)
            cx, cy = C[:, 0], C[:, 1]; ww, hh_ = self.hw + 0.15, h
            corners = np.stack([np.stack([cx - ww, cy - hh_], -1), np.stack([cx + ww, cy - hh_], -1),
                                np.stack([cx + ww, cy + hh_], -1), np.stack([cx - ww, cy + hh_], -1)], 1)   # (n,4,2)
            e2 = np.roll(corners, -1, axis=1)
            edge = seg_cross(P[:, None, None, :], self.T[:, None, None, :], corners[None], e2[None]).any(axis=2)
            hit = ins_p | ins_t | edge
            hit &= self.Lo[:, None] != np.arange(n)[None, :]
            c += 15.0 * np.count_nonzero(hit)
        if n > 1:
            d = np.linalg.norm(C[self.pi] - C[self.pj], axis=1)
            c += 1.2 * float(np.sum((3.6 - d)[d < 3.6]))
        return c

    def feasible(self, assign) -> bool:
        C = self.slots[np.asarray(assign, int)]
        if self.n < 2:
            return True
        dx = np.abs(C[self.pi, 0] - C[self.pj, 0]); dy = np.abs(C[self.pi, 1] - C[self.pj, 1])
        return not bool(np.any((dx < self.hw[self.pi] + self.hw[self.pj] + 0.35) & (dy < 2 * self.hh + 0.3)))

    def run(self, seed, restarts, iters):
        rng = np.random.default_rng(seed)
        n, S = self.n, self.S
        best, bc = None, 1e18
        for _ in range(restarts):
            cur = list(rng.choice(S, n, replace=False))
            tries = 0
            while not self.feasible(cur) and tries < 4000:
                cur = list(rng.choice(S, n, replace=False)); tries += 1
            if not self.feasible(cur):
                continue
            cc = self.cost(cur)
            for _ in range(iters):
                i = int(rng.integers(n)); new = list(cur); new[i] = int(rng.integers(S))
                if len(set(new)) < n or not self.feasible(new):
                    continue
                nc = self.cost(new)
                if nc < cc:
                    cur, cc = new, nc
            if cc < bc:
                best, bc = cur, cc
        return best, bc


def build_scene(lig: Ligand, interactions: Sequence[Interaction], style: Style, title: str = "", name: str = "") -> Scene:
    t0 = time.time()
    mol = lig.mol
    warns = list(lig.warnings)
    contacts, w2 = contacts_from_interactions(interactions, lig)
    warns += w2
    FS, PPU = style.fs, style.ppu
    r_at = style.atom_pt / style.fs                          # las etiquetas de átomo ocupan más si la letra es mayor
    ATOM_HALF = (0.38 * r_at, 0.32 * r_at)

    conf = mol.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])
    centre = xy.mean(0)
    rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
    ring_c = [xy[r].mean(0) for r in rings]
    nb = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()}
    bonds_a = np.array([xy[b.GetBeginAtomIdx()] for b in mol.GetBonds()], float).reshape(-1, 2)
    bonds_b = np.array([xy[b.GetEndAtomIdx()] for b in mol.GetBonds()], float).reshape(-1, 2)
    labelled = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "C"}

    def outward(i):
        """Bisectriz del mayor hueco angular entre enlaces, favoreciendo apuntar lejos del centro y no hacia anillos."""
        away = (xy[i] - centre) / (np.linalg.norm(xy[i] - centre) + 1e-9)
        if not nb[i]:
            return away
        angs = sorted(math.atan2(*(xy[j] - xy[i])[::-1]) for j in nb[i])
        if len(angs) == 1:
            return -np.array([math.cos(angs[0]), math.sin(angs[0])])
        best, bsc = None, -1e9
        for k in range(len(angs)):
            g0 = angs[k]; g1 = angs[(k + 1) % len(angs)] + (2 * math.pi if k == len(angs) - 1 else 0)
            m = (g0 + g1) / 2; bis = np.array([math.cos(m), math.sin(m)])
            sc = (g1 - g0) + 1.2 * float(np.dot(bis, away))
            for r in rings:
                if i in r:
                    rc = xy[r].mean(0) - xy[i]
                    sc -= 6 * max(0.0, float(np.dot(bis, rc / (np.linalg.norm(rc) + 1e-9))))
            if sc > bsc:
                best, bsc = bis, sc
        return best

    # ---- residuos, zonas hidrofóbicas (centro: anillo o átomo) y halo ----
    resid: Dict[tuple, List[Contact]] = {}
    for c in contacts:
        resid.setdefault(c.key + (0,), []).append(c)
    keys = sorted(resid)
    n = len(keys)
    hclusters: Dict[tuple, List[Cluster]] = {}
    halo_atoms: set = set()
    for k in keys:
        hat = sorted({i for c in resid[k] if c.kind == "hydrophobic" for i in c.idx})
        used, cl = set(), []
        for a in hat:
            if a in used:
                continue
            cand = [r for r in range(len(rings)) if a in rings[r]]
            if cand:
                r = max(cand, key=lambda q: len(set(rings[q]) & set(hat)))
                at = list(rings[r])
                used |= (set(rings[r]) & set(hat)) | {a}
                rad = float(np.mean([np.linalg.norm(xy[x] - ring_c[r]) for x in at]))
                cl.append(Cluster(ring_c[r], rad + 0.5, at, True, ("r", r)))
                halo_atoms |= set(at)
            else:
                used.add(a); halo_atoms.add(a)
                cl.append(Cluster(xy[a], 0.6, [a], False, ("a", a)))
        hclusters[k] = cl

    # ---- ranuras sobre la envolvente del ligando (dos capas) ----
    M = 1.75
    NS = 80 if style.legacy else 120                        # ranuras angulares por capa (más finas que el original)
    ang_s = 2 * np.pi * np.arange(NS) / NS
    U_dir = np.stack([np.cos(ang_s), np.sin(ang_s)], 1)                       # (NS,2)
    d = xy - centre
    proj = d @ U_dir.T
    perp = np.abs(d[:, [0]] * U_dir[None, :, 1] - d[:, [1]] * U_dir[None, :, 0])
    val = np.where(perp < M, proj + np.sqrt(np.maximum(M * M - perp ** 2, 0)), 0.0)
    hull = np.maximum(val.max(axis=0), 0.0)
    slots, layer = [], []
    for lay in (0, 1):
        for s_ in range(NS):
            u = U_dir[s_]
            slots.append(centre + u * (hull[s_] + 0.55 + 0.8 * abs(u[0]) + 0.35 * abs(u[1]) + 1.6 * lay))
            layer.append(lay)
    slots = np.array(slots); layer = np.array(layer)
    S = len(slots)
    def target_of(c: Contact):
        return xy[c.idx].mean(0) if c.kind in ("pipi", "picat") and len(c.idx) > 1 else xy[c.idx[0]]
    bond_ends = np.array([[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()], int).reshape(-1, 2)

    def skip_of(c: Contact):
        """Enlaces del propio anillo de un contacto π: la línea llega al borde del anillo, no a su centro."""
        if c.kind in ("pipi", "picat") and len(c.idx) > 1 and len(bond_ends):
            s_ = set(c.idx)
            return np.array([a in s_ and b in s_ for a, b in bond_ends])
        return None

    def count_cross(P, T_, excl=True, skip=None):
        """Cruces de la línea P(S,2)→T(2,) con cada enlace del ligando, por ranura."""
        if not len(bonds_a):
            return np.zeros(len(P))
        p = P[:, None, :]; t = np.broadcast_to(T_, (len(P), 1, 2))
        cr = seg_cross(p, t, bonds_a[None], bonds_b[None])
        if excl:
            cr &= (np.linalg.norm(bonds_a - T_, axis=1) > 0.05)[None] & (np.linalg.norm(bonds_b - T_, axis=1) > 0.05)[None]
        if skip is not None:
            cr &= ~skip[None]
        return cr.sum(axis=1).astype(float)

    if not style.legacy:
        # Un residuo con contactos polares en lados opuestos del ligando no puede llegar a todos sin cruzarlo:
        # se dibuja una caja por contacto (convención habitual en los diagramas 2D).
        for k in list(resid):
            pol = [c for c in resid[k] if c.kind != "hydrophobic"]
            if len(pol) < 2:
                continue
            tot = sum(count_cross(slots, target_of(c), skip=skip_of(c)) for c in pol)
            if tot.min() == 0:
                continue
            alone = [min(count_cross(slots, target_of(c), skip=skip_of(c))) == 0 for c in pol]
            if not all(alone):
                continue
            rest = [c for c in resid[k] if c.kind == "hydrophobic"] + [pol[0]]
            resid[k] = rest
            for sub, c in enumerate(pol[1:], 1):
                resid[k[:3] + (sub,)] = [c]
                hclusters[k[:3] + (sub,)] = []
            keys = sorted(resid)
            n = len(keys)

    hw = {k: (0.58 * FS * len(box_label(k, keys)) + 12) / PPU / 2 for k in keys}
    hh = (FS * 1.45 + 6) / PPU / 2

    targets = {k: [target_of(c) for c in resid[k] if c.kind != "hydrophobic"] for k in keys}
    tskips = {k: [None if style.legacy else skip_of(c) for c in resid[k] if c.kind != "hydrophobic"] for k in keys}

    anchors = {}
    for k in keys:
        a_ = []
        for cl in hclusters[k]:
            dirv = (cl.center - centre) / (np.linalg.norm(cl.center - centre) + 1e-9) if cl.is_ring else outward(cl.atoms[0])
            a_.append(cl.center + dirv * (cl.r0 + 1.1))
        anchors[k] = a_

    # tablas unarias U[k, s]
    XW = 12.0 if style.legacy else 40.0                     # penalización de cada cruce con un enlace del ligando
    lab_xy = np.array([xy[i] for i in sorted(labelled)], float).reshape(-1, 2)
    U = np.zeros((max(n, 1), S))
    for ki, k in enumerate(keys):
        u = 0.6 * layer.astype(float)
        for t, sk in zip(targets[k], tskips[k]):
            u += 0.15 * np.linalg.norm(slots - t, axis=1) + XW * count_cross(slots, t, skip=sk)
        for an in anchors[k]:
            u += 0.9 * np.linalg.norm(slots - an, axis=1)
        if style.improved:                                # las líneas grises caja→arcos no deben cruzar el ligando
            for cl in hclusters[k]:
                v = slots - cl.center
                nv = np.linalg.norm(v, axis=1, keepdims=True) + 1e-9
                tips = cl.center + (cl.r0 + 0.32) * v / nv
                if len(bonds_a):
                    cr = seg_cross(tips[:, None, :], slots[:, None, :], bonds_a[None], bonds_b[None])
                    u += 2.5 * XW * cr.sum(axis=1)
                if len(lab_xy):
                    u += 10.0 * _seg_rect_arr(tips, slots, lab_xy, ATOM_HALF).sum(axis=1)
            for t in targets[k]:                          # ni las líneas polares las etiquetas de otros átomos (OH, O⁻…)
                if len(lab_xy):
                    h_ = _seg_rect_arr(slots, np.broadcast_to(t, slots.shape), lab_xy, ATOM_HALF)
                    h_ &= (np.linalg.norm(lab_xy - t, axis=1) > 0.05)[None, :]
                    u += 10.0 * h_.sum(axis=1)
        U[ki] = u

    cost_final, evals, pos = 0.0, 0, {}
    if n:
        srch = _Search(slots, layer, keys, hw, hh, targets, U)
        best, bc = srch.run(style.seed, style.restarts, style.iterations)
        evals = srch.evals
        if best is None:                                  # ninguna combinación factible: reparto uniforme en la capa interior
            warns.append("No se encontró una colocación sin solapes de cajas; se usa un reparto uniforme (puede haber solapes).")
            best = [int(round(i * NS / n)) % NS for i in range(n)]
            bc = srch.cost(best)
        pos = {k: slots[s] for k, s in zip(keys, best)}
        cost_final = bc

    # ---- límites ----
    pts = [xy] + [np.array([pos[k] - [hw[k], hh], pos[k] + [hw[k], hh]]) for k in keys]
    allp = np.vstack(pts)
    x0, y0 = allp.min(0) - 1.0
    x1, y1 = allp.max(0) + 1.0

    # ---- líneas polares ----
    polar: List[PolarLine] = []
    for k in keys:
        pol = [c for c in resid[k] if c.kind != "hydrophobic"]
        for j, c in enumerate(pol):
            tgt = target_of(c)
            off = (j - (len(pol) - 1) / 2) * (0.24 if style.legacy else 0.40)
            v = tgt - pos[k]; nrm = np.array([-v[1], v[0]]) / (np.linalg.norm(v) + 1e-9)
            start = clip_box(pos[k], hw[k], hh, tgt) + nrm * off
            gap = 0.30 if c.idx[0] in labelled and c.kind in ("hbond", "salt") else 0.05
            end = tgt - (tgt - start) / (np.linalg.norm(tgt - start) + 1e-9) * gap + nrm * off * 0.7
            ck = {"salt": "purple", "pipi": "cyan", "picat": "orange"}.get(c.kind)
            if c.kind == "hbond":
                ck = "red" if c.prot_donor else "green"
            polar.append(PolarLine(k, c.kind, c.prot_donor, c.dist, j, c.idx[0], start, end, ck))

    # ---- arcos hidrofóbicos (separación angular entre juegos que comparten centro) y línea gris ----
    by_center: Dict[tuple, List[Tuple[tuple, Cluster]]] = {}
    for k in keys:
        for cl in hclusters[k]:
            by_center.setdefault(cl.ident, []).append((k, cl))
    arcs: List[ArcSet] = []

    def _vdw_text(k, cl):
        ds = sorted({c.dist for c in resid[k] if c.kind == "hydrophobic" and set(c.idx) & set(cl.atoms)})
        return ", ".join(f"{d:.1f}" for d in ds)
    for ck, lst in by_center.items():
        m_ = len(lst)
        base = 78 if lst[0][1].is_ring else 100
        span = min(base, 340.0 / m_ - 10)
        th = [math.degrees(math.atan2(*(pos[k] - cl.center)[::-1])) for k, cl in lst]
        order = sorted(range(m_), key=lambda i: th[i])
        ang = [th[i] for i in order]
        for _ in range(400):
            moved = False
            for i in range(m_):
                j = (i + 1) % m_
                gap = (ang[j] - ang[i]) % 360 if m_ > 1 else 360
                if m_ > 1 and gap < span + 8:
                    push = (span + 8 - gap) / 2 + 0.05
                    ang[i] -= push; ang[j] += push; moved = True
            if not moved:
                break
        for idx, i in enumerate(order):
            k, cl = lst[i]
            a_ = math.radians(ang[idx])
            tip = cl.center + (cl.r0 + 0.2 + 0.12) * np.array([math.cos(a_), math.sin(a_)])
            arcs.append(ArcSet(k, cl.center, cl.r0, ang[idx], span, cl.is_ring, tip, clip_box(pos[k], hw[k], hh, tip), _vdw_text(k, cl)))
            if not style.legacy and len(bonds_a):
                # el reparto angular puede desviar el arco hacia el ligando: se busca el ángulo libre de cruces más cercano
                def crosses(a_deg, k=k, cl=cl):
                    r_ = math.radians(a_deg)
                    t_ = cl.center + (cl.r0 + 0.32) * np.array([math.cos(r_), math.sin(r_)])
                    return int(seg_cross(t_[None], clip_box(pos[k], hw[k], hh, t_)[None], bonds_a, bonds_b).sum())
                if crosses(ang[idx]):
                    others = [ang[j] for j in range(m_) if j != idx]
                    for dlt in sorted(range(-130, 131, 5), key=abs):
                        a2 = ang[idx] + dlt
                        if crosses(a2) == 0 and all(abs((a2 - o + 180) % 360 - 180) >= 0.62 * span for o in others):
                            ang[idx] = a2
                            r_ = math.radians(a2)
                            tip = cl.center + (cl.r0 + 0.2 + 0.12) * np.array([math.cos(r_), math.sin(r_)])
                            arcs[-1] = ArcSet(k, cl.center, cl.r0, a2, span, cl.is_ring, tip, clip_box(pos[k], hw[k], hh, tip), _vdw_text(k, cl))
                            break

    sc = Scene(name or lig.name, title, mol, xy, rings, ring_c, halo_atoms, labelled, keys, pos, hw, hh, polar, arcs,
               contacts, (x0, y0, x1, y1), {c.kind for c in contacts},
               {p.prot_donor for p in polar if p.kind == "hbond"}, cost_final, warnings=warns, style=style)
    place_labels(sc, bonds_a, bonds_b)
    sc.stats = {"search_s": time.time() - t0, "evals": float(evals), "n_boxes": float(n)}
    return sc


# ---------------------------------------------------------------- etiquetas de distancia
def place_labels(sc: Scene, bonds_a: np.ndarray, bonds_b: np.ndarray) -> None:
    """Fija ``label_pos`` de cada línea con distancia (enlaces de H y puentes salinos)."""
    st = sc.style
    lab = [p for p in sc.polar if p.kind in ("hbond", "salt")]
    if not st.dist:
        return
    if st.legacy:                                         # comportamiento del script original
        for p in lab:
            f = 0.40 + 0.2 * (p.j % 3)
            m = p.start + (p.end - p.start) * f
            ln = (p.end - p.start) / (np.linalg.norm(p.end - p.start) + 1e-9)
            p.label_pos = m + np.array([-ln[1], ln[0]]) * 0.30
        return
    fs = st.fs - 5
    boxes = [(sc.pos[k], (sc.hw[k] + 0.1, sc.hh + 0.1)) for k in sc.keys]
    r_at = st.atom_pt / st.fs
    atom_rects = [(sc.xy[i], (0.38 * r_at, 0.32 * r_at)) for i in sc.labelled]
    segs = [(p.start, p.end) for p in sc.polar]
    placed: List[Tuple[np.ndarray, Tuple[float, float]]] = []
    npol: Dict[tuple, int] = {}
    for q in sc.polar:
        npol[q.key] = npol.get(q.key, 0) + 1
    for p in lab:
        txt = f"{p.dist:.1f}"
        half = ((0.6 * fs * len(txt) + 5) / st.ppu / 2, (fs * 1.3 + 3) / st.ppu / 2)
        ln = (p.end - p.start) / (np.linalg.norm(p.end - p.start) + 1e-9)
        nv = np.array([-ln[1], ln[0]])
        best, bs = None, 1e9
        for fi, f in enumerate((0.5, 0.4, 0.6, 0.32, 0.68, 0.25, 0.75, 0.18, 0.82)):
            for side in (1, -1):
                pos = p.start + (p.end - p.start) * f + nv * side * 0.30
                sc_ = 0.4 * abs(f - 0.5) + (0.05 if side < 0 else 0.0)
                if npol.get(p.key, 1) > 1:                  # líneas paralelas del mismo residuo: cada distancia por su lado exterior
                    sc_ += 2.0 * (side != (-1 if p.j < npol[p.key] / 2 else 1))
                sc_ += 10 * sum(rect_overlap(pos, half, c2, h2) for c2, h2 in placed)
                sc_ += 10 * sum(rect_overlap(pos, half, c2, h2) for c2, h2 in boxes)
                sc_ += 6 * sum(rect_overlap(pos, half, c2, h2) for c2, h2 in atom_rects)
                sc_ += 3 * sum(seg_rect_hit(a, b, pos, half) for a, b in segs if a is not p.start)
                if len(bonds_a):
                    sc_ += 2 * sum(seg_rect_hit(a, b, pos, half) for a, b in zip(bonds_a, bonds_b))
                if sc_ < bs:
                    best, bs = pos, sc_
        p.label_pos = best
        placed.append((best, half))
