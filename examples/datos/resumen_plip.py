"""Resume los report.xml de PLIP 2.4.0 de las poses de docking y de los fotogramas representativos de la MD.
Salida: resumen_plip.json (todas las interacciones) y manuscrito/notas/plip_2d_datos.md (tablas descriptivas).
Uso: python3 resumen_plip.py"""
import os, json
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
MAN = '/home/alecastro/calc/dsuazo/manuscrito/'
SYSTEMS = ['mol1', 'mol3', 'mol8', 'mol31', 'MRTX1133']
RESN = {s: 'UNK' for s in SYSTEMS}; RESN['MRTX1133'] = '6IC'
TYPES = [('hydrophobic', 'Hidrofóbicas'), ('hbond', 'Enlaces de H'), ('salt', 'Puentes salinos'), ('pipi', 'π–π'),
         ('picat', 'π–catión'), ('water', 'Puentes de agua')]
TAGS = {'hydrophobic_interaction': 'hydrophobic', 'hydrogen_bond': 'hbond', 'salt_bridge': 'salt', 'pi_stack': 'pipi',
        'pi_cation_interaction': 'picat', 'water_bridge': 'water'}


def parse(path, het):
    root = ET.parse(path).getroot()
    for bs in root.iter('bindingsite'):
        if bs.findtext('identifiers/hetid') != het:
            continue
        out = []
        for tag, key in TAGS.items():
            for it in bs.iter(tag):
                d = it.findtext('dist') or it.findtext('centdist') or it.findtext('dist_a-w') or it.findtext('dist_d-a')
                if key == 'hbond':
                    d = it.findtext('dist_d-a')
                out.append(dict(tipo=key, residuo=f"{it.findtext('restype')}{it.findtext('resnr')}",
                                dist=float(d) if d else None,
                                dador=('proteína' if it.findtext('protisdon') == 'True' else 'ligando') if key == 'hbond' else None))
        return out
    raise RuntimeError(path)


def fmt(x):
    return f'{x:.1f}'.replace('.', ',')


def res_dist(items):
    by = {}
    for i in items:
        by.setdefault(i['residuo'], []).append(i)
    parts = []
    for r in sorted(by, key=lambda s: (int(''.join(c for c in s if c.isdigit())), s)):
        ds = ', '.join(fmt(i['dist']) + (f" {i['dador'][0]}→{'l' if i['dador'][0] == 'p' else 'p'}" if i['dador'] else '') for i in by[r])
        parts.append(f'{r} ({ds})')
    return '; '.join(parts) if parts else '—'


if __name__ == '__main__':
    rep = json.load(open(MAN + 'pymol/md/representativos.json'))
    data = {}
    for s in SYSTEMS:
        for st, f in (('docking', 'docking/%s/report.xml'), ('md_rep', 'md_rep/%s/report.xml')):
            p = os.path.join(HERE, f % s)
            if os.path.exists(p):
                data.setdefault(s, {})[st] = parse(p, RESN[s])
    json.dump(data, open(os.path.join(HERE, 'resumen_plip.json'), 'w'), indent=1, ensure_ascii=False)

    L = ['# Interacciones PLIP en las poses de docking y en los fotogramas representativos de la dinámica',
         '',
         '**Alcance.** PLIP 2.4.0 es una detección geométrica sobre **un solo fotograma** (o una sola pose): no mide ocupación ni '
         'persistencia, y no indica afinidad. Las diferencias entre docking y fotograma representativo describen dos estructuras '
         'concretas. Los datos y los `report.xml` están en `/home/alecastro/Documentos/proyecto_2Deku_scheme/ejemplos/datos/` '
         '(`docking/<sistema>/`, `md_rep/<sistema>/`); `resumen_plip.py` regenera este documento.',
         '',
         '## Estructuras analizadas',
         '',
         '- **Docking:** mol1 = `best_BREED/md_7rpz_mol1_g2-1_200ns/xundrug/bestpose1_g2-1.pdb`; mol3, mol8 y mol31 = '
         '`manuscrito/pymol/docking/inputs/complex_<sistema>.pdb`; MRTX1133 = pose de re-docking de FRED '
         '(`KRAS/g12d/mrtx1133/redock/fred_docked.sdf`) sobre el receptor `7rpz_clean.pdb`, que no contiene GDP ni Mg²⁺ y cuya pose '
         'de FRED no tiene hidrógenos (PLIP/Open Babel los añadió); el complejo se armó aquí (`preparar_datos.py`).',
         '- **Fotograma representativo de la MD:** medoide del grupo más poblado (agrupamiento tipo GROMOS, umbral de 1,0 Å, RMSD de '
         'los átomos pesados del ligando tras alinear los Cα 1–169) de los últimos 100 ns de la única trayectoria de 200 ns de cada sistema '
         '(`pymol/scripts/cluster_representative.py`). En MRTX1133 el sitio de unión incluye GDP y Mg²⁺.',
         '',
         '| Sistema | Fotograma | Tiempo (ns) | Grupo más poblado | Ventana (fotogramas) |', '|---|---|---|---|---|']
    for s in SYSTEMS:
        r = rep[s]
        L.append(f"| {s} | {r['frame']} | {fmt(r['t_ns'])} | {r['tam_cluster1']}/{r['n_ventana']} ({fmt(100 * r['frac_cluster1'])} %) | {r['n_ventana']} |")
    L += ['', '## Número de interacciones por tipo', '',
          '| Sistema | Estructura | ' + ' | '.join(n for _, n in TYPES) + ' |', '|---|---|' + '---|' * len(TYPES)]
    for s in SYSTEMS:
        for st, nm in (('docking', 'Docking'), ('md_rep', 'MD rep.')):
            if st in data.get(s, {}):
                it = data[s][st]
                cells = []
                for k, _ in TYPES:
                    sel = [i for i in it if i['tipo'] == k]
                    cells.append(f"{len(sel)} ({len({i['residuo'] for i in sel})})" if sel else '0')
                L.append(f'| {s} | {nm} | ' + ' | '.join(cells) + ' |')
    L += ['', 'Entre paréntesis: número de residuos distintos. Distancias en Å: dador–aceptor en los enlaces de H, centro–centro en π–π, '
          'carga–centro en π–catión y puentes salinos, átomo de C a átomo de C en los contactos hidrofóbicos. '
          'En los enlaces de H, «p→l» indica dador en la proteína y «l→p» dador en el ligando.', '']
    for s in SYSTEMS:
        L += [f'## {s}', '']
        for st, nm in (('docking', 'Docking'), ('md_rep', 'Fotograma representativo')):
            if st not in data.get(s, {}):
                continue
            L += [f'**{nm}**', '']
            for k, name in TYPES:
                sel = [i for i in data[s][st] if i['tipo'] == k]
                if sel:
                    L.append(f'- {name}: {res_dist(sel)}')
            L.append('')
        d, m = data[s].get('docking', []), data[s].get('md_rep', [])
        L += ['**Docking frente a fotograma representativo (por residuo y tipo)**', '']
        for k, name in TYPES:
            a = {i['residuo'] for i in d if i['tipo'] == k}; b = {i['residuo'] for i in m if i['tipo'] == k}
            if not (a or b):
                continue
            key = lambda r: int(''.join(c for c in r if c.isdigit()))
            L.append(f"- {name}: en ambos: {', '.join(sorted(a & b, key=key)) or '—'}; solo en docking: {', '.join(sorted(a - b, key=key)) or '—'}; "
                     f"solo en el fotograma: {', '.join(sorted(b - a, key=key)) or '—'}")
        L.append('')
    L += ['## Limitaciones', '',
          '- Un fotograma por sistema y una trayectoria por sistema: no hay réplicas ni estadística sobre la variación entre fotogramas.',
          '- Las poses de docking y los fotogramas de la MD parten de entornos distintos (en la MD, solvente, iones y GDP–Mg²⁺ en la simulación); '
          'la comparación es entre dos estructuras, no entre dos poblaciones.',
          '- En el docking de MRTX1133 el receptor no tiene GDP ni Mg²⁺ y el ligando no traía hidrógenos; las diferencias con el fotograma '
          'de la MD no deben atribuirse solo a la dinámica.',
          '- PLIP aplica criterios geométricos propios (umbrales de distancia y ángulo), distintos de los usados en la tabla de ocupación '
          'de contactos de la dinámica (C–C ≤ 4,0 Å; N/O···N/O ≤ 3,5 Å) y de los diagramas de Schrödinger.', '']
    open(MAN + 'notas/plip_2d_datos.md', 'w', encoding='utf-8').write('\n'.join(L))
    print('escrito', MAN + 'notas/plip_2d_datos.md', len(L), 'líneas')
