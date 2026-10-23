"""Prepara y verifica los datos de ligandos para los esquemas 2D (ejemplos KRAS G12D / Switch-II).
 - ligandos.json: SMILES, fórmula y carga de cada ligando (mol1, mol3, mol8, mol31 desde manuscrito/scripts/plip_2d.py;
   MRTX1133 desde KRAS/g12d/mrtx1133/redock/fred_docked.sdf).
 - verificacion_ligandos.json: compara fórmula/átomos pesados/H del SMILES con los HETATM del ligando en cada PDB
   (pose de docking y fotograma representativo de la MD) y prueba AssignBondOrdersFromTemplate.
 - docking/MRTX1133/complejo_origen.pdb: receptor de 7rpz_clean.pdb (sin el ligando cristalográfico) + pose de FRED (6IC).
 - md_rep/MRTX1133/ligand.sdf: ligando del fotograma representativo (coordenadas del PDB) con enlaces y cargas de Maestro
   tomados de t155_lig.sdf (mmqmeku) como plantilla.
Uso: python3 preparar_datos.py"""
import os, sys, json, collections
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

HERE = os.path.dirname(os.path.abspath(__file__))
R = '/home/alecastro/calc/dsuazo/'
MAN = R + 'manuscrito/'
sys.path.insert(0, MAN + 'scripts')
import plip_2d as P  # noqa: E402  (solo el dict SMILES)

REDOCK = R + 'KRAS/g12d/mrtx1133/redock/'
MD_T155 = MAN + 'scripts/_out/mmqmeku/frames/MRTX1133/t155_lig.sdf'
SYSTEMS = ['mol1', 'mol3', 'mol8', 'mol31', 'MRTX1133']
RESN = {s: 'UNK' for s in SYSTEMS}; RESN['MRTX1133'] = '6IC'


def het_lines(pdb, resn):
    out = []
    for l in open(pdb):
        if l.startswith('HETATM') and l[17:20].strip() == resn:
            out.append(l)
    return out


def elem(l):
    e = l[76:78].strip()
    return (e or l[12:16].strip()[0]).capitalize()


def smiles_table():
    mr = Chem.RemoveHs(Chem.SDMolSupplier(REDOCK + 'fred_docked.sdf', removeHs=False)[0])
    tab = {s: dict(smiles=P.SMILES[s], fuente='manuscrito/scripts/plip_2d.py (dict SMILES)') for s in P.SMILES}
    tab['mol31']['notas'] = ('Ligando simulado (carga +1, sin S): NO coincide con la molécula n.º 31 de la tabla de las 50 mejores '
                             '(C16H14N6O2S); misma constitución plana que las moléculas 18, 24, 33, 39 y 41.')
    tab['MRTX1133'] = dict(smiles=Chem.MolToSmiles(mr), fuente='KRAS/g12d/mrtx1133/redock/fred_docked.sdf (pose de FRED, mismo SMILES que t155_lig.sdf de la MD)',
                           notas='Protonación de la simulación: +2 (dos aminas protonadas).')
    for s, d in tab.items():
        m = Chem.MolFromSmiles(d['smiles'])
        d.update(formula=CalcMolFormula(m), carga_formal=Chem.GetFormalCharge(m), atomos_pesados=m.GetNumHeavyAtoms())
    return tab


def check(system, struct, pdb, resn, smi):
    lines = het_lines(pdb, resn)
    heavy = [l for l in lines if elem(l) != 'H']
    n_h = len(lines) - len(heavy)
    tpl = Chem.MolFromSmiles(smi)
    tpl_el = collections.Counter(a.GetSymbol() for a in tpl.GetAtoms())
    pdb_el = collections.Counter(elem(l) for l in heavy)
    tpl_h = sum(a.GetTotalNumHs() for a in tpl.GetAtoms())
    row = dict(sistema=system, estructura=struct, pdb=os.path.relpath(pdb, HERE), resname=resn,
               hetatm_total=len(lines), pesados_pdb=len(heavy), H_pdb=n_h,
               pesados_smiles=tpl.GetNumHeavyAtoms(), H_smiles=tpl_h,
               elementos_pdb=dict(pdb_el), elementos_smiles=dict(tpl_el),
               pesados_ok=len(heavy) == tpl.GetNumHeavyAtoms(), elementos_ok=dict(pdb_el) == dict(tpl_el),
               H_ok=(n_h == tpl_h) if n_h else None)
    try:
        mol = Chem.MolFromPDBBlock(''.join(heavy), removeHs=False, sanitize=False, proximityBonding=True)
        flat = Chem.Mol(tpl); Chem.RemoveStereochemistry(flat)
        new = AllChem.AssignBondOrdersFromTemplate(flat, mol)
        Chem.SanitizeMol(new)
        row.update(bond_orders_ok=True, formula_asignada=CalcMolFormula(new))
    except Exception as e:  # noqa: BLE001
        row.update(bond_orders_ok=False, error=str(e)[:200])
    return row


def build_mrtx_docking_complex():
    d = os.path.join(HERE, 'docking', 'MRTX1133'); os.makedirs(d, exist_ok=True)
    prot = [l for l in open(REDOCK + '7rpz_clean.pdb') if l.startswith('ATOM')]
    lig = Chem.SDMolSupplier(REDOCK + 'fred_docked.sdf', removeHs=False)[0]
    for i, a in enumerate(lig.GetAtoms()):
        inf = Chem.AtomPDBResidueInfo()
        inf.SetName(f' {a.GetSymbol():<2}{i + 1:<1}'[:4].ljust(4)); inf.SetResidueName('6IC'); inf.SetResidueNumber(202)
        inf.SetChainId('A'); inf.SetIsHeteroAtom(True); inf.SetSerialNumber(100000 + i)
        a.SetMonomerInfo(inf)
    blk = [l for l in Chem.MolToPDBBlock(lig).splitlines() if l.startswith(('HETATM', 'CONECT'))]
    path = os.path.join(d, 'complejo_origen.pdb')
    with open(path, 'w') as f:
        f.writelines(prot); f.write('TER\n'); f.write('\n'.join(blk) + '\nEND\n')
    return path


def build_mrtx_md_sdf(pdb):
    tpl = Chem.SDMolSupplier(MD_T155, removeHs=False)[0]
    heavy = ''.join(l for l in het_lines(pdb, '6IC') if elem(l) != 'H')
    mol = Chem.MolFromPDBBlock(heavy, removeHs=False, sanitize=False, proximityBonding=True)
    new = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(tpl), mol)
    Chem.SanitizeMol(new)
    out = os.path.join(HERE, 'md_rep', 'MRTX1133', 'ligand.sdf')
    w = Chem.SDWriter(out); w.write(new); w.close()
    return out, CalcMolFormula(new), Chem.GetFormalCharge(new)


if __name__ == '__main__':
    tab = smiles_table()
    json.dump(tab, open(os.path.join(HERE, 'ligandos.json'), 'w'), indent=1, ensure_ascii=False)
    cx = build_mrtx_docking_complex()
    print('complejo MRTX1133 docking ->', cx)
    rows = []
    for s in SYSTEMS:
        for struct, pdb in (('docking', os.path.join(HERE, 'docking', s, 'complejo_origen.pdb')),
                            ('md_rep', os.path.join(HERE, 'md_rep', s, 'frame_rep.pdb'))):
            if os.path.exists(pdb):
                rows.append(check(s, struct, pdb, RESN[s], tab[s]['smiles']))
    json.dump(rows, open(os.path.join(HERE, 'verificacion_ligandos.json'), 'w'), indent=1, ensure_ascii=False)
    sdf, f, q = build_mrtx_md_sdf(os.path.join(HERE, 'md_rep', 'MRTX1133', 'frame_rep.pdb'))
    print('SDF MRTX1133 MD ->', sdf, f, 'carga', q)
    for r in rows:
        print(f"{r['sistema']:9s} {r['estructura']:7s} {r['resname']:3s} pesados {r['pesados_pdb']}/{r['pesados_smiles']} "
              f"H {r['H_pdb']}/{r['H_smiles']} elementos_ok={r['elementos_ok']} bond_orders={r['bond_orders_ok']} "
              f"{r.get('formula_asignada', r.get('error'))}")
