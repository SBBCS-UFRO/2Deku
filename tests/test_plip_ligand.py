import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collections import Counter  # noqa: E402

from common import COMPLEX, FIXED, SMILES_MOL3, XML  # noqa: E402
from deku2d import ligand, plip_io  # noqa: E402


class TestPlipIO(unittest.TestCase):
    def test_parse_counts(self):
        site = plip_io.parse_report(XML, "UNK")
        self.assertEqual((site.hetid, site.chain, site.position), ("UNK", "Z", "900"))
        self.assertEqual(Counter(i.kind for i in site.interactions), {"hbond": 8, "hydrophobic": 5})

    def test_hbond_direction_and_fields(self):
        site = plip_io.parse_report(XML, "UNK")
        hb = [i for i in site.interactions if i.kind == "hbond"]
        self.assertTrue(any(i.prot_donor for i in hb) and any(not i.prot_donor for i in hb))
        for i in site.interactions:
            self.assertTrue(i.res.isupper() and i.num > 0 and i.chain == "A" and i.dist > 0 and i.serials)

    def test_missing_ligand_gives_clear_error(self):
        with self.assertRaises(ValueError) as cm:
            plip_io.parse_report(XML, "XXX")
        self.assertIn("UNK", str(cm.exception))

    def test_list_sites(self):
        self.assertIn("UNK", [s.hetid for s in plip_io.list_sites(XML)])

    def test_find_fixed_pdb(self):
        self.assertEqual(os.path.basename(plip_io.find_fixed_pdb(os.path.dirname(XML))), os.path.basename(FIXED))

    def test_missing_plip_exe(self):
        with self.assertRaises(FileNotFoundError):
            plip_io.run_plip(COMPLEX, os.path.join(os.path.dirname(XML), "_no_existe"), exe="/no/existe/plip", reuse=False)


class TestLigand(unittest.TestCase):
    def test_atom_mapping_with_smiles(self):
        lig = ligand.build_ligand(FIXED, "UNK", SMILES_MOL3)
        self.assertEqual(lig.mol.GetNumAtoms(), 27)
        self.assertEqual(len(set(lig.serial_of.values())), 27)
        site = plip_io.parse_report(XML, "UNK")
        s2i = lig.serial_to_idx
        for it in site.interactions:
            self.assertTrue(all(s in s2i for s in it.serials), it)

    def test_smiles_heavy_atom_mismatch(self):
        with self.assertRaises(ligand.LigandError):
            ligand.build_ligand(FIXED, "UNK", "CCO")

    def test_invalid_smiles(self):
        with self.assertRaises(ligand.LigandError):
            ligand.build_ligand(FIXED, "UNK", "esto-no-es-smiles")

    def test_2d_coordinates_present(self):
        lig = ligand.build_ligand(FIXED, "UNK", SMILES_MOL3)
        conf = lig.mol.GetConformer()
        self.assertFalse(conf.Is3D())

    def test_pdb_only_fallback_maps_atoms(self):
        lig = ligand.build_ligand(FIXED, "UNK")
        self.assertEqual(lig.mol.GetNumAtoms(), 27)
        self.assertEqual(len(lig.serial_of), 27)

    def test_unknown_resname(self):
        with self.assertRaises(ligand.LigandError):
            ligand.build_ligand(FIXED, "ZZZ", SMILES_MOL3)

    def test_sdf_route_matches_smiles_route(self):
        import tempfile
        from rdkit import Chem
        from rdkit.Chem import AllChem
        lines, _ = ligand.read_pdb_ligand(FIXED, "UNK")
        m = Chem.MolFromPDBBlock("".join(lines), removeHs=False, sanitize=False, proximityBonding=True)
        flat = Chem.MolFromSmiles(SMILES_MOL3)
        Chem.RemoveStereochemistry(flat)
        m = AllChem.AssignBondOrdersFromTemplate(flat, m)
        Chem.SanitizeMol(m)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "lig.sdf")
            w = Chem.SDWriter(path); w.write(m); w.close()
            a = ligand.build_ligand(FIXED, "UNK", sdf=path)
        b = ligand.build_ligand(FIXED, "UNK", SMILES_MOL3)
        self.assertEqual(sorted(a.serial_of.values()), sorted(b.serial_of.values()))
        self.assertEqual(a.mol.GetNumBonds(), b.mol.GetNumBonds())
        self.assertEqual(a.warnings, [])

    def test_two_copies_warns(self):
        import tempfile
        lines, _ = ligand.read_pdb_ligand(FIXED, "UNK")
        dup = [l[:21] + "B" + l[22:] for l in lines]            # segunda copia en la cadena B
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "dos.pdb")
            with open(path, "w") as fh:
                fh.write("".join(lines) + "".join(dup))
            lig = ligand.build_ligand(path, "UNK", SMILES_MOL3, chain="B")
            self.assertEqual(lig.mol.GetNumAtoms(), 27)
            _, warns = ligand.read_pdb_ligand(path, "UNK")
            self.assertTrue(warns and "copias" in warns[0])


if __name__ == "__main__":
    unittest.main()
