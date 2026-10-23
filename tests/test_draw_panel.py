import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import COMPLEX, DATOS, MINI, ROOT, SMILES_MOL3, XML  # noqa: E402
from deku2d import Style, legacy_style, render_panel, render_single, save_figure, scene_from_spec  # noqa: E402
from deku2d.pipeline import SystemSpec  # noqa: E402

SPEC = SystemSpec(ligand_resname="UNK", plip_xml=XML, smiles=SMILES_MOL3, name="mol3", title="mol3")


class TestDraw(unittest.TestCase):
    def test_draw_png_not_empty(self):
        st = Style()
        sc = scene_from_spec(SPEC, st)
        self.assertEqual(len(sc.keys), 10)
        self.assertEqual(len(sc.contacts), 13)
        with tempfile.TemporaryDirectory() as d:
            paths = save_figure(render_single(sc, st), os.path.join(d, "x.png"), None, 100)
            self.assertTrue(os.path.getsize(paths[0]) > 20000)

    def test_multiple_formats(self):
        st = Style()
        sc = scene_from_spec(SPEC, st)
        with tempfile.TemporaryDirectory() as d:
            paths = save_figure(render_single(sc, st), os.path.join(d, "x.png"), ["png", "svg", "pdf"], 72)
            self.assertEqual([os.path.splitext(p)[1] for p in paths], [".png", ".svg", ".pdf"])
            self.assertTrue(all(os.path.getsize(p) > 1000 for p in paths))

    def test_no_box_overlap_and_distance_labels(self):
        st = Style()
        sc = scene_from_spec(SPEC, st)
        keys = sc.keys
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = sc.pos[keys[i]], sc.pos[keys[j]]
                over = abs(a[0] - b[0]) < sc.hw[keys[i]] + sc.hw[keys[j]] and abs(a[1] - b[1]) < 2 * sc.hh
                self.assertFalse(over, (keys[i], keys[j]))
        self.assertTrue(all(p.label_pos is not None for p in sc.polar if p.kind in ("hbond", "salt")))

    def test_no_dist_and_spanish(self):
        st = Style(dist=False, lang="es")
        sc = scene_from_spec(SPEC, st)
        self.assertTrue(all(p.label_pos is None for p in sc.polar))
        with tempfile.TemporaryDirectory() as d:
            paths = save_figure(render_single(sc, st), os.path.join(d, "es.png"), None, 72)
            self.assertTrue(os.path.getsize(paths[0]) > 20000)

    def test_seed_is_reproducible(self):
        a = scene_from_spec(SPEC, Style(seed=5))
        b = scene_from_spec(SPEC, Style(seed=5))
        self.assertEqual({k: tuple(v) for k, v in a.pos.items()}, {k: tuple(v) for k, v in b.pos.items()})

    def test_legacy_matches_reference_cost(self):
        sc = scene_from_spec(SPEC, legacy_style())
        self.assertAlmostEqual(sc.cost, 15.3, delta=0.1)       # coste del script original para mol3

    def test_pdb_only_runs_plip_if_available(self):
        from deku2d.config import DEFAULT_PLIP_EXE
        if not os.path.exists(DEFAULT_PLIP_EXE):
            self.skipTest("PLIP no disponible")
        with tempfile.TemporaryDirectory() as d:
            spec = SystemSpec(ligand_resname="UNK", pdb=COMPLEX, smiles=SMILES_MOL3, workdir=os.path.join(d, "plip"), name="e2e")
            sc = scene_from_spec(spec, Style(restarts=3, iterations=300))
            self.assertGreaterEqual(len(sc.contacts), 10)


class TestPanel(unittest.TestCase):
    def test_panel_two_systems(self):
        st = Style(restarts=4, iterations=400)
        s1 = scene_from_spec(SPEC, st)
        s2 = scene_from_spec(SystemSpec(ligand_resname="UNK", plip_xml=XML, smiles=SMILES_MOL3, name="b", title="otro"), st)
        fig = render_panel([s1, s2], st, ncols=2)
        w, h = fig.get_size_inches()
        self.assertGreater(w, h * 0.8)
        with tempfile.TemporaryDirectory() as d:
            paths = save_figure(fig, os.path.join(d, "panel.png"), None, 80)
            self.assertTrue(os.path.getsize(paths[0]) > 30000)

    def test_panel_single_column_three_panels(self):
        st = Style(restarts=3, iterations=300)
        sc = [scene_from_spec(SystemSpec(ligand_resname="UNK", plip_xml=XML, smiles=SMILES_MOL3, name=str(i), title=f"s{i}"), st)
              for i in range(3)]
        fig = render_panel(sc, st, ncols=1)
        w, h = fig.get_size_inches()
        self.assertGreater(h, w)


class TestCLI(unittest.TestCase):
    def run_cli(self, *args):
        env = dict(os.environ, PYTHONPATH=ROOT)
        return subprocess.run([sys.executable, "-m", "deku2d", *args], capture_output=True, text=True, env=env, cwd=tempfile.gettempdir())

    def test_draw_cli(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "f.png")
            r = self.run_cli("draw", "--plip-xml", XML, "--ligand-resname", "UNK", "--smiles", SMILES_MOL3, "-o", out,
                             "--dpi", "80", "--restarts", "3", "--iterations", "300", "--format", "png,svg", "-q")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue(os.path.getsize(out) > 20000 and os.path.exists(os.path.join(d, "f.svg")))

    def test_draw_cli_error_is_readable(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_cli("draw", "--plip-xml", XML, "--ligand-resname", "XXX", "-o", os.path.join(d, "f.png"))
            self.assertEqual(r.returncode, 2)
            self.assertIn("Error:", r.stderr)

    def test_panel_and_batch_cli(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = {"ncols": 2, "style": {"restarts": 3, "iterations": 300, "dpi": 80},
                   "systems": [{"title": "A", "name": "a", "plip_xml": XML, "ligand_resname": "UNK", "smiles": SMILES_MOL3},
                               {"title": "B", "name": "b", "plip_xml": XML, "ligand_resname": "UNK", "smiles": SMILES_MOL3}],
                   "panel": {"output": "panel", "ncols": 2}}
            cp = os.path.join(d, "cfg.json")
            with open(cp, "w") as fh:
                json.dump(cfg, fh)
            r = self.run_cli("panel", "--config", cp, "-o", os.path.join(d, "p.png"), "-q")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue(os.path.getsize(os.path.join(d, "p.png")) > 30000)
            r = self.run_cli("batch", "--config", cp, "--outdir", os.path.join(d, "out"), "-q")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue(all(os.path.exists(os.path.join(d, "out", f)) for f in ("a.png", "b.png", "panel.png")))

    def test_unknown_style_key(self):
        with tempfile.TemporaryDirectory() as d:
            cp = os.path.join(d, "cfg.json")
            with open(cp, "w") as fh:
                json.dump({"style": {"colorin": 1}, "systems": []}, fh)
            r = self.run_cli("panel", "--config", cp, "-o", os.path.join(d, "p.png"))
            self.assertEqual(r.returncode, 2)
            self.assertIn("desconocidos", r.stderr)


if __name__ == "__main__":
    unittest.main()
