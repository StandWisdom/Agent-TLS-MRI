"""Synthetic-only regression tests. No clinical files are accessed."""
import importlib.util
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

import numpy as np
import pandas as pd
import SimpleITK as sitk

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT/'skills/tls-feature-reduction/scripts'
sys.path.insert(0, str(SCRIPTS))
from statistics_core import bh_fdr, pearson_p_values, vectorized_pearson


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


selection = module('tls_selection', SCRIPTS/'select_tls_features.py')
extraction = module('extraction', ROOT/'skills/radiomics-extraction/scripts/extract.py')


def arguments(root, stage='reduced'):
    return selection.parser().parse_args([
        '--features', str(root/'features'), '--tls', str(root/'tls.xlsx'),
        '--output', str(root/stage), '--stage', stage, '--bootstrap', '20'])


class StatisticsTests(unittest.TestCase):
    def test_statistics_against_scipy(self):
        from scipy.stats import false_discovery_control, pearsonr
        rng = np.random.default_rng(61)
        for n in (5,20,89):
            y = rng.normal(size=n)
            x = rng.normal(size=(n,8)) + y[:,None]*np.linspace(-2,2,8)
            r,p = vectorized_pearson(x,y)
            expected = [pearsonr(x[:,j],y) for j in range(x.shape[1])]
            np.testing.assert_allclose(r,[a.statistic for a in expected],atol=1e-12)
            np.testing.assert_allclose(p,[a.pvalue for a in expected],rtol=1e-10,atol=1e-14)
            np.testing.assert_allclose(bh_fdr(p),false_discovery_control(p,method='bh'))

    def test_cli_entrypoints(self):
        for script in (SCRIPTS/'select_tls_features.py',
                       ROOT/'skills/radiomics-extraction/scripts/extract.py'):
            result = subprocess.run([sys.executable, str(script), '--help'],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('usage:', result.stdout)

    def test_bh_known_vector_and_nan(self):
        np.testing.assert_allclose(bh_fdr([.01, .04, .03, .002, np.nan]),
                                   [.02, .04, .04, .008, 1.])

    def test_pearson_analytic_pvalues_and_constants(self):
        # With n=4 the exact two-sided p-value is 1-|r|.
        np.testing.assert_allclose(pearson_p_values(np.array([0., .2, -.5, 1.]), 4),
                                   [1., .8, .5, 0.], atol=1e-12)
        x = np.arange(10.)
        r, p = vectorized_pearson(np.column_stack([x, -x, np.ones(10)]), x)
        np.testing.assert_allclose(r[:2], [1., -1.])
        np.testing.assert_allclose(p[:2], [0., 0.])
        self.assertTrue(np.isnan(r[2]))
        r, p = vectorized_pearson(x[:, None], np.ones(10))
        self.assertTrue(np.isnan(r).all())
        self.assertEqual(p[0], 1.)

    def test_strict_legacy_boundary_and_reduction_modes(self):
        args = arguments(Path('.'))
        r = np.array([.1, .3, .3, .3])
        p = np.array([.01, .05, .01, .01])
        q = np.array([.01, .01, .01, .2])
        freq = np.array([1., 1., .2, .9])
        legacy, reduced = selection.gates(r, p, q, freq, np.ones(4), args)
        np.testing.assert_array_equal(legacy, [False, False, True, True])
        np.testing.assert_array_equal(reduced, legacy)
        args.combine = 'and'
        _, reduced = selection.gates(r, p, q, freq, np.ones(4), args)
        self.assertFalse(reduced.any())

    def test_seeded_bootstrap(self):
        args = arguments(Path('.'))
        rng = np.random.default_rng(5)
        y = np.arange(30.)
        x = np.column_stack([y, rng.normal(size=30)])
        r, _ = vectorized_pearson(x, y)
        a = selection.stability(x, y, r, args, 12)
        b = selection.stability(x, y, r, args, 12)
        np.testing.assert_array_equal(a[0], b[0])
        self.assertEqual(a[0][0], 1.)
        self.assertEqual(a[2], 20)


class SelectionIntegrationTests(unittest.TestCase):
    def create_data(self, root):
        (root/'features').mkdir()
        rng = np.random.default_rng(10)
        ids = [f'synthetic{i:03}' for i in range(40)]
        y = np.arange(40.)
        tls = pd.DataFrame({'ID': ids, **{m: y for m in selection.MARKERS}})
        tls.to_excel(root/'tls.xlsx', index=False)
        for seq in selection.SEQUENCES:
            pd.DataFrame({'ID': ids, 'original_firstorder_Mean': y,
                          'wavelet-LLH_firstorder_Mean': rng.normal(size=40),
                          'original_firstorder_Minimum': np.ones(40)}).to_excel(
                              root/'features'/f'{seq}_features.xlsx', index=False)

    def test_full_workbooks_subset_and_original_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.create_data(root)
            source_hash = selection.sha(root/'tls.xlsx')
            for stage in ('legacy', 'reduced'):
                self.assertEqual(selection.run(arguments(root, stage)), 0)
            self.assertEqual(source_hash, selection.sha(root/'tls.xlsx'))
            for seq in selection.SEQUENCES:
                old = pd.read_excel(root/'legacy'/f'{seq}_legacy_features.xlsx')
                path = root/'reduced'/f'{seq}_reduced_features.xlsx'
                new = pd.read_excel(path)
                self.assertTrue(set(new).issubset(old.columns))
                self.assertIn('original_firstorder_Mean', new)
                self.assertNotIn('original_firstorder_Minimum', new)
                self.assertEqual(len(new), 40)
                tests = pd.read_excel(path, sheet_name='AllTests')
                self.assertEqual(len(tests), 12*3)
                self.assertTrue((tests.groupby('TLSMarker').size() == 3).all())

    def test_missing_values_fail_closed_and_duplicate_ids_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.create_data(root)
            tls = pd.read_excel(root/'tls.xlsx')
            tls.loc[0, selection.MARKERS[0]] = np.nan
            tls.to_excel(root/'tls.xlsx', index=False)
            self.assertEqual(selection.run(arguments(root)), 1)
            self.assertEqual(len(list((root/'reduced').glob('*.xlsx'))), 0)
            tls.loc[1, 'ID'] = tls.loc[0, 'ID']
            tls.to_excel(root/'tls.xlsx', index=False)
            with self.assertRaises(ValueError):
                selection.table(root/'tls.xlsx')


class ExtractionIntegrationTests(unittest.TestCase):
    def create_data(self, root):
        rng = np.random.default_rng(22)
        image = sitk.GetImageFromArray(rng.normal(100, 15, (16, 16, 16)).astype('float32'))
        roi = np.zeros((16, 16, 16), dtype='uint8')
        roi[3:13, 3:13, 3:13] = 1
        mask = sitk.GetImageFromArray(roi)
        sitk.WriteImage(image, str(root/'image.mha'))
        sitk.WriteImage(mask, str(root/'roi.mha'))
        manifest = {'cases': [{'id': 'synthetic001', 'reviewed': True,
                               'sequences': {s: {'image': 'image.mha', 'mask': 'roi.mha'}
                                             for s in extraction.SEQUENCES}}]}
        (root/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')

    def test_legacy_yaml_four_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.create_data(root)
            yaml = ROOT/'skills/radiomics-extraction/assets/legacy_mr.yaml'
            self.assertEqual(extraction.run(root/'manifest.json', yaml, root/'out'), 0)
            for seq in extraction.SEQUENCES:
                frame = pd.read_excel(root/'out'/f'{seq}_features.xlsx')
                self.assertEqual(len(frame), 1)
                self.assertTrue(any(c.startswith('original_') for c in frame))
                self.assertTrue(any(c.startswith('wavelet-') for c in frame))
                self.assertTrue(any(c.startswith('log-') for c in frame))
                self.assertFalse(any(c.startswith('diagnostics_') for c in frame))
            report = json.loads((root/'out/extraction_qc.json').read_text())
            self.assertEqual(report['effective_settings']['geometryTolerance'], 1e-5)

    def test_geometry_mismatch_does_not_auto_correct(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.create_data(root)
            mask = sitk.ReadImage(str(root/'roi.mha'))
            mask.SetOrigin((2., 0., 0.))
            sitk.WriteImage(mask, str(root/'roi.mha'))
            yaml = ROOT/'skills/radiomics-extraction/assets/legacy_mr.yaml'
            self.assertEqual(extraction.run(root/'manifest.json', yaml, root/'out'), 1)
            report = json.loads((root/'out/extraction_qc.json').read_text())
            self.assertTrue(all('Geometry mismatch' in r['reason'] for r in report['cases']))


if __name__ == '__main__':
    unittest.main()
