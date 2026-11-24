import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import SimpleITK as sitk

SCRIPT = Path(__file__).resolve().parents[1] / 'skills/mri-preparation/scripts/prepare_mri.py'
spec = importlib.util.spec_from_file_location('preparation', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PreparationTests(unittest.TestCase):
    def test_known_adc_and_invalid_signal(self):
        low = np.array([100., 100., 0., 100.])
        high = np.array([100 * np.exp(-0.8), 120., 1., np.nan])
        adc, valid = module.adc_array(low, high, 0, 800)
        self.assertAlmostEqual(float(adc[0]), 0.001, places=8)
        self.assertLess(adc[1], 0)
        self.assertTrue(np.isnan(adc[2:]).all())
        self.assertEqual(valid.tolist(), [True, True, False, False])
        with self.assertRaises(ValueError):
            module.adc_array(low, high, 800, 0)

    def test_geometry_checks_origin_not_just_shape(self):
        a = sitk.Image([8, 8, 8], sitk.sitkUInt8)
        b = sitk.Image(a)
        b.SetOrigin((20, 0, 0))
        self.assertFalse(module.same_grid(a, b))
        with self.assertRaises(ValueError):
            module.on_grid(a, b, False, True)

    def test_resampling_preserves_labels(self):
        a = np.zeros((4, 8, 8), dtype=np.uint8)
        a[:, 2:6, 2:6] = 2
        mask = sitk.GetImageFromArray(a)
        ref = sitk.Image([16, 16, 4], sitk.sitkFloat32)
        ref.SetSpacing((0.5, 0.5, 1))
        aligned, changed = module.on_grid(mask, ref, True, True)
        self.assertTrue(changed)
        self.assertEqual(set(np.unique(sitk.GetArrayFromImage(aligned))), {0, 2})

    def test_end_to_end_synthetic_and_failure_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            low = np.full((6, 12, 12), 100, dtype=np.float32)
            roi = np.zeros(low.shape, dtype=np.uint8)
            roi[1:5, 3:9, 3:9] = 1
            for name, data in [('low', low), ('high', low*np.exp(-0.8)), ('mask', roi)]:
                sitk.WriteImage(sitk.GetImageFromArray(data), str(root / (name+'.nii.gz')))
            case = dict(id='synthetic', mode='adc', low='low.nii.gz', high='high.nii.gz',
                        b_low=0, b_high=800, acquisition_note='synthetic phantom')
            manifest = root/'manifest.json'
            manifest.write_text(json.dumps(dict(cases=[case, dict(case, id='missing', low='missing.nii.gz')])) )
            self.assertEqual(module.run(manifest, root/'run'), 1)
            qc = json.loads((root/'run/synthetic/qc.json').read_text())
            self.assertEqual(qc['status'], 'REVIEW_REQUIRED')
            self.assertFalse((root/'run/synthetic/overlay.png').exists())
            self.assertFalse((root/'run/synthetic/ROI_aligned.nii.gz').exists())
            raw = sitk.GetArrayFromImage(sitk.ReadImage(str(root/'run/synthetic/ADC_raw.nii.gz')))
            np.testing.assert_allclose(raw, 0.001, rtol=1e-5)
            with self.assertRaises(FileExistsError):
                module.run(manifest, root/'run')

    def test_existing_image_needs_no_roi_and_combined_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sitk.WriteImage(sitk.Image([8,8,8],sitk.sitkFloat32),str(root/'image.nii.gz'))
            case = dict(id='synthetic',mode='image',image='image.nii.gz',acquisition_note='synthetic')
            path = root/'manifest.json'
            path.write_text(json.dumps({'cases':[case]}))
            self.assertEqual(module.run(path,root/'image_run'),0)
            case['mask'] = 'unused.nii.gz'
            path.write_text(json.dumps({'cases':[case]}))
            self.assertEqual(module.run(path,root/'old_format'),1)
            qc = json.loads((root/'old_format/synthetic/qc.json').read_text())
            self.assertIn('roi-preparation',qc['reason'])


if __name__ == '__main__':
    unittest.main()
