"""Separate MRI-to-ROI handoff and overlay checks using synthetic volumes."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import SimpleITK as sitk

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'skills/roi-preparation/scripts/prepare_roi.py'
spec = importlib.util.spec_from_file_location('roi_preparation',SCRIPT)
roi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(roi)
import prepare_mri as mri


class ROIPreparationTests(unittest.TestCase):
    def setup_images(self,root):
        low = np.full((6,12,12),100,dtype=np.float32)
        mask = np.zeros(low.shape,dtype=np.uint8)
        mask[1:5,3:9,3:9] = 1
        for name,data in [('low',low),('high',low*np.exp(-.8)),('mask',mask)]:
            sitk.WriteImage(sitk.GetImageFromArray(data),str(root/(name+'.nii.gz')))
        return dict(id='synthetic',image='high.nii.gz',mask='mask.nii.gz',acquisition_note='synthetic phantom')

    def test_mri_adc_to_roi_overlay_handoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = self.setup_images(root)
            mri_manifest = root/'mri.json'
            mri_manifest.write_text(json.dumps({'cases':[dict(id='synthetic',mode='adc',
                low='low.nii.gz',high='high.nii.gz',b_low=0,b_high=800,acquisition_note='synthetic')]}))
            self.assertEqual(mri.run(mri_manifest,root/'mri'),0)
            case['image'] = 'mri/synthetic/ADC_raw.nii.gz'
            case['valid_mask'] = 'mri/synthetic/ADC_valid.nii.gz'
            manifest = root/'roi.json'
            manifest.write_text(json.dumps({'cases':[case]}))
            self.assertEqual(roi.run(manifest,root/'roi'),0)
            image = sitk.ReadImage(str(root/'roi/synthetic/overlay.png'))
            self.assertEqual(image.GetNumberOfComponentsPerPixel(),3)
            self.assertEqual(image.GetSize(),(24,36))
            rgb = sitk.GetArrayFromImage(image)
            self.assertGreater(int((rgb[:,:,0]>rgb[:,:,1]).sum()),0)
            qc = json.loads((root/'roi/synthetic/qc.json').read_text())
            self.assertEqual(qc['status'],'REVIEW_REQUIRED')
            self.assertEqual(qc['tumor_voxels'],144)
            self.assertFalse((root/'roi/synthetic/ADC_raw.nii.gz').exists())

    def test_roi_rejects_invalid_adc_inside_tumor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = self.setup_images(root)
            im = sitk.ReadImage(str(root/'high.nii.gz'))
            data = sitk.GetArrayFromImage(im)
            data[2,4,4] = np.nan
            sitk.WriteImage(sitk.GetImageFromArray(data),str(root/'nonfinite.mha'))
            case['image'] = 'nonfinite.mha'
            manifest = root/'roi.json'
            manifest.write_text(json.dumps({'cases':[case]}))
            self.assertEqual(roi.run(manifest,root/'roi'),1)
            qc = json.loads((root/'roi/synthetic/qc.json').read_text())
            self.assertIn('Nonfinite intensity inside ROI',qc['reason'])

    def test_calculated_adc_validity_survives_file_handoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = self.setup_images(root)
            low = sitk.GetArrayFromImage(sitk.ReadImage(str(root/'low.nii.gz')))
            low[2,4,4] = 0
            sitk.WriteImage(sitk.GetImageFromArray(low),str(root/'low.nii.gz'))
            manifest = root/'mri.json'
            manifest.write_text(json.dumps({'cases':[dict(id='synthetic',mode='adc',
                low='low.nii.gz',high='high.nii.gz',b_low=0,b_high=800,acquisition_note='synthetic')]}))
            self.assertEqual(mri.run(manifest,root/'mri'),0)
            case['image'] = 'mri/synthetic/ADC_raw.nii.gz'
            manifest.write_text(json.dumps({'cases':[case]}))
            self.assertEqual(roi.run(manifest,root/'missing_valid'),1)
            case['valid_mask'] = 'mri/synthetic/ADC_valid.nii.gz'
            manifest.write_text(json.dumps({'cases':[case]}))
            self.assertEqual(roi.run(manifest,root/'roi'),1)
            qc = json.loads((root/'roi/synthetic/qc.json').read_text())
            self.assertIn('invalid ADC signal',qc['reason'])

    def test_empty_roi_and_unapproved_geometry_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = self.setup_images(root)
            manifest = root/'roi.json'
            manifest.write_text(json.dumps({'cases':[case]}))
            mask = sitk.ReadImage(str(root/'mask.nii.gz'))
            mask.SetOrigin((2.,0.,0.))
            sitk.WriteImage(mask,str(root/'mask.nii.gz'))
            self.assertEqual(roi.run(manifest,root/'mismatch'),1)
            empty = sitk.Image(mask.GetSize(),sitk.sitkUInt8)
            sitk.WriteImage(empty,str(root/'mask.nii.gz'))
            self.assertEqual(roi.run(manifest,root/'empty'),1)

    def test_separate_cli_entrypoints(self):
        for path in (SCRIPT,ROOT/'skills/mri-preparation/scripts/prepare_mri.py'):
            result = subprocess.run([sys.executable,str(path),'--help'],capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('--output',result.stdout)


if __name__=='__main__':
    unittest.main()
