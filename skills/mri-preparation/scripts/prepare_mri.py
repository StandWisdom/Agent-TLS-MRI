"""Explicit MRI preparation with local provenance; no automatic registration."""
import argparse
import hashlib
import json
import platform
import re
from pathlib import Path

import numpy as np
import SimpleITK as sitk


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def geometry(im):
    return dict(size=im.GetSize(), spacing=im.GetSpacing(),
                origin=im.GetOrigin(), direction=im.GetDirection())


def same_grid(a, b):
    return a.GetSize() == b.GetSize() and all(
        np.allclose(getattr(a, key)(), getattr(b, key)(), atol=1e-5, rtol=0)
        for key in ('GetSpacing', 'GetOrigin', 'GetDirection'))


def on_grid(im, ref, allowed, label=False):
    if same_grid(im, ref):
        return im, False
    if not allowed:
        raise ValueError('Grid mismatch: confirm shared physical coordinates before resampling')
    result = sitk.Resample(im, ref, sitk.Transform(3, sitk.sitkIdentity),
                          sitk.sitkNearestNeighbor if label else sitk.sitkLinear,
                          0, im.GetPixelID())
    return result, True


def adc_array(low, high, b_low, b_high):
    if not (np.isfinite(b_low) and np.isfinite(b_high) and 0 <= b_low < b_high):
        raise ValueError('Require finite b_high > b_low >= 0 in s/mm2')
    valid = np.isfinite(low) & np.isfinite(high) & (low > 0) & (high > 0)
    adc = np.full(low.shape, np.nan, dtype=np.float32)
    adc[valid] = (np.log(low[valid].astype(float)) - np.log(high[valid].astype(float))) / (b_high - b_low)
    return adc, valid


def save_volume(array, reference, path):
    im = sitk.GetImageFromArray(array)
    im.CopyInformation(reference)
    sitk.WriteImage(im, str(path))


def resampling_policy(case):
    if not case.get('acquisition_note', '').strip():
        raise ValueError('acquisition_note is required')
    allowed = case.get('allow_resample', False)
    if not isinstance(allowed, bool):
        raise ValueError('allow_resample must be boolean')
    if allowed and not case.get('resample_reason', '').strip():
        raise ValueError('resample_reason is required')
    return allowed


def read_inputs(case, base, names):
    images, provenance = {}, {}
    for name in names:
        path = (base / case[name]).resolve()
        im = sitk.ReadImage(str(path))
        if im.GetDimension() != 3 or im.GetNumberOfComponentsPerPixel() != 1:
            raise ValueError('Only scalar 3D images are supported')
        images[name] = im if name == 'mask' else sitk.Cast(im, sitk.sitkFloat32)
        provenance[name] = dict(path=str(path), sha256=digest(path), geometry=geometry(im))
    return images, provenance


def process(case, base, dest):
    allowed = resampling_policy(case)
    if 'mask' in case:
        raise ValueError('Move mask and overlay processing to the roi-preparation skill')
    mode = case['mode']
    if mode not in ('adc', 'image'):
        raise ValueError('MRI mode must be adc or image; overlays belong to roi-preparation')
    names = ['low', 'high'] if mode == 'adc' else ['image']
    images, provenance = read_inputs(case, base, names)
    qc = dict(id=case['id'], status='REVIEW_REQUIRED', parameters=case, inputs=provenance)
    if mode == 'adc':
        ref = images['high']
        low, changed = on_grid(images['low'], ref, allowed)
        array, valid = adc_array(sitk.GetArrayFromImage(low), sitk.GetArrayFromImage(ref),
                                 float(case['b_low']), float(case['b_high']))
        if not np.any(valid):
            raise ValueError('No valid positive paired diffusion signal')
        qc.update(low_resampled=changed, adc_unit='mm2/s', valid_fraction=float(valid.mean()),
                  negative_fraction=float(np.mean(array[valid] < 0)))
        save_volume(array, ref, dest / 'ADC_raw.nii.gz')
        save_volume(valid.astype(np.uint8), ref, dest / 'ADC_valid.nii.gz')
        if 'adc_clip_max' in case:
            maximum = float(case['adc_clip_max'])
            if not np.isfinite(maximum) or maximum <= 0:
                raise ValueError('adc_clip_max must be finite and positive')
            save_volume(np.clip(array, 0, maximum), ref, dest / 'ADC_clipped.nii.gz')
            qc['above_clip_fraction'] = float(np.mean(array[valid] > maximum))
    else:
        array = sitk.GetArrayFromImage(images['image'])
        valid = np.isfinite(array)
        if not valid.any():
            raise ValueError('No finite image intensities')
        qc.update(finite_fraction=float(valid.mean()),
                  intensity_min=float(array[valid].min()), intensity_max=float(array[valid].max()))
    qc['output_sha256'] = {p.name: digest(p) for p in dest.iterdir() if p.is_file()}
    return qc


def run_cases(manifest, output, processor, script_path):
    manifest, output = Path(manifest).resolve(), Path(output).resolve()
    spec = json.loads(manifest.read_text(encoding='utf-8'))
    cases = spec['cases']
    if not cases:
        raise ValueError('Empty case list')
    ids = [c['id'] for c in cases]
    if len(set(ids)) != len(ids) or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', i) for i in ids):
        raise ValueError('IDs must be unique path-safe pseudonyms')
    output.mkdir(parents=True, exist_ok=False)
    results = []
    for case in cases:
        dest = output / case['id']
        dest.mkdir()
        try:
            qc = processor(case, manifest.parent, dest)
        except Exception as exc:
            qc = dict(id=case['id'], status='FAILED', reason=str(exc), partial_outputs_not_approved=True)
        (dest / 'qc.json').write_text(json.dumps(qc, indent=2, allow_nan=False), encoding='utf-8')
        results.append(qc)
        print(case['id'], qc['status'])
    summary = dict(python=platform.python_version(), numpy=np.__version__,
                   SimpleITK=sitk.Version_VersionString(), script_sha256=digest(script_path), utilities_sha256=digest(__file__),
                   manifest_sha256=digest(manifest), results=results)
    (output / 'run_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False), encoding='utf-8')
    return int(any(r['status'] == 'FAILED' for r in results))


def run(manifest, output):
    return run_cases(manifest, output, process, __file__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.manifest, args.output))
