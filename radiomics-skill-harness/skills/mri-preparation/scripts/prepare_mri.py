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


def overlay(array, mask, path):
    finite = array[np.isfinite(array)]
    if not finite.size:
        raise ValueError('No finite image intensities')
    lo, hi = np.percentile(finite, [1, 99])
    gray = (255 * np.clip((np.nan_to_num(array, nan=float(lo)) - lo) / max(hi-lo, 1e-12), 0, 1)).astype(np.uint8)
    occupied = np.flatnonzero(np.any(mask > 0, axis=(1, 2)))
    slices = [int(occupied[int(round(q * (len(occupied)-1)))]) for q in (0.25, 0.5, 0.75)]
    panels = []
    for z in slices:
        raw = np.repeat(gray[z, :, :, None], 3, axis=2)
        tint = raw.copy()
        roi = mask[z] > 0
        tint[roi] = (0.6 * raw[roi] + 0.4 * np.array([255, 0, 0])).astype(np.uint8)
        panels.append(np.concatenate([raw, tint], axis=1))
    sitk.WriteImage(sitk.GetImageFromArray(np.concatenate(panels, axis=0), isVector=True), str(path))
    return dict(slice_indices_z=slices, window=[float(lo), float(hi)], orientation='native grid')


def process(case, base, dest):
    if not case.get('acquisition_note', '').strip():
        raise ValueError('acquisition_note is required')
    allowed = case.get('allow_resample', False)
    if not isinstance(allowed, bool):
        raise ValueError('allow_resample must be boolean')
    if allowed and not case.get('resample_reason', '').strip():
        raise ValueError('resample_reason is required')
    mode = case['mode']
    if mode not in ('adc', 'overlay'):
        raise ValueError('mode must be adc or overlay')
    names = ['low', 'high', 'mask'] if mode == 'adc' else ['image', 'mask']
    images, provenance = {}, {}
    for name in names:
        path = (base / case[name]).resolve()
        im = sitk.ReadImage(str(path))
        if im.GetDimension() != 3 or im.GetNumberOfComponentsPerPixel() != 1:
            raise ValueError('Only scalar 3D images are supported')
        images[name] = im if name == 'mask' else sitk.Cast(im, sitk.sitkFloat32)
        provenance[name] = dict(path=str(path), sha256=digest(path), geometry=geometry(im))
    original_mask = sitk.GetArrayFromImage(images['mask'])
    if not np.all(np.isfinite(original_mask)) or np.any(original_mask < 0) or not np.all(original_mask == np.floor(original_mask)):
        raise ValueError('Mask must contain finite nonnegative integer labels')
    if not np.any(original_mask > 0):
        raise ValueError('Empty ROI')
    ref = images['high'] if mode == 'adc' else images['image']
    mask_im, mask_changed = on_grid(images['mask'], ref, allowed, True)
    mask = sitk.GetArrayFromImage(mask_im)
    if not np.any(mask > 0):
        raise ValueError('ROI has no overlap with reference')
    if int(np.count_nonzero(mask)) * np.prod(ref.GetSpacing()) < 0.9 * int(np.count_nonzero(original_mask)) * np.prod(images['mask'].GetSpacing()):
        raise ValueError('ROI physical volume loss exceeds 10%; investigate coverage')
    qc = dict(id=case['id'], status='REVIEW_REQUIRED', parameters=case,
              inputs=provenance, mask_resampled=mask_changed, tumor_voxels=int(np.count_nonzero(mask)))
    array = sitk.GetArrayFromImage(ref)
    if mode == 'adc':
        low, changed = on_grid(images['low'], ref, allowed)
        array, valid = adc_array(sitk.GetArrayFromImage(low), array, float(case['b_low']), float(case['b_high']))
        if not np.all(valid[mask > 0]):
            raise ValueError('ROI includes invalid/nonpositive DWI signal or missing low-b coverage')
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
    elif not np.all(np.isfinite(array[mask > 0])):
        raise ValueError('Nonfinite intensity inside ROI')
    sitk.WriteImage(mask_im, str(dest / 'ROI_aligned.nii.gz'))
    qc['preview'] = overlay(array, mask, dest / 'overlay.png')
    qc['output_sha256'] = {p.name: digest(p) for p in dest.iterdir() if p.is_file()}
    return qc


def run(manifest, output):
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
            qc = process(case, manifest.parent, dest)
        except Exception as exc:
            qc = dict(id=case['id'], status='FAILED', reason=str(exc), partial_outputs_not_approved=True)
        (dest / 'qc.json').write_text(json.dumps(qc, indent=2, allow_nan=False), encoding='utf-8')
        results.append(qc)
        print(case['id'], qc['status'])
    summary = dict(python=platform.python_version(), numpy=np.__version__,
                   SimpleITK=sitk.Version_VersionString(), script_sha256=digest(__file__),
                   manifest_sha256=digest(manifest), results=results)
    (output / 'run_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False), encoding='utf-8')
    return int(any(r['status'] == 'FAILED' for r in results))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.manifest, args.output))
