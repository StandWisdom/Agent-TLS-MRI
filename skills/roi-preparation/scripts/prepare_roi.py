"""Align existing tumor masks and generate image/ROI overlays."""
import argparse
from pathlib import Path
import sys

import numpy as np
import SimpleITK as sitk

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'mri-preparation/scripts'))
from prepare_mri import digest, on_grid, read_inputs, resampling_policy, run_cases, same_grid


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
    allowed = resampling_policy(case)
    calculated_adc = (case.get('adc_source') == 'calculated' or
                      Path(case['image']).name in ('ADC_raw.nii.gz', 'ADC_clipped.nii.gz'))
    if calculated_adc and 'valid_mask' not in case:
        raise ValueError('Calculated ADC requires its ADC_valid.nii.gz as valid_mask')
    names = ['image', 'mask'] + (['valid_mask'] if 'valid_mask' in case else [])
    images, provenance = read_inputs(case, base, names)
    original_mask = sitk.GetArrayFromImage(images['mask'])
    if not np.all(np.isfinite(original_mask)) or np.any(original_mask < 0) or not np.all(original_mask == np.floor(original_mask)):
        raise ValueError('Mask must contain finite nonnegative integer labels')
    if not np.any(original_mask > 0):
        raise ValueError('Empty ROI')
    ref = images['image']
    mask_im, changed = on_grid(images['mask'], ref, allowed, True)
    mask = sitk.GetArrayFromImage(mask_im)
    if not np.any(mask > 0):
        raise ValueError('ROI has no overlap with reference')
    before = np.count_nonzero(original_mask) * np.prod(images['mask'].GetSpacing())
    after = np.count_nonzero(mask) * np.prod(ref.GetSpacing())
    if after < 0.9 * before:
        raise ValueError('ROI physical volume loss exceeds 10%; investigate coverage')
    array = sitk.GetArrayFromImage(ref)
    if 'valid_mask' in images:
        validity = images['valid_mask']
        if not same_grid(validity,ref):
            raise ValueError('ADC validity mask must share the reference grid')
        valid = sitk.GetArrayFromImage(validity)
        if not np.isin(valid,[0,1]).all():
            raise ValueError('ADC validity mask must be binary')
        if not np.all(valid[mask > 0] == 1):
            raise ValueError('ROI includes invalid ADC signal or missing diffusion coverage')
    if not np.isfinite(array[mask > 0]).all():
        raise ValueError('Nonfinite intensity inside ROI; check ADC valid signal and coverage')
    sitk.WriteImage(mask_im, str(dest / 'ROI_aligned.nii.gz'))
    preview = overlay(array, mask, dest / 'overlay.png')
    return dict(id=case['id'], status='REVIEW_REQUIRED', parameters=case, inputs=provenance,
                mask_resampled=changed, tumor_voxels=int(np.count_nonzero(mask)),
                validity_checked='valid_mask' in images,
                roi_volume_mm3=float(after), preview=preview,
                output_sha256={p.name:digest(p) for p in dest.iterdir() if p.is_file()})


def run(manifest, output):
    return run_cases(manifest, output, process, __file__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.manifest, args.output))
