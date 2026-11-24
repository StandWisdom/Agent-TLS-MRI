# MRI input and output contract

The manifest contains `cases` with unique path-safe IDs. Paths are relative to the manifest location or absolute local paths. Each case has `mode`, `sequence`, and a nonempty `acquisition_note`.

| Mode | Required inputs | Operation |
|---|---|---|
| adc | low, high, b_low, b_high | Compute two-point ADC |
| image | image | Inspect geometry, finite intensities and intensity range |

Supply b-values in s/mm2 with b_high > b_low >= 0. Check same patient/time point, signal scaling, reconstruction, diffusion direction averaging, motion and field-of-view overlap before execution. Supply confirmed sequence identities; DCE or T1+C filenames alone do not establish enhancement phase.

Optional ADC fields: `allow_resample` (default false), `resample_reason` (required if true), `adc_clip_max` (positive upper limit in mm2/s). Resampling uses an identity transform in physical coordinates. Raw ADC retains invalid voxels as NaN; cases with no valid paired signal fail.

Outputs:
- ADC mode: `ADC_raw.nii.gz`, `ADC_valid.nii.gz`, optional `ADC_clipped.nii.gz`.
- Image mode: metadata QC referencing the original image, without rewriting the image.
- Per-case `qc.json`: parameters, geometry, source hashes, signal statistics and review state.
- `run_summary.json`: environment, script/manifest hashes and case outcomes.

Next, use the original MRI or chosen ADC derivative as the `image` input to ROI preparation. For calculated ADC, also pass `ADC_valid.nii.gz` as `valid_mask` and set `adc_source: calculated`. The ROI step checks tumor coverage, signal validity and finite intensities. Record approval separately in a local review record. Invalid ADC outside the tumor requires an explicit finite-image policy before radiomics extraction.

## Updating from the combined skill

Remove `mask` from MRI manifests. Move former `mode: overlay` entries to the ROI manifest with `image` and `mask`. Call the two runners separately. The MRI runner rejects legacy combined entries to avoid silently skipping ROI checks.

