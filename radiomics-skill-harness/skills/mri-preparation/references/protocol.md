# Input and output contract

The manifest contains a `cases` array with unique path-safe IDs. Paths are relative to the manifest location or absolute local paths. Each case has `mode` (`adc` or `overlay`), `mask`, and a nonempty `acquisition_note`. For ADC supply `low`, `high`, `b_low`, and `b_high`; for overlay supply `image`. Supply `sequence` to describe the confirmed acquisition, not an inferred filename mapping.

Required acquisition checks (performed before running): same patient/time point; actual b-values and units; matched signal scaling and reconstruction; diffusion direction averaging compatibility; motion/anatomical alignment; field-of-view overlap. ADC is computed as ln(S_low/S_high)/(b_high-b_low), in mm2/s.

Optional fields: `allow_resample` (default false), `resample_reason` (required if true), `adc_clip_max` (positive upper limit in mm2/s; opt-in legacy derivative). Resampling uses an identity transform in physical coordinates, not array resizing, header replacement or registration. The runner checks coverage of the ROI after resampling and fails when any tumor voxel lacks valid positive low/high signal. Non-overlap or motion requires investigation, not another interpolation attempt.

Outputs for each successful case:

- `ADC_raw.nii.gz` and `ADC_valid.nii.gz` for ADC mode; optional `ADC_clipped.nii.gz` when requested.
- `ROI_aligned.nii.gz`, preserving nonnegative integer labels; original image grid is preserved for overlay mode.
- `overlay.png`: three native-grid slices near ROI extent quartiles, each showing grayscale image beside red ROI overlay. Slice indices and display window are in QC. This is not a standardized anatomical orientation or complete-volume review.
- `qc.json`: input geometry, hashes, parameters, invalid/negative ADC fractions, coverage and output review status.
- `run_summary.json`: software versions, manifest/script hashes and all case outcomes.

Review each volume in a medical image viewer in addition to the limited PNG. A suggested `review.json` entry is `{ "id": "case001", "decision": "approved|rejected|pending", "reviewer": "pseudonym", "reason": "..." }`. Only approved cases advance to feature extraction. NaN ADC voxels outside the ROI require an explicit, documented downstream extraction policy; never silently treat them as measured ADC zero.

