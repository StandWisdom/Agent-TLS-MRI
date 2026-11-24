# ROI input and output contract

The manifest contains `cases`, each with a unique path-safe `id`, `image`, `mask`, `sequence` and `acquisition_note`. Paths are relative to the manifest or absolute local paths. Image and mask must be scalar 3D NIfTI/MHA volumes.

Use an existing MRI or selected ADC derivative as the reference image. Calculated ADC requires `adc_source: calculated` and `valid_mask` pointing to its binary `ADC_valid.nii.gz` on the same grid. The valid-signal map preserves invalid voxel locations independently of NIfTI NaN handling. Use a tumor mask drawn on the corresponding sequence or an independently verified registered image. Matching voxel dimensions alone are insufficient: origin, spacing, direction and anatomy must match.

Optional fields: `allow_resample` (default false) and `resample_reason` (required when true). Resampling maps labels by nearest neighbor using physical coordinates; it is not registration. Foreground means all positive mask labels. The subsequent radiomics YAML determines which label is extracted (default 1).

Outputs per case:
- `ROI_aligned.nii.gz`: validated labels on the reference image grid.
- `overlay.png`: native-grid slices at ROI extent quartiles, grayscale beside red foreground overlay.
- `qc.json`: input geometry/hashes, resampling state, ROI voxel count/volume, overlay slice indices/window and review status.
- Run-level `run_summary.json`: cases, parameters and software provenance.

Review the full volume in addition to the preview. A local review entry may contain id, reviewer, approved/rejected decision and reason. Do not advance failed or unreviewed cases. ADC NaN values inside the tumor are rejected here; NaNs outside the tumor still need an explicit downstream policy before extraction.
