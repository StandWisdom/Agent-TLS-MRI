---
name: mri-preparation
description: Prepare MRI inputs and compute two-b-value ADC maps with geometry, signal checks and provenance before ROI processing.
---

# MRI preprocessing

1. Read [protocol](references/protocol.md). Match patient, sequence and acquisition time. Establish b-values from metadata or a documented acquisition protocol; confirm signal scaling and anatomical alignment.
2. Create a manifest using [example.json](assets/example.json). Use `mode: adc` for paired low/high-b inputs or `mode: image` to inspect an existing ADC/T1/T2 volume. Inputs are scalar 3D NIfTI/MHA; resolve DICOM series and 4D time points upstream.
3. Run `python skills/mri-preparation/scripts/prepare_mri.py data/mri_manifest.json --output runs/mri` from the repository root.
4. Review image identity, geometry and ADC signal QC. Forward the original image or chosen ADC derivative to `roi-preparation`.

ADC uses ln(S_low/S_high)/(b_high-b_low) and is saved in mm2/s. Invalid signal voxels remain NaN and are recorded in a separate valid-signal mask. Raw negative/high ADC values are retained; optional clipping produces a separate derivative.

Grid mismatch fails by default. Low-b resampling requires `allow_resample: true` and `resample_reason` documenting a shared physical coordinate system. Linear interpolation uses the high-b reference. Resampling does not correct motion.

This skill takes no tumor ROI. ROI coverage, alignment and overlays are handled by `roi-preparation`. Preserve source images and local QC. All technically completed cases remain `REVIEW_REQUIRED` until acquisition review.
