---
name: roi-preparation
description: Validate existing MRI tumor masks, align their physical grids and generate image/ROI overlays for radiomics quality control.
---

# ROI preparation and quality control

Read [protocol](references/protocol.md) and use [example.json](assets/example.json). Install beside `mri-preparation`, which provides shared image geometry and batch-QC utilities.

1. Obtain an existing tumor mask for the confirmed patient, acquisition and MRI sequence. This runner processes supplied segmentations; it does not generate a tumor segmentation.
2. Create a manifest with the MRI/ADC `image`, `mask` and acquisition note. For computed ADC, set `adc_source: calculated` and supply its `ADC_valid.nii.gz` as `valid_mask` along with the chosen ADC image.
3. Run `python skills/roi-preparation/scripts/prepare_roi.py data/roi_manifest.json --output runs/roi`.
4. Inspect full volumes and the overlay. Check tumor coverage, label values, anatomy and sequence identity. Record approval in a local review record.
5. Pass the original image and `ROI_aligned.nii.gz` to radiomics extraction only after approval.

Masks must contain finite, nonnegative integer labels and a nonempty foreground. Grid mismatch fails unless resampling is explicitly permitted with a reason. Resampling uses nearest-neighbor interpolation on the image's physical grid; labels are preserved. More than 10% ROI physical-volume loss, a nonfinite image value, or invalid ADC signal inside the aligned ROI fails the case. ADC validity is checked using the saved binary mask because NIfTI readers may replace NaN values.

Outputs include aligned masks, three-slice grayscale/ROI overlays, geometry QC and source/output hashes. Overlays are previews, not inputs for radiomics. Completed cases remain `REVIEW_REQUIRED` until anatomical review. Preserve source masks and keep patient-level artifacts local.
