---
name: mri-preparation
description: Prepare research MRI and existing tumor ROIs, compute two-b-value ADC maps, inspect physical geometry, and generate image/ROI overlays with provenance and review flags. Use before radiomics feature extraction.
---

# MRI preparation

Use this skill for ADC calculation, existing ROI alignment, and visual quality control. This module does not perform tumor segmentation, deformable registration, feature extraction, or clinical interpretation.

## Execution

1. Inventory the requested inputs locally. Match pseudonymous patient IDs, acquisition time points and sequences. Read [protocol.md](references/protocol.md) before preparing a manifest.
2. Establish b-values from acquisition metadata or a documented acquisition protocol. Establish enhancement timing independently: a filename such as DCE or T1+C does not prove phase 3. Stop affected cases when these are unresolved.
3. Create a JSON manifest following [example.json](assets/example.json). Only include de-identified NIfTI/MHA scalar 3D inputs in this runner. Resolve DICOM series and 4D volumes explicitly upstream; never choose the largest series or first time point silently.
4. Run `python scripts/prepare_mri.py manifest.json --output new_run_directory`. Use a fresh output directory. Each case emits QC and provenance; the command returns a nonzero exit status if any case failed.
5. Inspect the raw image, mask and overlay, not only the QC numbers. All technically completed cases remain `REVIEW_REQUIRED` until a reviewer confirms anatomical alignment and acquisition identity. Record approval in `review.json` separately from machine QC.
6. Handoff only reviewed outputs. Preserve the manifest, checksums, software versions and failed-case reasons. Never overwrite source images or replace failed cases with zeros.

## Agent decisions

- Matching grids do not prove matching anatomy. Resampling is not registration.
- A geometry mismatch defaults to failure. `allow_resample` is permitted only after documented confirmation of a shared physical coordinate system. Use linear interpolation for low-b intensity and nearest-neighbor for ROI labels; high-b is the ADC reference.
- ADC inputs must be comparable diffusion-weighted acquisitions with known b-values in s/mm2. Record intensity scaling compatibility, motion assessment and spatial coverage; arbitrary MRI normalization must not precede the signal ratio.
- Save ADC in mm2/s and a valid-signal mask. Undefined voxels remain NaN. Retain negative/high values in the raw ADC and flag them. Clipping requires an explicit manifest parameter and creates a separate derivative.
- The overlay is a visual derivative, never a radiomics input. Display percentile windowing affects the PNG only.
- Use patient-level files locally; do not publish their images, identifiers, manifests, paths or hashes without the project's data-sharing authorization.

## Verification and reporting

Run `python -m unittest discover -s tests -v` from the repository root. Report processing counts, errors, review state and output paths.
