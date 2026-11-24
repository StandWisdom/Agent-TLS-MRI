---
name: radiomics-extraction
description: Extract ADC, T1C_ph1, T1C_ph3 and T2 radiomics from reviewed 3D MRI and tumor masks using an explicit PyRadiomics YAML protocol, with per-sequence Excel outputs and provenance.
---

# Four-sequence radiomics extraction

Read the repository README and the extraction manifest example. Use reviewed images from MRI preparation and matching aligned masks from ROI preparation. Confirm enhancement phases; T1+C and DCE names alone do not establish phase identity. Preserve ADC units and note whether the input is vendor ADC, calculated raw ADC or a clipped derivative.

Run `python skills/radiomics-extraction/scripts/extract.py manifest.json --yaml protocol.yaml --output runs/extraction001` from the repository root.

`assets/legacy_mr.yaml` specifies normalization on all sequences including ADC, scale 100, bin width 5, voxel shift 300, forced 2D along dimension 0, Original/LoG sigma 3/Wavelet and the listed GLCM features. Confirm image axis orientation before using forced 2D. Record ADC normalization because it changes the quantitative intensity scale.

The old YAML has geometryTolerance 1e5. The runner first checks size, spacing, origin and direction strictly, then overrides this tolerance to 1e-5 and disables automatic mask correction. No geometry repair occurs during extraction. Send image/mask mismatches back to ROI preparation for review. Both requested YAML hash and effective settings are recorded.

Use an explicitly documented finite-image derivative when ADC has NaN outside the tumor; global normalization/filtering can be affected by outside-ROI values. The runner rejects all nonfinite intensities rather than filling silently.

Outputs: one `*_features.xlsx` per sequence (ID and numerical features only), an extraction QC JSON and frozen YAML. Missing sequences and failed cases are recorded; do not replace missing feature rows with fabricated values.
