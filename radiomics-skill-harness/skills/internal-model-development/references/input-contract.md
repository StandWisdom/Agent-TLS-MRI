# Local-only input contract

All patient tables remain outside the source repository (recommended: ignored `data/` and `runs/`). Each row is one patient, one pre-treatment observation. Repeat scans/slides must be consolidated before modelling. ID must be a globally reconciled pseudonymous patient key: center-specific numbering or prefixes alone cannot prove that two records are different people. Keep the reidentification mapping private.

1. `signature.json`: `origin` is `independent_tls_discovery`; `features` is the ordered list of frozen sequence-prefixed radiomics column names. See the synthetic schema example. Build the list from upstream selected feature names, not outcome association in development/external centers. For Agg/FL1/FL2-specific models, select the corresponding upstream `AllTests` rows; do not accidentally use the union of all markers.
2. Discovery ID workbook: column `ID`, covering ALL patients whose values influenced feature screening/reduction, not only the final survivors. No TLS measurements are needed at model development time.
3. Development workbook: `ID`, outcome column (default `TRG`), and exactly the signature feature columns. Join the sequence tables by ID and rename, e.g. `original_firstorder_Mean` in ADC to `ADC::original_firstorder_Mean`. The signature fixes which sequences are required; it is not restricted to four-sequence models.
4. External prediction workbook: `ID` plus exactly the SAME feature columns. It must contain NO outcomes or clinical covariates. Missing/extra columns or nonfinite features cause failure rather than a silent change of model.
5. External outcome workbook: `ID` plus the configured outcome column only, with the same patient set as the frozen predictions. TRG mode accepts integers 0/1/2/3; binary mode accepts 0/1 only. Missing labels are never interpreted as non-pCR.

Raw TRG and already-binarized outcomes are not interchangeable modes. Default output probability means P(pCR) because positive label=0; this convention is recorded in the private contract. Change the positive label only as a prespecified endpoint decision, consistently for all centers.

Keep train/validation assignment fixed across old/agent comparisons by using identical eligible patient IDs, labels and seed. Compare split workbooks to confirm equality; using the same seed with different eligibility sets is insufficient. Freeze eligibility before outcome analysis. Document exclusions rather than silently intersecting a different patient subset for each method.
