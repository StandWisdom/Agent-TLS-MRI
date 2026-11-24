---
name: tls-correlation-selection
description: Match private patient-level TLS measurements with MRI radiomics and perform fixed-threshold legacy Pearson screening, retaining all tests for subsequent correction.
---

# Legacy TLS association screening

Read [TLS schema](references/tls-schema.md). Accept a local private XLSX with one ID and 12 TLS columns. Never include that workbook, real rows, counts or identifiers in the public repository. Use the shared runner in `../tls-feature-reduction/scripts/select_tls_features.py` with `--stage legacy`. Install both sibling skills together; they share the numerical implementation.

Evaluate each sequence and each TLS marker separately. Match IDs exactly as strings, reject duplicates and invalid cells, report unmatched IDs locally. Do not collapse `001` to `1` without a verified mapping. Use Pearson two-sided p-values; binary markers correspond to point-biserial correlation. Ordinal labels require documented coding; Pearson assumes the specified numeric spacing.

Legacy candidates satisfy STRICT `p < 0.05` AND `abs(r) > 0.1`. No FDR or bootstrap gate is applied in this stage. Export per-sequence candidate feature values and every tested feature-marker pair, including unselected tests, for the reduction stage. The selected feature table is the union across markers, deduplicated by feature name within sequence; association counts and unique feature counts differ.

For supervised TLS prediction, run screening inside training partitions. For downstream response prediction, freeze the signature in a separate, nonoverlapping TLS discovery cohort.
