---
name: external-model-validation
description: Apply a frozen development-center model directly to independent centers, with feature-only prediction separated from label-based scoring and no external fitting or adaptation.
---

# Direct external validation

Install this skill beside `internal-model-development`; the executable shares its data and metric implementation. Read its `references/validation-protocol.md` and `references/input-contract.md` before execution.

Use a completed frozen 70%-training bundle. For EACH external center, run `scripts/validate_external.py predict` against a label-free feature workbook. The exact feature names/order, model, preprocessor, positive class and threshold are inherited from development. The runner rejects extra label columns, missing features and patient overlap with the entire development or TLS discovery cohort. Do not rename an unverified enhancement phase to satisfy the schema, and do not drop missing sequences post hoc. Develop a prespecified common-sequence model if necessary.

After saving predictions, run `score` with that center's outcome workbook. The scorer verifies the prediction checksum and uses the frozen threshold. Keep each center's results separate.

Apply the development preprocessor directly; do not fit transformations, resample cases or recalibrate using the external cohort.

Load trusted local bundles only. Keep the original checkpoint at its recorded location and restore the recorded environment. Keep fitted archives private because they can contain patient training context.
