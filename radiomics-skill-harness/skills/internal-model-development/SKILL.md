---
name: internal-model-development
description: Develop a local TabPFN response model using a fixed patient-level 70/30 training-validation split, training-only preprocessing and a frozen model bundle for independent centers.
---

# Internal development: 70/30

Read [validation protocol](references/validation-protocol.md) and [input contract](references/input-contract.md). Use one stratified 70% training / 30% internal validation split, with seed and threshold fixed before evaluation.

Use the frozen independent-TLS feature signature. Require globally consistent patient keys and check discovery/development disjointness. This runner does not support evaluating TLS classification in the same patients used for TLS screening; such tasks require selection inside training partitions. For response development with an overlapping TLS cohort, stop and redesign selection rather than bypassing the guard.

Run `scripts/develop_model.py` using the repository-root example in README. It fits optional scaling/PCA only on the 70% training set; by default both are off. Finite input values are required. TabPFN receives only training cases in its `fit` call.

The default probability threshold is fixed at 0.5. Explicitly choose raw four-grade TRG or already-binary labels. With `--label-mode trg`, map TRG=0 to label 0 and TRG=1/2/3 to label 1. Default positive label is 0, meaning pCR; never assume probability column 1 represents pCR.

Persist the SAME 70%-trained estimator and preprocessor used for internal validation. Do not refit on the combined 100% development data. Export private split tables and held-out predictions, metrics including Kappa/Brier, and descriptive calibration/DCA tables. These predictions are held-out validation predictions, not OOF. No independent test set is created inside the development center.

Require an explicit local checkpoint and exact installed TabPFN version. Keep fitted state and patient-level outputs private as specified in the validation protocol. The `logistic-smoke` backend runs synthetic software tests.
