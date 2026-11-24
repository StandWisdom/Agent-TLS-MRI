# Validation protocol

## Data flow

1. Freeze the feature signature in an independent TLS discovery cohort. Reconcile patient identities across discovery, development and external centers; the runner rejects overlapping IDs.
2. Split the development center by patient with stratification: 70% training, 30% validation. Fix the seed and parameters before evaluation. Sample counts follow scikit-learn rounding.
3. Fit the preprocessing and classifier on the 70% training partition. Evaluate on the 30% once, then save that same estimator.
4. Apply the saved estimator and preprocessor directly to each external center. Save label-free predictions before joining outcomes for scoring.

There is no full-development refit, external resampling, center-specific adjustment or recalibration. Any future tuning belongs inside the 70% training partition; the 30% and external centers remain evaluation sets. A signature derived from overlapping patients requires training-partition feature selection instead of this frozen-signature workflow.

## Outputs

Internal predictions are a held-out validation set, not OOF. Kappa compares binary pCR/non-pCR labels. AUC/AUPRC are missing for single-class cohorts. Calibration uses fixed probability bins; DCA uses thresholds from 0.05 to 0.95. The exported tables contain point estimates and no confidence intervals. External centers are reported separately.

## Private artifacts

Keep patient tables, split assignments, predictions, fitted archives, logs and private contracts local. Fitted TabPFN state can retain training context. Use the local TabPFN package with a licensed checkpoint and restricted network access for patient-data runs. Only load trusted bundles through `--trust-local-model`; serialization can execute code. Retain the recorded checkpoint location and environment for loading. SHA256 records detect artifact changes.

## References

- [scikit-learn preprocessing and leakage](https://scikit-learn.org/stable/common_pitfalls.html)
- [scikit-learn train_test_split](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html)
- [TabPFN fitted-model persistence](https://github.com/PriorLabs/TabPFN/blob/main/examples/save_and_load_model.py)
