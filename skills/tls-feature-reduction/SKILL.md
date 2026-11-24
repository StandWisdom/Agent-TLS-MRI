---
name: tls-feature-reduction
description: Reduce legacy TLS-associated MRI candidates using BH-adjusted association evidence and bootstrap selection/direction stability, retaining original feature meanings and a complete audit trail.
---

# Stability-aware feature reduction

Run `scripts/select_tls_features.py --stage reduced` after or together with legacy screening. This is feature subset selection (dimensionality reduction by dropping columns), not PCA or learned feature transformation.

The runner calculates legacy candidates at p<0.05 and |r|>0.1. BH-adjusted p-values are calculated from ALL finite feature tests for each sequence-marker family, never only legacy survivors. Bootstrap draws ceil(0.8*n) rows WITH replacement for 1000 iterations by default. An iteration selects a feature when it has the same direction as the full estimate, p<0.05 and |r|>=0.20. Degenerate marker resamples are skipped; denominators are recorded.

Default reduction: original p<0.05, |r|>=0.20, direction consistency>=0.80 AND (BH q<0.10 OR bootstrap frequency>=0.65). `--combine and` requires both BH and bootstrap gates. Set parameters explicitly through the CLI. Two-sided Pearson p-values use the beta-distribution implementation in both full-sample and bootstrap tests.

BH is applied separately to each sequence-marker family. The default OR rule can retain bootstrap-stable features with q>=0.10; report adjusted p-values and selection frequencies as separate quantities. Bootstrap measures selection variability within the discovery cohort.

Keep thresholds fixed before outcome evaluation. For TLS prediction perform the entire selection inside training folds. For downstream response prediction, freeze the signature from a disjoint TLS cohort and intersect feature names, not patients, with other centers. No automatic threshold relaxation to reach a desired feature count. Output an empty signature if no features qualify.
