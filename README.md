# Radiomics Skill Harness

Six reusable agent skills connect multiparametric MRI, pathology-derived TLS measurements and response prediction. Each skill combines task instructions with executable analysis code.

## Workflow

| Skill | Operation | Output |
|---|---|---|
| [MRI preparation](skills/mri-preparation/SKILL.md) | Two-b-value ADC, ROI grid alignment, image/ROI overlays | Images, overlays, QC |
| [Radiomics extraction](skills/radiomics-extraction/SKILL.md) | ADC / T1C_ph1 / T1C_ph3 / T2 with a shared YAML | Four feature workbooks |
| [TLS correlation selection](skills/tls-correlation-selection/SKILL.md) | Pearson p < 0.05 and absolute r > 0.1 | Candidate features and all tests |
| [TLS feature reduction](skills/tls-feature-reduction/SKILL.md) | BH-adjusted evidence and bootstrap stability | Reduced feature subsets |
| [Internal model development](skills/internal-model-development/SKILL.md) | Stratified 70% training / 30% validation, TabPFN | Frozen model and held-out predictions |
| [External model validation](skills/external-model-validation/SKILL.md) | Apply the same frozen model to each independent center | Center-specific predictions and metrics |

Keep the `skills/` directories together: correlation/reduction share a runner, and development/external validation share model utilities. An agent reads the relevant `SKILL.md`, prepares its local inputs, invokes the script and checks its outputs. MRI alignment review is a human handoff.

## Installation

Python 3.11 is the tested runtime. From the repository root:

```sh
python -m pip install numpy==1.26.4 "setuptools<81" wheel
python -m pip install --no-build-isolation -r requirements-tested.txt
python -m pip install -r requirements-model-tests.txt
python -m unittest discover -s tests -v
```

PyRadiomics may require a C compiler (Windows: C++ build tools and SDK). The requirement files record tested direct dependencies. Use a separate server environment for the selected TabPFN/PyTorch versions and local checkpoint; the model-test requirements install the synthetic-test backend. GPU TabPFN execution must be tested in that server environment.

## Running the workflow

Prepare local manifests from the skill assets. Read the [TLS schema](skills/tls-correlation-selection/references/tls-schema.md) and [model input contract](skills/internal-model-development/references/input-contract.md). Output directories must be new.

```sh
python skills/mri-preparation/scripts/prepare_mri.py data/mri_manifest.json --output runs/prepared
python skills/radiomics-extraction/scripts/extract.py data/extraction_manifest.json --yaml skills/radiomics-extraction/assets/legacy_mr.yaml --output runs/extracted
python skills/tls-feature-reduction/scripts/select_tls_features.py --stage legacy --features runs/extracted --tls data/private_tls.xlsx --output runs/legacy
python skills/tls-feature-reduction/scripts/select_tls_features.py --stage reduced --features runs/extracted --tls data/private_tls.xlsx --output runs/reduced
```

Both selection stages read the original feature tables. BH uses all finite tests within each sequence-marker family; the reduced set is restricted to legacy candidates. The default reduction uses FDR evidence **or** bootstrap stability, plus effect-size and direction gates. `--combine and` requires both. See the reduction skill for exact thresholds.

Each selection workbook contains `Features`, `AllTests` and `Counts`. `Features` is the union across TLS markers; a feature associated with multiple markers contributes to multiple association counts. Use the relevant `AllTests` rows to build an Agg-, FL1- or FL2-specific signature.

Create the frozen signature and join other centers' feature tables using the [model input contract](skills/internal-model-development/references/input-contract.md). Replace `YOUR_INSTALLED_VERSION` with the exact TabPFN package version:

```sh
python skills/internal-model-development/scripts/develop_model.py --data data/development.xlsx --signature data/signature.json --discovery-ids data/discovery_ids.xlsx --label-mode trg --positive-label 0 --checkpoint checkpoints/model.ckpt --tabpfn-version YOUR_INSTALLED_VERSION --device cuda --output runs/development
python skills/external-model-validation/scripts/validate_external.py predict --bundle runs/development --features data/center_a_features.xlsx --device cuda --trust-local-model --output runs/center_a_predictions
python skills/external-model-validation/scripts/validate_external.py score --predictions runs/center_a_predictions --labels data/center_a_outcomes.xlsx --output runs/center_a_evaluation
```

Repeat `predict` and `score` for each external center, keeping the same development bundle. The model is fitted on the 70% training partition and is not refitted after internal validation. Optional scaling/PCA is also fitted on that partition; both are off by default. Default positive label is 0 (pCR) and probability threshold is 0.5.

Model outputs include split workbooks, held-out predictions, AUC, AUPRC, ACC, sensitivity, specificity, PPV, NPV, F1, binary Kappa, Brier, calibration bins and DCA tables. See the [validation protocol](skills/internal-model-development/references/validation-protocol.md) for data isolation and model-loading requirements.

## Public package

Only code, skill instructions, parameter/schema examples and synthetic tests are included. Keep patient tables, images, results, fitted models and private contracts in local `data/` or `runs/`. Fitted TabPFN archives can contain patient training context. Load only trusted model artifacts.

Build the upload package using the explicit source-file list:

```sh
python scripts/build_public_release.py --output dist/radiomics-skill-harness-public.zip
```

License selection is pending team confirmation; no license is included.
