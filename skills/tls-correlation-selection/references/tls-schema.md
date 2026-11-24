# Private TLS table schema (no patient data)

One row represents a patient. MRI and pathology share a pseudonymous `ID`; joining uses this ID rather than row order. The public code accepts a user-provided `.xlsx`; no actual pathology workbook or real example rows are distributed.

| Column | Meaning | Contract |
|---|---|---|
| ID | Pseudonymous patient key | Unique, nonempty string |
| tls_area_total | Total annotated TLS area | Nonnegative; unit must be documented |
| tls_area_max | Largest annotated TLS area | Same unit |
| tls_area_mean | Mean TLS region area | Same unit |
| tls_area_std | TLS region area standard deviation | Document population vs sample SD |
| Agg_num | Count of Agg regions | Nonnegative integer |
| FL1_num | Count of FL1 regions | Nonnegative integer |
| FL2_num | Count of FL2 regions | Nonnegative integer |
| Total_num | Total TLS region count | Confirm how unclassified regions are counted |
| Agg_flag | Agg label | Confirm 0/1 mapping with data dictionary |
| FL1_flag | FL1 label | Confirm 0/1 mapping with data dictionary |
| FL2_flag | FL2 label | Confirm 0/1 mapping with data dictionary |
| TLS_012_flag | Combined TLS category | Confirm exact encoding; do not infer presence of all three classes |

Agg, FL1 and FL2 describe pathology-defined maturation categories. Record original pathology definitions and how multiple slides/regions are aggregated. Do not infer maturation classes from unlabeled XML contours, or infer mm2 from pixel areas without slide calibration. Missing values are not zero. Empty regions, absent TLS and unavailable pathology have distinct meanings.
