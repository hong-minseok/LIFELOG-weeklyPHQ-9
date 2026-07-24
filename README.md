# Analysis Code for *Structure, Stability, and Clinical Validity of the Weekly PHQ-9 in a Clinical and Community Mood Monitoring Cohort*

## Manuscript

**First author:** Minseok Hong

- Department of Neuropsychiatry, School of Medicine, Eulji University, Daejeon, Republic of Korea
- Department of Psychiatry, Uijeongbu Eulji Medical Center, Eulji University, Uijeongbu, Republic of Korea

**Corresponding author:** Sang Jin Rhee

- Department of Neuropsychiatry, Seoul National University Hospital, Seoul, Republic of Korea
- Department of Psychiatry, Seoul National University College of Medicine, Seoul, Republic of Korea
- Emails: hellojr1123@hanmail.net or najr0722@snu.ac.kr

**Status:** Under submission

## Analysis code

The scripts in [`analysis-code/`](analysis-code/) are intended to be read in numeric order.

| Script | Statistical content |
|---|---|
| `0_analysis_definitions.R` | Analysis-sample rules, PHQ model definitions, and lavaan syntax helpers |
| `1_factor_structure.R` | Even-week EFA, odd-week ordinal CFA, and ordinal omega |
| `2_measurement_invariance.R` | Longitudinal and clinic-community measurement invariance, including week-specific estimability checks |
| `3_multilevel_structure.R` | Multilevel CFA, item variance decomposition, and loading isomorphism |
| `4_cluster_bootstrap.R` | Participant-level cluster bootstrap for correlations, AUC, and cutoff estimates |
| `5_external_validity.R` | MADRS correlations, concurrent MDE discrimination, external-criterion correlations, and Youden cutoffs |

## Variable definition

| Variable | Definition |
|---|---|
| `studyID` | Participant identifier (character) |
| `studyWeek` | Study week (numeric) |
| `phq9_1`-`phq9_9` | PHQ-9 items (0-3) |
| `phq9_sum` | Complete-item PHQ-9 total |
| `cohort_source` | `Clinic` or `Community` |
| `madrs01`-`madrs10` | MADRS items |
| `madrs_sum` | Complete-item MADRS total |
| `cgis` | CGI-S score |
| `hama01`-`hama14` | HAMA items |
| `ymrs01`-`ymrs10` | YMRS items |
| `gad2_1`-`gad2_2` | GAD-2 items |
| `gad2_sum` | Complete-item GAD-2 total |
| `pss01`-`pss10` | PSS items |
| `ruls1`-`ruls6` | ULS-6 items |
| `mhc01`-`mhc14` | MHC-SF items |
| `mde_base` | Baseline MDE status |
| `mde_followup` | Follow-up MDE status |
| `MDE` | Concurrent MDE status used in validity analyses |

## Analysis environment

The analyses used R 4.6.0 with the following direct dependencies:

- lavaan 0.6-21
- semTools 0.5-8
- psych 2.6.5
- pROC 1.19.0.1
- dplyr 1.2.1 and tidyr 1.3.2 for data handling

Scientific rationale, complete methods, results, tables, and figures are provided in the manuscript and supplement.

## License

Code is released under the [MIT License](LICENSE).
