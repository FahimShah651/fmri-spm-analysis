# fMRI Brain Data Analysis using SPM and MATLAB

**Author:** Fahim Ur Rehman Shah (EE2629)

**Course:** CSE532 — Signal and Image Processing (MS Level)

**Supervisor:** Dr. Adnan Shah, FCSE, GIKI

**Institution:** GIK Institute of Engineering Sciences & Technology

---

## Overview

End-to-end fMRI analysis pipeline including:
- Preprocessing (Slice Timing, Realign, Coreg, Segment, Normalise, Smooth)
- General Linear Model (GLM) estimation and contrast mapping
- Functional connectivity analysis (seed-based)
- Cross-dataset comparative analysis

## Datasets

1. **ds000114** — Test-retest fMRI (Motor task: finger/foot/lips)
   - Subject: sub-01, TR: 2.5s, 40 slices
   - URL: https://openneuro.org/datasets/ds000114

2. **ds000105** — Visual object recognition
   - Subject: sub-1, TR: 2.5s, 40 slices
   - URL: https://openneuro.org/datasets/ds000105

## Pipeline Steps

| # | Script | Description |
|---|--------|-------------|
| 01 | `step01_setup_and_download.m` | Environment check, directory setup, data verification |
| 02 | `step02_dataset1_preprocessing.m` | Full preprocessing — Dataset 1 |
| 03 | `step03_dataset2_preprocessing.m` | Full preprocessing — Dataset 2 |
| 04 | `step04_dataset1_glm.m` | GLM and contrasts — Dataset 1 |
| 05 | `step05_dataset2_glm.m` | GLM and contrasts — Dataset 2 |
| 06 | `step06_functional_connectivity.m` | Seed-based connectivity — Both |
| 07 | `step07_comparison_analysis.m` | Cross-dataset comparison |

## Requirements

- MATLAB R2022b or later
- SPM25 (Statistical Parametric Mapping)
- Git

## How to Run

```matlab
>> run('step01_setup_and_download.m')
>> run('step02_dataset1_preprocessing.m')
```

## Output Structure

```
figures/          ← PNG + EPS figures (IEEE naming)
results/          ← MATLAB .mat results and CSV tables
latex/            ← IEEE-format LaTeX report
project_log.txt   ← Running log of all steps
```

---
*Generated on: 22-May-2026 14:01:12*
