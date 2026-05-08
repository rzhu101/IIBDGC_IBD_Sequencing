# IIBDGC_IBD_Sequencing
This repository contains analysis code supporting the manuscript:

> **Exome sequencing directly implicates 68 genes in inflammatory bowel disease**  


---

## Overview

This study performed large-scale whole-exome sequencing meta-analysis across IBD cases and controls (82,142 IBD cases total) to identify independent rare coding associations with Crohn's disease (CD), ulcerative colitis (UC), and IBD. The repository covers scripts for quality control (Broad datasets), variant annotation, association analysis, meta-analysis, LD-based conditional analysis.

---


## Dependencies

| Tool / Package | Version | Usage |
|---|---|---|
| Hail | 0.2.133 | QC, variant annotation |
| VEP | 112 | Variant consequence annotation |
| LOFTEE | — | High-confidence PTV flagging |
| AlphaMissense | — | Missense pathogenicity scoring |
| REGENIE | 3.2.2 / 4.1 | Single-variant and burden association testing |
| METAL | — | Inverse-variance weighted meta-analysis |
| PLINK | 1.9 | Pairwise LD estimation |
| CELLEX | 1.2.1 | Cell-type expression specificity |
| Python | ≥ 3.8 | Analysis scripting |
| R | ≥ 4.0 | Statistical modeling, visualization |

---

## Data Availability
Summary statistics from the meta-analysis are deposited at [].

---

## Citation