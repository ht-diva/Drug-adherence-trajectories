# Drug-adherence-trajectories

Pipelines (Finregistry / Finland and Regione Lombardia / Italy) to construct long-term medication purchase histories, compute adherence time-series (e.g., MPR), smooth trajectories, and run functional trajectory modelling with downstream clustering and association analyses.

> **Important:** This repository contains code only. The underlying registry data are not included. Many scripts use environment-specific/hard-coded paths that you will need to adapt.

---

## Data inputs (high level)

You will need individual-level longitudinal purchase/dispensing data with variables of this type (names differ by source):

- person identifier (e.g., `ID`)
- dispensing/purchase date (e.g., `DATE`)
- medication identifier (e.g., ATC code)
- quantity / pack information to derive pills and days supplied
- optional daily dose assumptions/derivations (some scripts filter to 1 tablet/day)

The Lombardia pipeline includes logic to merge nearby purchases (e.g., within 7 days), compute pills per purchase, and derive both **point** and **cumulative** MPR (including capped versions at 1).

---

## Workflow (analysis scripts)

Scripts are numbered to reflect the intended order.

1. **Process purchase histories**
   - `Finregistry/Analysis/1_process_trajectories.R`
   - `RegioneLombardia/Analysis/1_process_purchases.R`

2. **Select cohort & summarise**
   - `Finregistry/Analysis/2_select&summarise.R`
   - `RegioneLombardia/Analysis/2_select_summarise.R` (and `_4covid.R` variant)

3. **Smooth adherence trajectories**
   - `Finregistry/Analysis/3_smooth_adherence.R`
   - `RegioneLombardia/Analysis/3_smooth_traj.R`

4. **Functional data analysis**
   - `*/Analysis/4_fda.R` with helper functions in `*/Analysis/fda_funs.R`

5. **Create covariates / cluster trajectories**
   - `5.1_create_endpoints.R`, `5.2_create_socioeco.R`, `5.3_create_medication.R`
   - Lombardia additionally: `5.4_create_visits.R`
   - Clustering scripts:
     - `Finregistry/Analysis/5_cluster_curves.R`
     - `RegioneLombardia/Analysis/5_cluster.R`

6. **Assemble covariate table**
   - `*/Analysis/6_covariates.R`

7. **Fit models**
   - `*/Analysis/7_models.R`
   - Lombardia also: `models_PC.R` and `summary_stats.R`

---


