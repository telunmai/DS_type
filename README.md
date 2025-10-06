# Detecting Phenotype-Filtered Differential Co-expression (DS) shifts in small-n RNA-seq

This repository contains two R scripts to compute per-condition gene–gene correlations and classify **differential scenarios (DS1/DS2/DS3)** under a phenotype-anchored design (e.g., DMSO → GL24) with **small-n replicates per group**.

- **Script 1** builds per-group correlation statistics and permutation empirical P-values.
- **Script 2** calls **S/NS/NA** per condition, derives per-line transitions, and classifies **DS1/DS2/DS3**.

> Designed for small-n perturbation studies where edge-level changes are enforced to track a phenotype (e.g., “effective” vs “non-effective” cell lines).

---

## Quick start

```bash
# 1) Build per-group stats (V_, Cor_, P_, empP_)
Rscript Step1.Correlation_for_DS.R

# 2) Classify DS1/DS2/DS3 from the stats table
Rscript Step2.Identify_DS_types.R
```

Outputs:
- `DS_pairwise_cor_stats.tsv` (from Script 1)
- `DS_types.tsv` (from Script 2)

---

## Inputs & file formats

### Expression matrix
- **File:** `TPM_matrix_sample.txt`
- **Shape:** genes × samples, with one **gene ID column**.

### Group design
- Example cell_lines/tags used downstream:  
  `MDAMB231_DMSO / MDAMB231_GL24`, `MDAMB157_DMSO / MDAMB157_GL24`, `Hs578T_DMSO / Hs578T_GL24`

---

## What Script 1 computes

For every gene pair and group:

- `V_<group>`: 1 if a valid Pearson `cor.test()` ran and returned finite `r` and `P`; 0 otherwise.
- `Cor_<group>`: Pearson correlation coefficient (*r*).
- `P_<group>`: uncorrected P from `cor.test()`.
- `empP_<group>`: permutation empirical P = mean(`P_perm < P_obs`) over `B_perm` shuffles (default **200**).

**Edge cases handled**
- All-zero expression in either gene ⇒ NA (with `V=0`).
- Requires ≥3 finite points and ≥3 unique (x,y) pairs.
- Any permutation failure for a pair/group ⇒ `empP=NA`.

**Reproducibility**
- `set.seed(b)` is called **inside** the permutation loop; this yields deterministic, shared permutation streams across pairs for a given `b`.

---

## What Script 2 produces

Reads `DS_pairwise_cor_stats.tsv` and:

1. **Labels state per condition** (per gene pair):
   - **S**: `P < 0.005` **and** `empP < 0.005` (both finite)
   - **NS**: `P > 0.05` (finite; `empP` not required)
   - **NA**: everything else (including `V==0` or ambiguous)
2. **Defines per-line transitions** (before → after):
   - `(NS/NA) → S` = **gain**
   - `S → (NS/NA)` = **loss**
   - `NS/NA ↔ NS/NA` or `S → S` = **unchanged**
3. **Calls DS categories** using effective vs control lines:
   - **DS1 (gain):** all **effective** lines show **gain**, and **control** lines are **unchanged** (or, if relaxed, do not show gain).
   - **DS2 (loss):** all **effective** lines show **loss**, and **control** lines are **unchanged** (or, if relaxed, do not show loss).
   - **DS3 (change):** all **effective** lines are **S→S** with  
     \[
       |Δr| = \big||r_\text{after}| - |r_\text{before}|\big| \ge \texttt{delta_r_thresh}
     \]
     Controls must be **NS/NA↔NS/NA** or **S→S** with **unchanged** |Δr|.
   - Default **`delta_r_thresh = 0.10`**.
   - Control behavior toggled by **`require_NE_unchanged`** (default `TRUE`).

**Outputs**
- `DS_types.tsv`: rows for DS1/DS2/DS3 hits with:
  - `GeneA`, `GeneB`
  - per-line **before/after states** and **transition** (`T_<Line>`)
  - For DS3 display, `T_<Line>` is replaced by `"changed"/"unchanged"` when S→S is applicable.

---

## Parameters

- **Permutation size:** `B_perm` (default 200). Increase for smoother `empP` at higher compute cost.
- **DS3 sensitivity:** `delta_r_thresh` (default 0.10). Larger = more conservative.

---






## Citation

If you use this code, please cite the manuscript describing the phenotype-filtered DS framework:

---

Jung-Chen Su, Chen-Ling Lee, Fan-Wei Yang, Yan-Chih Chen, and Te-Lun Mai*
"Network-based exploration of 4-(Phenylsulfonyl)morpholine molecules for metastatic triple-negative breast cancer suppression", under review.

---

