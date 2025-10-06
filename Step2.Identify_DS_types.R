# ===============================================================================
# Classify S / NS / NA per condition, then call DS1/DS2/DS3 (parametric labels)
#   Input : file with columns V_/Cor_/P_/empP_ + <CellLine>_<TAG>
#           Example: MDAMB231_DMSO / MDAMB231_GL24
#   Output: DS_types.tsv
#
# Purpose:
#   Given per-group correlation statistics (from the previous script), assign
#   per-condition states S / NS / NA for every gene pair, derive transition
#   labels between before/after within each cell line, and then classify pairs
#   into DS1 (gain), DS2 (loss), or DS3 (S->S with |Δr| change), with optional
#   constraints using a non-effective control line.
#
# Key definitions (as implemented here):
#   - S   : significant only when Both P < 0.005 AND empP < 0.005 are finite.
#   - NS  : non-significant when P > 0.05 (empP not required for NS).
#   - NA  : any other case (including V==0, missing/indeterminate tests).
#   - Transitions within a line (before_tag→after_tag):
#       * (NS/NA) -> S   : "gain"
#       * S -> (NS/NA)   : "loss"
#       * All other combos (NS/NA↔NS/NA, S->S) : "unchanged"
#
# Notes:
#   - DS3 uses a magnitude change on |r|: abs(|r_after| - |r_before|).
#     This treats sign flips symmetrically by operating on |r|.
#   - Control usage is configurable via `require_NE_unchanged`:
#       TRUE  -> controls must be strictly "unchanged"
#       FALSE -> controls must simply avoid the same-direction change
#               (no "gain" for DS1; no "loss" for DS2)
#   - This script assumes that each line has both <before_tag> and <after_tag>
#     columns present in the input file for V_/Cor_/P_/empP_.
#   - `Phenotpye` is a typo kept as-is in a variable comment to avoid code edits.
# ===============================================================================

library(data.table)
library(stringr)

# ----------- I/O -----------
infile <- "DS_pairwise_cor_stats.tsv"   # Input: output of the previous step
out_ds <- "DS_types.tsv"                # Output: DS calls + per-line transitions/states

# ----------- Parameters -----------
before_tag <- "DMSO"    # "before" condition suffix in column names
after_tag  <- "GL24"    # "after"  condition suffix in column names

# Phenotype A (effective lines must satisfy DS logic in ALL listed lines)
effective_lines    <- c("MDAMB231", "MDAMB157") # e.g., "MDAMB231", "MDAMB157"

# Phenotpye B (non-effective control lines; may be character(0) to disable)
noneffective_lines <- c("Hs578T")

delta_r_thresh <- 0.10         # DS3 threshold on |Δr| = ||r_after| - |r_before||
require_NE_unchanged <- TRUE   # If TRUE, control lines must be "unchanged"

# ----------- Read -----------
dt <- fread(infile)
id_cols <- names(dt)[1:4]  # Assumes first 4 columns are identifiers (NumA, NumB, GeneA, GeneB)

rx_escape <- function(x) gsub("([][{}()+*^$|?.\\-])", "\\\\\\1", x)
sanitize_lines <- function(x) unique(x[is.finite(nchar(x)) & nzchar(x)])

before_tag_rx <- rx_escape(before_tag)
after_tag_rx  <- rx_escape(after_tag)

# ----------- Infer <Line> Name -----------
# Extract condition identifiers of the form "<Line>_<TAG>" from V_/Cor_/P_/empP_ columns
cond_cols <- grep("^(V_|Cor_|P_|empP_)", names(dt), value = TRUE)
conds <- unique(gsub("^(V_|Cor_|P_|empP_)", "", cond_cols))

# Derive all <Line> names that have BOTH before and after tags
lines <- unique(gsub(paste0("_(?:", before_tag_rx, "|", after_tag_rx, ")$"), "", conds, perl=TRUE))
has_pair <- function(ln) all(paste0(ln, c(paste0("_", before_tag), paste0("_", after_tag))) %in% conds)
stopifnot(all(vapply(lines, has_pair, logical(1))))  # ensure each line has before & after

# Check/clean requested line sets against inferred lines
effective_lines    <- sanitize_lines(effective_lines)
noneffective_lines <- sanitize_lines(noneffective_lines)
if (!require_NE_unchanged) noneffective_lines <- character(0)   # No-control mode when relaxing NE

stopifnot(all(effective_lines %in% lines))
if (length(noneffective_lines) > 0) stopifnot(all(noneffective_lines %in% lines))

# Ensure required per-condition columns exist (V_, P_, empP_; Cor_ is optional downstream)
pick <- function(prefix, cond) paste0(prefix, cond)
need_ok <- all(unlist(lapply(conds, function(cn){
  all(c(pick("V_",cn), pick("P_",cn), pick("empP_",cn)) %in% names(dt))
})))
stopifnot(need_ok)

# ----------- S / NS / NA Types -----------
# State labeling rule:
#   NA     : V is NA or 0, or ambiguous (P ≤ 0.05 but not < 0.005 with empP < 0.005)
#   S      : P < 0.005 & empP < 0.005 (both finite)
#   NS     : P > 0.05 (finite), regardless of empP
#   (Else) : NA
lab_one_cond <- function(V, P, empP){
  if (is.na(V) || V==0) return("NA")
  if (is.finite(P) && is.finite(empP) && P < 0.005 && empP < 0.005) return("S")
  if (is.finite(P) && P > 0.05) return("NS")
  return("NA")
}

# Long format assembly: one row per pair × condition with V/Cor/P/empP and derived state
rows <- vector("list", nrow(dt)*length(conds))
k <- 0L
for (i in seq_len(nrow(dt))) for (cn in conds) {
  k <- k + 1L
  rows[[k]] <- data.table(
    NumA = dt$NumA[i], NumB = dt$NumB[i],
    GeneA = dt$GeneA[i], GeneB = dt$GeneB[i],
    condition = cn,
    V    = dt[[pick("V_",cn)]][i],
    Cor  = if (pick("Cor_",cn) %in% names(dt)) dt[[pick("Cor_",cn)]][i] else NA_real_,
    P    = dt[[pick("P_",cn)]][i],
    empP = dt[[pick("empP_",cn)]][i]
  )
}
long <- rbindlist(rows, use.names = TRUE)
long[, state := mapply(lab_one_cond, V, P, empP)]

# Wide views for states and correlations
state_wide <- dcast(long[, .(NumA,NumB,GeneA,GeneB,condition,state)],
                    NumA + NumB + GeneA + GeneB ~ condition, value.var = "state")
cor_wide   <- dcast(long[, .(NumA,NumB,GeneA,GeneB,condition,Cor)],
                    NumA + NumB + GeneA + GeneB ~ condition, value.var = "Cor")

# ----------- Transition -----------
# Define per-line transitions between before/after using the state_wide table.
#   NS/NA↔NS/NA → "unchanged"
#   (NS/NA)->S  → "gain"
#   S->(NS/NA)  → "loss"
#   S->S        → "unchanged" here; DS3 distinguishes "changed"/"unchanged" through |Δr|
trans_one_vec <- function(bef, aft){
  in_nons_bef <- bef %in% c("NS","NA")
  in_nons_aft <- aft %in% c("NS","NA")
  out <- rep("unchanged", length(bef))
  out[in_nons_bef & !in_nons_aft & aft=="S"] <- "gain"
  out[!in_nons_bef & bef=="S" & in_nons_aft] <- "loss"
  out
}

trans <- copy(state_wide)
for (ln in lines){
  bcol <- paste0(ln, "_", before_tag)
  acol <- paste0(ln, "_", after_tag)
  tcol <- paste0("T_", ln)
  trans[[tcol]] <- trans_one_vec(trans[[bcol]], trans[[acol]])
}

# ----------- DS3：S->S changed / unchanged -----------
# Use absolute-r change: abs(|r_after| - |r_before|) ≥ delta_r_thresh → "changed"
# Only defined when both states are S (S->S); otherwise left as NA.
delta_state <- function(r_before, r_after, thr){
  if (!is.finite(r_before) || !is.finite(r_after)) return(NA_character_)
  d <- abs(abs(r_after) - abs(r_before))
  if (d >= thr) "changed" else "unchanged"
}

# Join state and correlation views to compute per-line S->S change calls (C_<line>)
tw <- merge(state_wide, cor_wide, by = c("NumA","NumB","GeneA","GeneB"),
            suffixes = c("", "_cor"), all.x = TRUE)

for (ln in lines){
  sb <- paste0(ln, "_", before_tag)
  sa <- paste0(ln, "_", after_tag)
  rb <- paste0(ln, "_", before_tag, "_cor")
  ra <- paste0(ln, "_", after_tag,  "_cor")
  ch <- paste0("C_", ln)  # C_<line> is meaningful only when S->S
  tw[[ch]] <- NA_character_
  idx_ss <- which(tw[[sb]]=="S" & tw[[sa]]=="S")
  if (length(idx_ss) > 0){
    tw[[ch]][idx_ss] <- mapply(delta_state, tw[[rb]][idx_ss], tw[[ra]][idx_ss],
                               MoreArgs = list(thr = delta_r_thresh))
  }
}

# all_equal_cols: row-wise test that ALL requested columns equal a target value
all_equal_cols <- function(DT, cols, val){
  cols <- intersect(cols, names(DT))          
  if (length(cols)==0) return(rep(TRUE, nrow(DT)))
  DT[, rowSums(.SD == val) == length(cols), .SDcols = cols]
}
# none_equal_cols: row-wise test that NONE of the requested columns equal a "bad" value
none_equal_cols <- function(DT, cols, bad){
  cols <- intersect(cols, names(DT))        
  if (length(cols)==0) return(rep(TRUE, nrow(DT)))
  DT[, rowSums(.SD == bad) == 0, .SDcols = cols]
}

# ----------- DS1 / DS2 -----------
# Effective lines must ALL show "gain" (DS1) or "loss" (DS2).
# Control logic:
#   - If no control columns are present: NE masks default to TRUE (no restriction).
#   - If require_NE_unchanged == TRUE: all controls must be "unchanged".
#   - Else (relaxed): controls must not show the same-direction change (no "gain" for DS1; no "loss" for DS2).
T_eff_cols <- paste0("T_", effective_lines)
T_ne_cols_all <- paste0("T_", noneffective_lines)
T_ne_cols_valid <- intersect(T_ne_cols_all, names(trans))   

mask_eff_gain <- all_equal_cols(trans, T_eff_cols, "gain")
mask_eff_loss <- all_equal_cols(trans, T_eff_cols, "loss")

if (length(T_ne_cols_valid) == 0) {
  # if no control
  mask_ne_for_gain <- rep(TRUE, nrow(trans))
  mask_ne_for_loss <- rep(TRUE, nrow(trans))
} else if (require_NE_unchanged) {
  mask_ne_for_gain <- all_equal_cols(trans, T_ne_cols_valid, "unchanged")
  mask_ne_for_loss <- all_equal_cols(trans, T_ne_cols_valid, "unchanged")
} else {
  mask_ne_for_gain <- none_equal_cols(trans, T_ne_cols_valid, "gain")
  mask_ne_for_loss <- none_equal_cols(trans, T_ne_cols_valid, "loss")
}

hit_gain <- trans[mask_eff_gain & mask_ne_for_gain]
hit_loss <- trans[mask_eff_loss & mask_ne_for_loss]

# ----------- DS3 -----------
# Effective lines: must be S->S with "changed" per C_<line>.
eff_ss_changed_list <- lapply(effective_lines, function(ln){
  sb <- paste0(ln, "_", before_tag); sa <- paste0(ln, "_", after_tag); cl <- paste0("C_", ln)
  (state_wide[[sb]]=="S" & state_wide[[sa]]=="S" & tw[[cl]]=="changed")
})
mask_eff_ds3 <- if (length(eff_ss_changed_list)==0) {rep(TRUE, nrow(state_wide))
} else {Reduce(`&`, eff_ss_changed_list)}

# Control lines:
#   A pair is acceptable for DS3 if, for every control line:
#     - it is (NS/NA ↔ NS/NA), OR
#     - it is S->S but "unchanged" by |Δr| (i.e., below delta_r_thresh).
valid_ne_lines <- noneffective_lines[noneffective_lines %in% lines]
if (length(valid_ne_lines) == 0) {
  mask_ne_ds3 <- rep(TRUE, nrow(state_wide))
} else {
  ne_unchanged_ds3_list <- lapply(valid_ne_lines, function(ln){
    sb <- paste0(ln,"_", before_tag); sa <- paste0(ln,"_", after_tag); cl <- paste0("C_", ln)
    ((state_wide[[sb]] %in% c("NS","NA") & state_wide[[sa]] %in% c("NS","NA")) |
        (state_wide[[sb]]=="S" & state_wide[[sa]]=="S" & tw[[cl]]=="unchanged"))
  })
  mask_ne_ds3 <- Reduce(`&`, ne_unchanged_ds3_list)
}

hit_ds3_mask <- mask_eff_ds3 & mask_ne_ds3
hit_ds3 <- state_wide[hit_ds3_mask]

# ----------- Summary & Output -----------
T_cols <- paste0("T_", lines)
C_cols <- paste0("C_", lines)

full_T <- Reduce(function(a,b) merge(a, b, by = id_cols, all = TRUE),
                 list(state_wide, trans[, c(id_cols, T_cols), with = FALSE]))

full_TC <- Reduce(function(a,b) merge(a, b, by = id_cols, all = TRUE),
                  list(full_T,     tw[, c(id_cols, C_cols), with = FALSE]))

ds3_view <- copy(full_TC)
for (ln in lines){
  tcol <- paste0("T_", ln)
  ccol <- paste0("C_", ln)
  ds3_view[[tcol]] <- ifelse(!is.na(ds3_view[[ccol]]), ds3_view[[ccol]], ds3_view[[tcol]])
}
ds3_view[, (C_cols) := NULL]

show_cols <- unique(c(
  id_cols,
  unlist(lapply(lines, function(ln) paste0(ln, c(paste0("_", before_tag), paste0("_", after_tag))))),
  T_cols
))

# Summary table (counts per DS)
summary_dt <- rbindlist(list(
  data.table(DS_type = "DS1 (gain)",   n = nrow(hit_gain)),
  data.table(DS_type = "DS2 (loss)",   n = nrow(hit_loss)),
  data.table(DS_type = "DS3 (change)", n = nrow(hit_ds3))
), use.names = TRUE, fill = TRUE)

hits <- rbindlist(list(
  cbind(data.table(DS_type="DS1 (gain)"),
        merge(hit_gain[, id_cols, with=FALSE], full_T[, ..show_cols], by=id_cols, all.x=TRUE)),
  cbind(data.table(DS_type="DS2 (loss)"),
        merge(hit_loss[, id_cols, with=FALSE], full_T[, ..show_cols], by=id_cols, all.x=TRUE)),
  cbind(data.table(DS_type="DS3 (change)"),
        merge(hit_ds3[, id_cols, with=FALSE], ds3_view[, ..show_cols], by=id_cols, all.x=TRUE))
), use.names=TRUE, fill=TRUE)

# Output（remove NumA / NumB）— keeps GeneA/GeneB labels plus per-line annotations
fwrite(hits[,-c("NumA","NumB")], out_ds, sep = "\t")

# ============================== END =========================================== #
