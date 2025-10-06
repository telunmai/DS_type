# =========================================================
# Build DS input file (V_/Cor_/P_/empP_ + <CellLine>_<TAG>)
# Outputs:
#   - DS_pairwise_cor_stats.tsv
#
# Description:
#   For each cell-line/treatment group, compute Pearson correlation for every
#   gene pair and a permutation-based empirical P-value. Results are emitted
#   as a wide table with per-group columns:
#     V_<group>     : indicator (1) that a valid cor.test() was computed; 0 otherwise
#     Cor_<group>   : Pearson correlation coefficient (r), NA if not computed
#     P_<group>     : uncorrected P-value from cor.test(), NA if not computed
#     empP_<group>  : empirical P = Pr(P_perm < P_obs) over B permutations, NA if not computed
#
#   This file is the input for downstream DS classification (DS1/DS2/DS3) used
#   in the manuscript’s phenotype-filtered differential co-expression framework.
#
# Reproducibility note:
#   - set.seed(b) is called inside the permutation loop. This makes the permutation stream
#     deterministic AND shared across all gene pairs for a given b (and group). If you prefer
#     independent permutation streams per pair, set the seed once globally before the outer loops
#
# Input format:
#   - infile is a tabular matrix: rows = genes, columns = samples; one column contains gene IDs.
#     IMPORTANT: gene_col is 1-based and is removed before correlation; mapping accounts for that.
#
# Output schema:
#   Columns:
#     NumA, NumB : integer indices of the two genes in the upper-tri enumeration
#     GeneA, GeneB : gene IDs from 'genes'
#     For each group G:
#       V_G, Cor_G, P_G, empP_G as defined above
#
# Edge cases handled:
#   - All-zero expression in either gene within a group => NA results, NAflag=TRUE internally
#   - Fewer than 3 finite points OR < 3 unique (x,y) pairs => NA results
#   - cor.test() failure or non-finite outputs => NA results
#   - Any permutation cor.test() failure => empP set to NA for that pair/group
#
# =========================================================

library(data.table)

# ===================== 0) User Setting =====================
infile   <- "TPM_matrix_sample.txt"            # Gene x Samples (including geneID)
gene_col <- 1                                  # Location of GeneID（1-based）
ngene    <- NA                                 # Number of genes; NA = all genes

group_cols <- list(
  "Hs578T_DMSO"      = "^Hs578T_DMSO-",
  "Hs578T_GL24"      = "^Hs578T_GL24-",
  "MDAMB157_DMSO"    = "^MDAMB157_DMSO-",
  "MDAMB157_GL24"    = "^MDAMB157_GL24-",
  "MDAMB231_DMSO"    = "^MDAMB231_DMSO-",
  "MDAMB231_GL24"    = "^MDAMB231_GL24-"
)
# Tip: Each group is expected to have exactly 5 replicates (checked below).

# number of permutation for empirical P 
B_perm <- 200  # Empirical null size per group (computationally significant for large G)

# ===================== 1) pre-processing =====================
Expr <- fread(infile)
if (!is.na(ngene)) Expr <- Expr[1:ngene, ]
stopifnot(gene_col >= 1, gene_col <= ncol(Expr))

genes    <- Expr[[gene_col]]
expr_mat <- as.matrix(Expr[, -gene_col, with = FALSE])   # removing geneID
colmap   <- data.table(orig_name = colnames(Expr), orig_pos = seq_len(ncol(Expr)))
colmap   <- colmap[orig_pos != gene_col]                
colmap[, expr_pos := seq_len(.N)]

resolve_group_to_expr_idx <- function(gdef) {
  if (is.numeric(gdef)) {
    pos <- as.integer(gdef)
    if (any(pos < 1 | pos > ncol(Expr))) {
      stop("Out of Index", paste(pos, collapse = ", "))
    }
  } else if (is.character(gdef)) {
    if (length(gdef) == 1L) {
      hits <- grep(gdef, colnames(Expr), value = FALSE)
      if (length(hits) == 0L) {
        stop("Can not identify the column：", gdef)
      }
      pos <- hits
    } else {
      hits <- match(gdef, colnames(Expr))
      if (any(is.na(hits))) {
        stop("Not found:", paste(gdef[is.na(hits)], collapse = ", "))
      }
      pos <- hits
    }
  } else {
    stop("group_cols should be numeric or character。")
  }
  # Adjust positions because gene_col was removed to form expr_mat:
  expr_idx <- pos - (pos > gene_col)
  if (any(expr_idx < 1 | expr_idx > ncol(expr_mat))) {
    stop("Out of Index", paste(expr_idx, collapse = ", "))
  }
  expr_idx
}

group_cols_em <- lapply(group_cols, resolve_group_to_expr_idx)
# Check the number of replicates
stopifnot(all(vapply(group_cols_em, length, 1L) == 5))  # Enforces n=5 per group design

message("== Group → expr_mat ==")
for (nm in names(group_cols_em)) {
  idx <- group_cols_em[[nm]]
  orig_names <- colmap[expr_pos %in% idx, orig_name]
  message(sprintf("%-18s : %s", nm, paste(orig_names, collapse = ", ")))
}

# Correlation and Permutation
pair_stats_one_group <- function(x, y, B = 200) {
  # Compute Pearson correlation (r), P, and permutation empirical P between x and y.
  # x, y are numeric vectors for one gene across the 5 replicates in a single group.
  #
  # Returns a list with:
  #   V      : 1 if cor.test() succeeded and r/P are finite; 0 otherwise
  #   Cor    : Pearson r (NA if not computed)
  #   P      : uncorrected P-value from cor.test() (NA if not computed)
  #   empP   : empirical P = mean(P_perm < P_obs) over B permutations (NA if not computed)
  #   NAflag : TRUE only for the all-zero-expression case; FALSE otherwise
  
  # 1) All values = 0
  if (sum(x > 0) == 0 || sum(y > 0) == 0) {
    return(list(V = 0L, Cor = NA_real_, P = NA_real_, empP = NA_real_, NAflag = TRUE))
  }
  
  # 2) Valided points >= 3 (x,y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(V = 0L, Cor = NA_real_, P = NA_real_, empP = NA_real_, NAflag = FALSE))
  xy <- cbind(x[ok], y[ok])
  if (nrow(unique(xy)) < 3) return(list(V = 0L, Cor = NA_real_, P = NA_real_, empP = NA_real_, NAflag = FALSE))
  
  # 3) Correlation
  ct <- try(suppressWarnings(cor.test(xy[,1], xy[,2], method = "pearson")), silent = TRUE)
  if (inherits(ct, "try-error") || !is.finite(ct$p.value) || !is.finite(ct$estimate)) {
    return(list(V = 1L, Cor = NA_real_, P = NA_real_, empP = NA_real_, NAflag = FALSE))
  }
  r_obs <- unname(ct$estimate)
  p_obs <- unname(ct$p.value)
  
  # 4) empirical P Calculation
  #    Null is generated by independently permuting sample indices of x and y.
  #    Empirical P is the fraction of permuted P-values that are STRICTLY smaller than P_obs.
  #    If any permutation test errors out, empP is set to NA for the pair.
  if (B > 0 && is.finite(p_obs)) {
    n <- nrow(xy)
    count_p <- 0L
    empP <- NA_real_
    for (b in seq_len(B)) {
      set.seed(b)                             # Deterministic across pairs; see note above
      perm1 <- sample.int(n)                  # independent shuffles for x and y
      perm2 <- sample.int(n)                   
      ctb <- try(suppressWarnings(
        cor.test(xy[perm1, 1], xy[perm2, 2], method = "pearson")
      ), silent = TRUE)
      
      if (inherits(ctb, "try-error") || !is.finite(ctb$p.value)) {
        empP <- NA_real_
        break
      }
      
      if (ctb$p.value < p_obs) count_p <- count_p + 1L  # strict inequality as implemented
      if (b == B) empP <- count_p / B
    }
  } else {
    empP <- NA_real_
  }
  
  list(V = 1L, Cor = r_obs, P = p_obs, empP = empP, NAflag = FALSE)
}


# ===================== 3) Main: pair × group =====================
G <- nrow(expr_mat)
if (G < 2) stop("Too few genes.")  # Requires at least 2 genes to form a pair

# Enumerate upper-triangular gene pairs (i < j)
ut_idx  <- which(upper.tri(matrix(NA, G, G)), arr.ind = TRUE)
n_pairs <- nrow(ut_idx)

group_names <- names(group_cols_em)

# Output Columns
out_cols <- c("NumA","NumB","GeneA","GeneB")
for (gn in group_names) out_cols <- c(out_cols, paste0(c("V_","Cor_","P_","empP_"), gn))

# Pre-allocate results; numeric matrix coerced to data.table; strings added below
res <- data.table(matrix(NA_real_, nrow = n_pairs, ncol = length(out_cols)))
setnames(res, out_cols)
res[, NumA := ut_idx[, "row"]]
res[, NumB := ut_idx[, "col"]]
res[, GeneA := genes[ut_idx[, "row"]]]
res[, GeneB := genes[ut_idx[, "col"]]]

message(sprintf("Total pairs: %d; groups: %d", n_pairs, length(group_names)))

# Core loop over pairs and groups (computational hot spot)
for (idx in seq_len(n_pairs)) {
  if (idx %% 1000 == 0) message(sprintf(".. processed %d / %d pairs", idx, n_pairs))
  i <- ut_idx[idx, "row"]; j <- ut_idx[idx, "col"]
  
  for (gn in group_names) {
    cols <- group_cols_em[[gn]]
    xi <- as.numeric(expr_mat[i, cols])
    yj <- as.numeric(expr_mat[j, cols])
    ps <- pair_stats_one_group(xi, yj, B = B_perm)
    
    res[idx, paste0("V_",    gn) := ps$V]
    res[idx, paste0("Cor_",  gn) := ps$Cor]
    res[idx, paste0("P_",    gn) := ps$P]
    res[idx, paste0("empP_", gn) := ps$empP]
  }
}

# ===================== 4) Output =====================
fwrite(res, file = "DS_pairwise_cor_stats.tsv", sep = "\t", na = "NA",quote = FALSE)

message("Done. Files written: DS_pairwise_cor_stats.tsv")
