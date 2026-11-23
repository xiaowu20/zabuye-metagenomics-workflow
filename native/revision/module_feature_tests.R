# Match the archived Figure 3 station-level Welch-studentized enumeration.
# Provenance: Fig3_Deep_Dive_20261020/code/02_analyse.py, test_features and welch_t.
welch_stat <- function(a, b) {
  delta <- colMeans(a) - colMeans(b)
  se2 <- apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b)
  statistic <- rep(0, length(delta))
  ok <- se2 > 1e-28
  statistic[ok] <- delta[ok] / sqrt(se2[ok])
  statistic[!ok & abs(delta) > 1e-14] <- sign(delta[!ok & abs(delta) > 1e-14]) * 1e100
  statistic
}

module_feature_tests <- function(root, group) {
  meta <- read.delim(file.path(root, "inputs/paired14.tsv"))
  output <- file.path(root, "catalogues", group, "statistics")
  effects <- list()
  directions <- list()
  for (cutoff in c(50, 25, 10)) for (response in c("GMP", "SCA", "DCP")) {
    mat <- as.matrix(read.delim(file.path(root, "catalogues", group, "responses",
      paste0(response, "_breadth", cutoff, ".tsv")), row.names = 1, check.names = FALSE))
    mat <- mat[meta$sample, , drop = FALSE]
    sm <- rowsum(mat, meta$station) / 2
    north <- grepl("^N", rownames(sm))
    stopifnot(sum(north) == 4, sum(!north) == 3, ncol(sm) == 49)
    assignments <- combn(nrow(sm), 4, simplify = FALSE)
    obs <- welch_stat(sm[north, , drop = FALSE], sm[!north, , drop = FALSE])
    perm <- sapply(assignments, function(ii) welch_stat(sm[ii, , drop = FALSE], sm[-ii, , drop = FALSE]))
    pval <- rowMeans(abs(perm) >= abs(obs) - 1e-10)
    invariant <- apply(sm, 2, function(x) diff(range(x)) < 1e-14)
    pval[invariant] <- 1
    # All equally likely station resampling sequences within each basin: 4^4 * 3^3.
    bn <- as.matrix(expand.grid(rep(list(which(north)), 4)))
    bs <- as.matrix(expand.grid(rep(list(which(!north)), 3)))
    mn <- t(apply(bn, 1, function(ii) colMeans(sm[ii, , drop = FALSE])))
    ms <- t(apply(bs, 1, function(ii) colMeans(sm[ii, , drop = FALSE])))
    ci <- sapply(seq_len(ncol(sm)), function(j) quantile(as.vector(outer(mn[, j], ms[, j], "-")), c(.025, .975)))
    effects[[length(effects) + 1]] <- data.frame(catalogue = group, response = response,
      breadth = cutoff / 100, module = colnames(sm), north_mean = colMeans(sm[north, , drop = FALSE]),
      south_mean = colMeans(sm[!north, , drop = FALSE]),
      north_minus_south = colMeans(sm[north, , drop = FALSE]) - colMeans(sm[!north, , drop = FALSE]),
      CI_low = ci[1, ], CI_high = ci[2, ], welch_t = obs, P_exact = pval,
      q_BH = p.adjust(pval, "BH"), invariant = invariant, n_assignments = 35, bootstrap_sequences = 6912)
    if (cutoff == 50) {
      scopes <- c("paired14", "spring", "autumn", paste0("omit_", rownames(sm)))
      for (scope in scopes) {
        if (scope %in% c("spring", "autumn")) {
          keep <- meta$season == scope
          z <- rowsum(mat[keep, , drop = FALSE], meta$station[keep])
        } else {
          z <- sm
          if (startsWith(scope, "omit_")) z <- z[rownames(z) != sub("^omit_", "", scope), , drop = FALSE]
        }
        nn <- grepl("^N", rownames(z))
        directions[[length(directions) + 1]] <- data.frame(catalogue = group, response = response,
          module = colnames(z), scope = scope, north_stations = sum(nn), south_stations = sum(!nn),
          north_minus_south = colMeans(z[nn, , drop = FALSE]) - colMeans(z[!nn, , drop = FALSE]))
      }
    }
  }
  tab <- do.call(rbind, effects)
  write.table(tab, file.path(output, "module_basin_contrasts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(do.call(rbind, directions), file.path(output, "module_direction_diagnostics.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  if (group == "archival_original269") {
    old <- read.delim(file.path(root, "inputs/original_feature_tests.tsv"))
    current <- tab[tab$breadth == .5, ]
    key <- paste(current$response, current$module)
    idx <- match(paste(old$family, old$feature), key)
    stopifnot(!anyNA(idx))
    check <- data.frame(family = old$family, module = old$feature,
      effect_difference = abs(current$north_minus_south[idx] - old$delta_NZ_minus_SZ),
      P_difference = abs(current$P_exact[idx] - old$P_exact), q_difference = abs(current$q_BH[idx] - old$q_BH))
    stopifnot(max(check$effect_difference) < 1e-10, max(check$P_difference) < 1e-10, max(check$q_difference) < 1e-10)
    write.table(check, file.path(output, "archived_feature_reproduction_QA.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
