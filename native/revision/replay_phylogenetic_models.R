suppressPackageStartupMessages({ library(ape); library(phangorn); library(phytools); library(tibble) })
root <- commandArgs(trailingOnly=TRUE)[1]
out <- file.path(root, "results/statistics")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.delim(file.path(root, "inputs/Table_S25_ZBY_pI_multivariable_input_corrected.tsv"), check.names=FALSE)
tree <- phangorn::midpoint(ape::read.tree(file.path(root, "inputs/ar53_iqtree.treefile")))
tree$edge.length <- pmax(tree$edge.length, 1e-8)
halo <- d[d$gtdb_phylum == "Halobacteriota", ]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, halo$mag_id))
halo <- halo[match(tree$tip.label, halo$mag_id), ]
stopifnot(nrow(halo)==84, identical(as.character(halo$mag_id), tree$tip.label))
fit_manual_pgls <- function(formula, data, phy) {
  data <- as.data.frame(data)
  rownames(data) <- data$mag_id
  data <- data[phy$tip.label, , drop = FALSE]
  X <- model.matrix(formula, data = data)
  y <- model.response(model.frame(formula, data = data))
  C <- vcv.phylo(phy, model = "Brownian", corr = FALSE)
  C <- C[rownames(data), rownames(data)]
  n <- length(y)
  p <- ncol(X)
  evaluate <- function(lambda, return_fit = FALSE) {
    V <- lambda * C
    diag(V) <- diag(C)
    V <- V + diag(1e-8, n)
    Vi <- solve(V)
    XtViX <- crossprod(X, Vi %*% X)
    beta <- solve(XtViX, crossprod(X, Vi %*% y))
    resid <- y - as.vector(X %*% beta)
    rss <- as.numeric(crossprod(resid, Vi %*% resid))
    sigma2 <- rss / n
    logdet <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
    nll <- 0.5 * (n * log(2 * pi) + n * log(sigma2) + logdet + n)
    if (!return_fit) return(nll)
    vc_beta <- sigma2 * solve(XtViX)
    se <- sqrt(diag(vc_beta))
    t_value <- as.vector(beta) / se
    p_value <- 2 * pt(abs(t_value), df = n - p, lower.tail = FALSE)
    tibble(
      term = colnames(X), estimate = as.vector(beta), se = se,
      statistic = t_value, p_value = p_value,
      conf_low = as.vector(beta) - qt(0.975, n - p) * se,
      conf_high = as.vector(beta) + qt(0.975, n - p) * se,
      n = n, lambda = lambda, sigma2 = sigma2, logLik = -nll
    )
  }
  opt <- optimize(evaluate, interval = c(0, 1), tol = 1e-6)
  evaluate(opt$minimum, return_fit = TRUE)
}

set.seed(20261110)
fit <- fit_manual_pgls(z_pI ~ z_salt + z_completeness + z_contamination + z_log_n_proteins, halo, tree)
signal <- phytools::phylosig(tree, setNames(halo$median_pI, halo$mag_id), method="lambda", test=TRUE)
write.table(fit, file.path(out,"archival_PGLS_replay.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(n=84, lambda=signal$lambda, logL=signal$logL, logL0=signal$logL0,
 LRT=2*(signal$logL-signal$logL0), p_chisq1=signal$P), file.path(out,"archival_signal_replay.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
ape::write.tree(tree, file.path(out,"archival_midpoint_halo84.nwk"))
C <- ape::vcv.phylo(tree, model="Brownian", corr=FALSE)
write.table(C, file.path(out,"archival_Brownian_covariance.tsv"), sep="\t", quote=FALSE, col.names=NA)
# Separate finite-residual-df sensitivity; do not relabel the archived ML-scale intervals.
f <- as.data.frame(fit)
f$se <- f$se*sqrt(84/79)
f$statistic <- f$estimate/f$se
f$p_value <- 2*pt(abs(f$statistic),df=79,lower.tail=FALSE)
f$conf_low <- f$estimate-qt(.975,79)*f$se
f$conf_high <- f$estimate+qt(.975,79)*f$se
write.table(f,file.path(out,"archival_PGLS_residual_df_sensitivity.tsv"), sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"replay_sessionInfo.txt"))
capture.output(phytools::phylosig,file=file.path(out,"phylosig_implementation.txt"))
saveRDS(list(data=halo,tree=tree,C=C,signal=signal,fit=fit),file.path(out,"archival_phylogenetic_replay.rds"))
print(signal);print(fit)
