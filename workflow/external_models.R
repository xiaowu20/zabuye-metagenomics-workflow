# Refit the published equal-weight, hierarchical and leave-one-lake-out models.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript external_models.R S4.xlsx|lake_traits.tsv OUTPUT_DIRECTORY")
data <- if (endsWith(tolower(args[1]), ".xlsx")) {
  as.data.frame(readxl::read_excel(args[1], sheet = "QXLSG_lake_MAG_traits"))
} else read.delim(args[1], check.names=FALSE)
data <- data[complete.cases(data[, c("Abbreviation", "log_sal", "median_pI")]), ]
stopifnot(!anyDuplicated(data$mag_id))
data$Abbreviation <- factor(data$Abbreviation)
data$Phylum[is.na(data$Phylum)] <- "Unclassified"
data$Phylum <- factor(data$Phylum)
lakes <- aggregate(median_pI ~ Abbreviation + log_sal, data = data, FUN = mean)
stopifnot(!anyDuplicated(lakes$Abbreviation))
row <- function(model, omission, beta, se, interval, p) {
  data.frame(model, omitted_lake = omission, estimate = unname(beta), se = unname(se),
             conf_low = unname(interval[1]), conf_high = unname(interval[2]), P = unname(p))
}
fit_pair <- function(omission) {
  x <- droplevels(data[as.character(data$Abbreviation) != omission, ])
  y <- lakes[as.character(lakes$Abbreviation) != omission, ]
  ols <- lm(median_pI ~ log_sal, data = y)
  co <- coef(summary(ols))["log_sal", ]
  lake <- nlme::lme(median_pI ~ log_sal, random = ~ 1 | Abbreviation, data = x,
                    method = "ML", control = nlme::lmeControl(opt = "optim", msMaxIter = 200))
  lc <- summary(lake)$tTable["log_sal", ]
  ci <- nlme::intervals(lake, which = "fixed")$fixed["log_sal", c("lower", "upper")]
  rbind(row("equal_weight_OLS", omission, co[1], co[2], confint(ols)["log_sal", ], co[4]),
        row("lake_LMM", omission, lc["Value"], lc["Std.Error"], ci, lc["p-value"]))
}
rows <- do.call(rbind, lapply(c("all", levels(data$Abbreviation)), fit_pair))
cross <- lme4::lmer(median_pI ~ log_sal + (1 | Abbreviation) + (1 | Phylum), data = data,
                    REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))
null <- update(cross, . ~ . - log_sal)
lrt <- anova(null, cross)
cross_co <- coef(summary(cross))["log_sal", ]
rows <- rbind(rows, row("lake_phylum_LMM", "all", cross_co[1], cross_co[2],
                      confint(cross, parm = "log_sal", method = "Wald"), lrt[2, "Pr(>Chisq)"]))
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
write.table(rows, file.path(args[2], "lake_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
