library(qPCRtools)

cq <- data.frame(
  Position = paste0("A", 1:12),
  Cq = c(22, 20, 22.2, 20.1, 21.8, 19.9,
         21, 20, 21.1, 20.0, 20.9, 19.9)
)
design <- data.frame(
  Position = paste0("A", 1:12),
  Group = rep(c("CK", "CK", "CK", "T", "T", "T"), each = 2),
  Gene = rep(c("G1", "REF"), 6),
  BioRep = rep(1:3, each = 2, times = 2)
)

ddct <- CalExp2ddCt(
  cq, design,
  ref.gene = "REF",
  ref.group = "CK",
  stat.method = "t.test",
  remove.outliers = FALSE,
  fig.type = "box"
)
tb <- ddct[["table"]]
ck.mean <- mean(tb$mean.expre[tb$Treatment == "CK"])
stopifnot(abs(ck.mean - 1) < 0.05)
stopifnot(all(abs(tb$se.expre - tb$sd.expre / sqrt(tb$n)) < 1e-12))
stopifnot(all(c("RawPValue", "PValue", "PAdjustMethod") %in% names(tb)))
stopifnot(is.finite(tb$PValue[tb$Treatment == "T"]))
stopifnot(unique(tb$PAdjustMethod[tb$Treatment == "T"]) == "none")

text.layers <- vapply(
  ddct[["figure"]]$layers,
  function(layer) inherits(layer$geom, "GeomText"),
  logical(1)
)
stopifnot(sum(text.layers) == 1)

bar <- CalExp2ddCt(
  cq, design,
  ref.gene = "REF",
  ref.group = "CK",
  remove.outliers = FALSE,
  fig.type = "bar"
)
expected <- unique(bar[["figure"]]$data[, c("gene", "Treatment", "mean.expre")])
expected <- expected[order(expected$gene, expected$Treatment), ]
bar.data <- ggplot2::ggplot_build(bar[["figure"]])$data[[1]]
stopifnot(max(abs(sort(bar.data$y) - sort(expected$mean.expre))) < 1e-8)

anova.design <- expand.grid(
  BioRep = 1:3,
  Group = c("CK", "T1", "T2"),
  Gene = c("G1", "REF"),
  stringsAsFactors = FALSE
)
anova.cq <- data.frame(
  Position = paste0("B", seq_len(nrow(anova.design))),
  Cq = ifelse(
    anova.design$Gene == "REF",
    20,
    22 + 2 * as.numeric(factor(anova.design$Group))
  )
)
anova.design$Position <- anova.cq$Position
anova.res <- CalExp2ddCt(
  anova.cq, anova.design,
  ref.gene = "REF",
  ref.group = "CK",
  stat.method = "anova",
  remove.outliers = FALSE,
  fig.type = "box"
)
anova.tb <- anova.res[["table"]]
nonref <- anova.tb$Treatment != "CK"
stopifnot(all(c("RawPValue", "PValue", "PAdjustMethod") %in% names(anova.tb)))
stopifnot(all(is.finite(anova.tb$RawPValue[nonref])))
stopifnot(all(is.finite(anova.tb$PValue[nonref])))
stopifnot(all(anova.tb$PValue[nonref] >= anova.tb$RawPValue[nonref] - 1e-12))
stopifnot(unique(anova.tb$PAdjustMethod[nonref]) == "Tukey")
invisible(ggplot2::ggplot_build(anova.res[["figure"]]))

rq.cq <- read.table(
  system.file("examples", "cal.expre.rqpcr.cq.txt", package = "qPCRtools"),
  header = TRUE
)
rq.design <- read.table(
  system.file("examples", "cal.expre.rqpcr.design.txt", package = "qPCRtools"),
  header = TRUE
)
rq <- CalExpRqPCR(
  rq.cq, rq.design,
  ref.gene = "OsUBQ",
  ref.group = "CK",
  stat.method = "t.test"
)
stopifnot(nrow(rq[["table"]]) > 0)
stopifnot(any(is.finite(rq[["table"]][["Expression"]])))
stopifnot(all(c("RawPValue", "PValue", "PAdjustMethod") %in% names(rq[["table"]])))

sc.cq <- read.table(
  system.file("examples", "calsc.cq.txt", package = "qPCRtools"),
  header = TRUE
)
sc.info <- read.table(
  system.file("examples", "calsc.info.txt", package = "qPCRtools"),
  header = TRUE
)
curve <- CalCurve(
  sc.cq, sc.info,
  lowest.concen = 4, highest.concen = 4096,
  dilution = 4, by.mean = TRUE
)
deprecated <- FALSE
invisible(withCallingHandlers(
  ggplot2::ggplot_build(curve[["figure"]]),
  warning = function(w) {
    if (grepl("deprecated", conditionMessage(w), ignore.case = TRUE)) {
      deprecated <<- TRUE
    }
    invokeRestart("muffleWarning")
  }
))
stopifnot(!deprecated)

cat("OK: CRAN release regression tests passed\n")
