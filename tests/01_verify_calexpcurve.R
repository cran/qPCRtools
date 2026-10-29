# tests/01_verify_calexpcurve.R
#
# 验证 CalExpCurve() 的两处关键口径（导师 2026-10-14 审稿意见第 5 条）：
#   ① 反对数：表达量应为**浓度** = dilution^((Cq - Intercept)/Slope)；
#      原实现直接把 (Cq - Intercept)/Slope（对数浓度）拿去除以参考基因的对数值，两对数相除无意义。
#   ② 技术重复层级：技术孔只用来估计同一生物学重复的测量值；组间检验必须以**生物学重复**为单位。
#      否则同一 BioRep 的多个技术孔被当成独立样本（伪重复），p 值被算小。
#      本脚本的判别用例：按生物学重复 p = 0.074（NS），按孔 p = 1.2e-4（***）。
#
# 运行：在包根目录 `Rscript tests/01_verify_calexpcurve.R`
# 也随 `R CMD check` 自动运行（失败即 check 失败）。

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(rstatix)
  library(multcomp)
  library(tibble)
  library(magrittr)
})

# 已安装包时用包内函数；否则直接 source 源文件（便于未安装全部依赖时本地跑）
if (requireNamespace("qPCRtools", quietly = TRUE)) {
  CalExpCurve <- qPCRtools::CalExpCurve
} else {
  r.file <- if (file.exists(file.path("R", "CalExpCurve.R"))) {
    file.path("R", "CalExpCurve.R")
  } else {
    file.path("..", "R", "CalExpCurve.R")
  }
  source(r.file)
}

stopifnot(is.function(CalExpCurve))

curve.table <- data.frame(
  Gene = "G1", Slope = -2, Intercept = 20, min.Cq = 10, max.Cq = 40
)

# --- 判别用例：CK 三个生物学重复的浓度 9/10/11，T 三个重复都是 12 -----------------
# 每个生物学重复 3 个技术孔，孔内完全一致（技术重复只提供测量精度，不提供重复数）
expre.ck <- c(9, 10, 11)
expre.t <- c(12, 12, 12)
cq.of <- function(expre) 20 - 2 * log(expre, base = 4)

build.inputs <- function() {
  cq <- data.frame(Position = character(0), Cq = numeric(0), stringsAsFactors = FALSE)
  design <- data.frame(Position = character(0), Treatment = character(0),
                       Gene = character(0), BioRep = character(0),
                       stringsAsFactors = FALSE)
  prefix <- c("A", "B")
  values <- list(expre.ck, expre.t)
  groups <- c("CK", "T")
  for (g in seq_along(groups)) {
    for (b in seq_along(values[[g]])) {
      for (t in 1:3) {
        position <- paste0(prefix[g], b, t)
        cq <- rbind(cq, data.frame(Position = position, Cq = cq.of(values[[g]][b]),
                                   stringsAsFactors = FALSE))
        design <- rbind(design, data.frame(Position = position, Treatment = groups[g],
                                           Gene = "G1", BioRep = as.character(b),
                                           stringsAsFactors = FALSE))
      }
    }
  }
  list(cq = cq, design = design)
}

inputs <- build.inputs()
res <- CalExpCurve(inputs$cq, curve.table, inputs$design,
                   correction = FALSE, stat.method = "t.test", ref.group = "CK")

tb <- res[["table"]]

# ① 反对数：CK 的生物学重复应还原成浓度 9/10/11，T 的应还原成 12
ck.values <- sort(round(tb$expre[tb$Treatment == "CK"], 6))
t.values <- sort(round(tb$expre[tb$Treatment == "T"], 6))
stopifnot(
  "expression must be back-transformed to concentration" =
    all(abs(ck.values - c(9, 10, 11)) < 1e-6),
  "expression must be back-transformed to concentration (T)" =
    all(abs(t.values - c(12, 12, 12)) < 1e-6)
)

# ② 技术重复层级：每个 (Gene, Treatment, BioRep) 只剩一行（3 个技术孔合并）
stopifnot(
  "technical replicates must be averaged within BioRep" = nrow(tb) == 6
)

# ③ 组间检验以生物学重复为单位：p = 0.074 → NS（按孔计算会得到 p = 1.2e-4 → ***）
#    参考组自身的行没有 p 值（rstatix 只给"其他组 vs 参考组"），故只检查非 NA 的行。
stopifnot(
  "reference group rows carry no p-value" = all(is.na(tb$signif[tb$Treatment == "CK"])),
  "between-group tests must use biological replicates" =
    all(tb$signif[!is.na(tb$signif)] == "NS"),
  "treatment rows must be labelled" =
    all(!is.na(tb$signif[tb$Treatment == "T"]))
)

# --- 缺 BioRep 时必须明确告警（不能静默按孔当独立样本） --------------------------
no.biorep <- inputs$design %>% dplyr::select(-BioRep)
warned <- FALSE
withCallingHandlers(
  CalExpCurve(inputs$cq, curve.table, no.biorep,
              correction = FALSE, stat.method = "t.test", ref.group = "CK"),
  warning = function(w) {
    if (grepl("BioRep", conditionMessage(w))) warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)
stopifnot("missing BioRep must warn about pseudoreplication" = warned)

cat("OK: CalExpCurve back-transforms to concentration and aggregates technical replicates within BioRep\n")
