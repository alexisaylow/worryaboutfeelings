## #########################################################################
## #########################################################################
##
##  WORRY IS UNDERPINNED BY PERSISTENCE IN SEARCHING UNINFORMATIVE AVENUES
##  Low, Friston & Moutoussis
##
##  ANALYSIS SCRIPT — reproduces every statistic, table and figure that
##  appears in the manuscript and its Supplement.
##
## -------------------------------------------------------------------------
##  HOW THIS FILE IS ORGANISED
##
##  Sections run top to bottom. Section 1 builds every object the rest of
##  the file needs, so run it first; after that, each analysis section is
##  self-contained and can be run on its own.
##
##  Every section carries a greppable tag. Search (Ctrl-F) for the tag to
##  jump straight to the code behind a given output:
##
##    @SETUP      packages and file paths
##    @DATA       data import, cleaning, derived variables, analysis subsets
##    @HELPERS    small functions used in several places
##
##    @S2.1       Table S2.1 Descriptive statistics (a-d)
##    @S2.2       Table S2.2 Statistical tests (a-d)
##    @S3         Table S3   Trait worry on sampling
##    @S4         Table S4   State worry on sampling
##    @S5         Table S5   Stimuli aversiveness on sampling
##    @S6         Table S6   Stimuli aversiveness on self-reported worry
##    @S7         Table S7   Stimuli uncertainty on sampling
##    @S8         Table S8   Stimuli uncertainty on worry
##    @S9         Table S9   Self-reported uncertainty on sampling
##    @S10        Table S10  Self-reported uncertainty on worry
##    @S11        Table S11  Both uncertainty types on sampling
##    @S12        Table S12  Both uncertainty types on worry
##    @S13        Table S13  Information environment, overall group
##    @S14        Table S14  Information environment, high worriers
##    @S15        Table S15  Information environment, low worriers
##    @S16        Table S16  Information environment x trait worry
##    @S17        Table S17  Stopping and declaring
##    @S18        Table S18  Stopping and declaring, interaction model
##    @S19        Table S19  Subgroup analyses
##
##    @ROPE       Region of practical equivalence analyses (main text)
##    @FIG2       Figure 2, panels a-f and panel assembly
##
## -------------------------------------------------------------------------
##  SOFTWARE
##
##  R 4.3.2. Mixed models for sampling counts use GLMMadaptive 0.9.7 with
##  the negative binomial family; models for self-reported worry and error
##  magnitude use nlme::lme (Gaussian). Estimated marginal means come from
##  emmeans. Bayesian models use rstanarm with default settings (4 chains,
##  2000 iterations, 1000 warmup, weakly-informative auto-scaled normal
##  priors). Figures use ggplot2 4.0.3 with patchwork for panel assembly.
##
##  Negative binomial coefficients exponentiate to incidence rate ratios
##  (IRR); logistic coefficients exponentiate to odds ratios (OR).
##  Confidence intervals are Wald intervals on the link scale, then
##  back-transformed.
##
## #########################################################################
## #########################################################################


## #########################################################################
## @SETUP — PACKAGES AND FILE PATHS
## #########################################################################

## ---- Packages ------------------------------------------------------------
## Install any that are missing, then load. Grouped by what they are for so
## it is clear which analyses depend on which package.

required <- c(
  ## data handling and plotting
  "dplyr", "ggplot2", "ggpubr", "patchwork", "grid",
  ## mixed models
  "nlme", "lme4", "GLMMadaptive", "emmeans",
  ## Bayesian estimation and ROPE
  "rstanarm", "bayestestR", "BayesFactor",
  ## rank-based effect sizes
  "effectsize",
  ## table output
  "knitr"
)

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

## ---- File paths ----------------------------------------------------------
## All three input files sit in the SAME folder as this script, so they are
## referenced by filename only (a relative path). This keeps the script
## portable: anyone who clones the repository can run it without editing
## paths, as long as the CSVs travel with the .R file.
##
## IMPORTANT: R must be pointed at the folder this script lives in, otherwise
## a bare filename will not be found. Two ways to guarantee that:
##   - In RStudio: Session > Set Working Directory > To Source File Location,
##     before running; or open the script via an .Rproj in this folder.
##   - Or set it yourself, e.g. setwd("C:/Users/aylow/Dropbox/PhD/PAPER WRITING")
## The commented line below works in RStudio and sets the working directory
## to wherever this script is saved; uncomment it to make the script
## self-locating rather than relying on the menu step above.
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

file_data         <- "cleanedData.csv"
file_demographics <- "cleanedData_demographics.csv"
file_lookup       <- "corrEmotion_fullset_lookup.csv"

## Centre plot titles globally, as used for all figures in the paper.
theme_update(plot.title = element_text(hjust = 0.5))


## #########################################################################
## @DATA — IMPORT, CLEANING, DERIVED VARIABLES, ANALYSIS SUBSETS
## #########################################################################
##
## The raw export is one row per trial and one block of columns per
## participant (wide). Analysis needs one row per participant-trial (long),
## so the first step reshapes it. Everything downstream uses the long frame
## `cleanedDataANOVA`.
## -------------------------------------------------------------------------

## ---- 1.1 Read and reshape wide to long -----------------------------------

cleanedData <- read.csv(file_data, stringsAsFactors = FALSE)

nTrial <- 32    ## trials per participant
nPt    <- 306   ## participants in the raw export (before worry-group filtering)
nVar   <- 40    ## columns per participant block

## Column names, in the order they appear within each participant block.
hd <- c('openings', 'infoRatio', 'newInfo', 'emotion', 'precision', 'interviewer',
        'filename', 'ptN', 'preoccupied_by_scream', 'uncertainty_before',
        'uncertainty_after', 'worry_before', 'worry_after', 'emotionChosen_before',
        'emotionChosen_after', 'emotionChosenExact_before', 'emotionChosenExact_after',
        'SES', 'BDI', 'PSWQ', 'STAI', 'STAI_state', 'STAI_trait', 'MFQ', 'ICI', 'RSES',
        'perfectionscale', 'perfectionscale_mistakes', 'perfectionscale_parents',
        'perfectionscale_standards', 'perfectionscale_precision', 'mostRecentImage',
        'mostRecentIsNew', 'mostRecentValue', 'combinedImagesList', 'status', 'ptID',
        'worryLevel', 'SOP_before', 'SOP_after')

## Walk down the long frame one row at a time. For each long row, work out
## which participant it belongs to and which trial, then pull that
## participant's column block out of the wide frame.
xwide <- cleanedData
xlong <- data.frame(matrix(NA, nPt * nTrial, nVar))

for (longrow in 1:dim(xlong)[1]) {
  ptN     <- floor((longrow - 1) / nTrial) + 1     ## which participant
  widerow <- longrow %% nTrial                     ## which trial
  if (widerow == 0) widerow <- nTrial              ## last trial wraps to 0
  widecols <- ((ptN - 1) * nVar + 1):(ptN * nVar)  ## that participant's block
  xlong[longrow, ] <- xwide[widerow, widecols]
}

colnames(xlong) <- hd
cleanedDataANOVA <- xlong

## ---- 1.2 Coerce types ----------------------------------------------------
## The reshape above loses column classes, so set the numeric ones back.

num_cols <- c('openings', 'worry_before', 'worry_after',
              'uncertainty_before', 'uncertainty_after',
              'SOP_before', 'SOP_after', 'preoccupied_by_scream',
              'PSWQ', 'STAI', 'STAI_state', 'STAI_trait', 'BDI', 'MFQ',
              'ICI', 'RSES', 'SES',
              'perfectionscale', 'perfectionscale_mistakes',
              'perfectionscale_parents', 'perfectionscale_standards',
              'perfectionscale_precision',
              'emotionChosen_before', 'emotionChosen_after')

cleanedDataANOVA[num_cols] <- lapply(cleanedDataANOVA[num_cols], as.numeric)

## The emotion-choice scale runs 1-9 but is stored 0-8, so shift by one.
cleanedDataANOVA$emotionChosen_before <- cleanedDataANOVA$emotionChosen_before + 1
cleanedDataANOVA$emotionChosen_after  <- cleanedDataANOVA$emotionChosen_after  + 1

## ---- 1.3 Task version and experimental factor labels ---------------------
## Two task versions were run. In v2 the high-precision condition was not
## presented, so every v2 trial is low precision; trials are labelled from
## `filename`, which is the authoritative record of what was shown.

cleanedDataANOVA$taskVersion <- ifelse(grepl("v2", cleanedDataANOVA$filename), "v2", "v1")

## The labels as they came out of the task export are kept alongside the
## analysis labels, so both are available.
cleanedDataANOVA$emotion_asExported   <- cleanedDataANOVA$emotion
cleanedDataANOVA$precision_asExported <- cleanedDataANOVA$precision

v2_highprec <- cleanedDataANOVA$taskVersion == "v2" &
               grepl("highprecision", cleanedDataANOVA$filename)
cleanedDataANOVA$emotion[v2_highprec &
  grepl("happyfroma", cleanedDataANOVA$filename)] <- "angry"
cleanedDataANOVA$precision[v2_highprec] <- "lowPrecision"

cat("\n@DATA  trials by version and aversiveness condition:\n")
print(table(version = cleanedDataANOVA$taskVersion,
            emotion = cleanedDataANOVA$emotion))

## ---- 1.4 Error magnitude (accuracy) --------------------------------------
## Error magnitude is the absolute difference between the correct answer for
## a stimulus set and the participant's final rating. Two definitions of
## "correct" are available and both are constructed here:
##
##   accuracy_hc       the condition's nominal correct answer (6 for the
##                     aversive condition, 4 for the appetitive condition)
##   accuracy_fullset  the mean valence of the full image pool the trial drew
##                     from, i.e. what an ideal observer with access to every
##                     available image would report. Read from a lookup keyed
##                     on `filename`.
##
## accuracy_fullset is the measure reported in the manuscript and behind
## Figure 2d and Table S2.2c.

cleanedDataANOVA$corrEmotion_hc <- ifelse(cleanedDataANOVA$emotion == "angry", 6, 4)

lookup <- read.csv(file_lookup, stringsAsFactors = FALSE)
stopifnot(nrow(lookup) == 64)   ## one row per image pool file

## Normalise the join key on both sides so the match survives differences in
## path prefix, file extension and the optional "initial_" prefix.
normkey <- function(x) {
  x <- tolower(as.character(x))
  x <- sub("^.*(v[12])/", "\\1/", x)   ## keep from v1/ or v2/ onward
  x <- sub("\\.xlsx?$", "", x)         ## drop extension
  sub("/initial_", "/", x)             ## tolerate the initial_ prefix
}

idx <- match(normkey(cleanedDataANOVA$filename), normkey(lookup$facesFile))
if (anyNA(idx)) {
  cat("Unmatched filename values:\n")
  print(head(unique(as.character(cleanedDataANOVA$filename)[is.na(idx)]), 10))
  stop("Join key must match for every row; a partial join would bias error magnitude.")
}
cat("@DATA  image-pool lookup joined for", sum(!is.na(idx)), "of",
    nrow(cleanedDataANOVA), "rows\n")

cleanedDataANOVA$corrEmotion_fullset <- lookup$corrEmotion_fullset[idx]
cleanedDataANOVA$nGreenAvailable     <- lookup$nGreenRows[idx]
cleanedDataANOVA$greenFlagMissing    <- lookup$greenFlagMissing[idx]

cleanedDataANOVA$accuracy_fullset <-
  abs(cleanedDataANOVA$corrEmotion_fullset - cleanedDataANOVA$emotionChosen_after)
cleanedDataANOVA$accuracy_hc <-
  abs(cleanedDataANOVA$corrEmotion_hc - cleanedDataANOVA$emotionChosen_after)

## ---- 1.5 Merge demographics ----------------------------------------------
## Age and gender are held in a separate export with the same row order.

demographics <- read.csv(file_demographics, stringsAsFactors = FALSE)
if (nrow(cleanedDataANOVA) != nrow(demographics))
  stop("Demographics file must have the same number of rows as the trial data.")
cleanedDataANOVA$Age <- demographics$Age
cleanedDataANOVA$Sex <- demographics$Sex

## ---- De-identify the external participant id -----------------------------
## ptID as exported is the Prolific ID, which is identifiable. It is replaced
## here with a plain integer, so that nothing downstream — and no file written
## from this script — carries the real id. factor() groups identical ptID
## values, so the participant grouping is preserved exactly, INCLUDING the one
## person who completed the task twice (a single ptID shared across two ptN).
## Every participant-level result is therefore unchanged. The mapping is not
## stored, so the de-identification is one-way.
##
## NOTE: the raw cleanedData.csv still contains the Prolific IDs. Do not
## deposit that file in a public repository; deposit a copy written from this
## (already de-identified) frame instead, or use a private repository.
cleanedDataANOVA$ptID <- as.integer(factor(cleanedDataANOVA$ptID))

## ---- 1.6 Derived variables and factor levels -----------------------------

## Trial index within participant.
cleanedDataANOVA$trialN <- rep(1:nTrial, times = nPt)

## State worry is rated 1-5. Centring at 1 puts the model intercept at the
## lowest value that actually occurs, rather than at an off-scale zero.
## Centring does not change any slope, SE or p-value.
cleanedDataANOVA$worry_before_c1 <- cleanedDataANOVA$worry_before - 1

## Reference levels. These set which group each coefficient is read against,
## so they determine the sign and interpretation of every table below:
##   infoRatio  "1 to 5"   = information rich   (reference)
##              "1 to 40"  = information poor
##   worryLevel "lowWorry" = reference, so trait-worry terms read as the
##              high-worry group relative to low
cleanedDataANOVA$infoRatio <- factor(cleanedDataANOVA$infoRatio,
                                     levels = c("1 to 5", "1 to 40"))
cleanedDataANOVA$worryLevel <- factor(cleanedDataANOVA$worryLevel,
                                      levels = c("lowWorry", "highWorry", "unknown"))
##   emotion    "happyWithAngry" = reference (less aversive), so the
##              coefficient is named `emotionangry` and reads as the AVERSIVE
##              condition relative to the appetitive one
##   precision  "highPrecision"  = reference (low stimulus variance), so the
##              coefficient is named `precisionlowPrecision` and reads as the
##              HIGH stimulus uncertainty condition relative to low
## Both are set so that the reported coefficient increases with the quantity
## named in the Supplement row ("Stimuli Aversiveness", "Stimuli
## Uncertainty"), rather than with its opposite.
cleanedDataANOVA$emotion   <- factor(cleanedDataANOVA$emotion,
                                     levels = c("happyWithAngry", "angry"))
cleanedDataANOVA$precision <- factor(cleanedDataANOVA$precision,
                                     levels = c("highPrecision", "lowPrecision"))

## Readable labels used by the Figure 2 panels.
cleanedDataANOVA$worryGroup <- factor(
  ifelse(cleanedDataANOVA$worryLevel == "lowWorry", "Low", "High"),
  levels = c("Low", "High"))
cleanedDataANOVA$infoWealth <- factor(
  ifelse(cleanedDataANOVA$infoRatio == "1 to 40", "Poor", "Rich"),
  levels = c("Rich", "Poor"))

## ---- 1.7 Analysis subsets ------------------------------------------------
## Participants screened into neither worry decile are labelled "unknown"
## and are excluded from all group comparisons.

ddR <- cleanedDataANOVA[cleanedDataANOVA$worryLevel != "unknown", ]
ddR$worryLevel <- droplevels(ddR$worryLevel)
ddR$ptN        <- factor(ddR$ptN)

cleanedDataANOVA_highworry <- subset(ddR, worryLevel == "highWorry")
cleanedDataANOVA_lowworry  <- subset(ddR, worryLevel == "lowWorry")

## Self-reported uncertainty: participants who never used the bottom of the
## scale. Used for the self-reported uncertainty models (Tables S9-S12).
## Derived from the worry-group sample (ddR), matching the published tables.
no_ones <- ddR %>%
  group_by(ptN) %>%
  filter(all(uncertainty_before != 1)) %>%
  ungroup()

cat("\n@DATA  analysis subsets\n")
cat("  all trials            :", nrow(cleanedDataANOVA), "\n")
cat("  worry groups only     :", nrow(ddR), "trials from",
    length(unique(ddR$ptN)), "participants\n")
cat("  high worriers         :", nrow(cleanedDataANOVA_highworry), "trials\n")
cat("  low worriers          :", nrow(cleanedDataANOVA_lowworry), "trials\n")
cat("  no bottom-of-scale    :", nrow(no_ones), "trials\n")

## ---- 1.8 Search-level data for the stopping models -----------------------
## The stopping models (Tables S17-S19) treat each sample within a trial as
## an observation: at every sample the participant either continued or
## stopped and declared. `combinedImagesList` is the ordered list of images
## seen in that trial, so one row per image gives one row per decision point.
##
##   search   how many samples have been made so far this trial
##   stopped  1 at the final sample of the trial, 0 otherwise
##   newInfo  1 if this sample returned a previously unseen image
##
## Trials with no samples contribute no decision points and are therefore
## not represented; they are retained in the sampling-count models above,
## where zero is a valid count.

df_long <- do.call(rbind, lapply(which(ddR$openings != 0), function(i) {
  images <- strsplit(ddR$combinedImagesList[i], ", ")[[1]]
  if (length(images) == 0 || (length(images) == 1 && images[1] == "")) return(NULL)
  data.frame(
    search       = seq_along(images),
    stopped      = as.integer(seq_along(images) == length(images)),
    ptN          = ddR$ptN[i],
    newInfo      = as.integer(grepl("newsampling", images)),
    infoRatio    = as.character(ddR$infoRatio[i]),
    worryLevel   = as.character(ddR$worryLevel[i]),
    Age          = ddR$Age[i],
    Sex          = ddR$Sex[i],
    SES          = ddR$SES[i],
    worry_before = ddR$worry_before[i],
    stringsAsFactors = FALSE
  )
}))

df_long$stopped      <- as.numeric(df_long$stopped)
df_long$newInfo      <- as.numeric(df_long$newInfo)
df_long$worry_before <- as.numeric(df_long$worry_before)
df_long$ptN          <- factor(df_long$ptN)

## Reference levels for the search-level models.
##
## NOTE ON infoRatio: the stopping models take "1 to 40" (information poor)
## as the reference, so their information-environment coefficient is named
## `infoRatio1 to 5` and reads as the information-RICH environment relative
## to the poor one. This is the opposite reference to the trial-level
## sampling models above, where "1 to 5" is the reference and the
## coefficient reads as information scarcity. Both are stated here
## explicitly so the direction of each table is unambiguous.
df_long$infoRatio  <- relevel(factor(df_long$infoRatio),  ref = "1 to 40")
df_long$worryLevel <- relevel(factor(df_long$worryLevel), ref = "lowWorry")

df_long_high <- df_long[df_long$worryLevel == "highWorry", ]
df_long_low  <- df_long[df_long$worryLevel == "lowWorry", ]

cat("\n@DATA  search-level frame\n")
cat("  decision points       :", nrow(df_long), "from",
    length(unique(df_long$ptN)), "participants\n")
cat("  infoRatio reference   :", levels(df_long$infoRatio)[1], "\n")
cat("  worryLevel reference  :", levels(df_long$worryLevel)[1], "\n")


## #########################################################################
## @HELPERS — SHARED FUNCTIONS
## #########################################################################

## ---- Negative binomial mixed model shorthand ----------------------------
## Every sampling-count model in the Supplement has the same shape: a
## negative binomial GLMM with a participant random intercept. This wrapper
## keeps the calls readable and guarantees they are specified identically.

NB <- function(formula, data) {
  mixed_model(fixed  = formula,
              random = ~ 1 | ptN,
              data   = data,
              family = GLMMadaptive::negative.binomial(),
              iter_EM = 30)
}

## ---- Confidence intervals from a GLMMadaptive fit -----------------------
## Wald 95% CI: estimate +/- 1.96 * SE on the link (log) scale, then
## exponentiated. For a negative binomial model the exponentiated quantity
## is an incidence rate ratio (IRR); for a logistic model it is an odds
## ratio (OR). Set exponentiate = FALSE to stay on the link scale.
##
## Returns one row per coefficient, plus a `formatted` column giving the
## estimate in the style used in the manuscript text.

ci_table <- function(model, exponentiate = TRUE, ratio_label = "OR", table_ref = "") {
  co  <- summary(model)$coef_table
  est <- co[, "Estimate"]
  se  <- co[, "Std.Err"]
  z   <- co[, "z-value"]
  p   <- co[, "p-value"]

  lo <- est - 1.96 * se
  hi <- est + 1.96 * se

  ratio    <- if (exponentiate) exp(est) else est
  ratio_lo <- if (exponentiate) exp(lo)  else lo
  ratio_hi <- if (exponentiate) exp(hi)  else hi

  out <- data.frame(
    term   = rownames(co),
    ratio  = round(ratio, 2),
    ci_low = round(ratio_lo, 2),
    ci_high = round(ratio_hi, 2),
    beta   = round(est, 2),
    z      = round(z, 1),
    p      = p,
    row.names = NULL
  )
  out$formatted <- sprintf(
    "(%s %.2f, 95%% CI (%.2f, %.2f), P %s, beta = %.2f, z = %.1f%s)",
    ratio_label, out$ratio, out$ci_low, out$ci_high,
    ifelse(out$p < 0.001, "< 0.001",
           paste0("= ", formatC(out$p, format = "f", digits = 3))),
    out$beta, out$z,
    ifelse(nzchar(table_ref), paste0(", ", table_ref), "")
  )
  out
}

## ---- Confidence intervals from an nlme::lme fit -------------------------
## lme objects use intervals(), not confint(). Models of self-reported worry
## and error magnitude are Gaussian, so these stay on the response scale.

lme_ci <- function(model) {
  tt <- summary(model)$tTable
  iv <- nlme::intervals(model, which = "fixed")$fixed
  data.frame(
    term    = rownames(tt),
    beta    = round(tt[, "Value"], 3),
    ci_low  = round(iv[, "lower"], 3),
    ci_high = round(iv[, "upper"], 3),
    se      = round(tt[, "Std.Error"], 3),
    t       = round(tt[, "t-value"], 2),
    p       = signif(tt[, "p-value"], 3),
    row.names = NULL
  )
}

## ---- p-value formatting for tables --------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001",
                            formatC(p, format = "f", digits = 3))

## ---- Mean / median / SD block used by the descriptive tables ------------
desc_block <- function(x) {
  c(Mean   = round(mean(x, na.rm = TRUE), 2),
    Median = round(median(x, na.rm = TRUE), 2),
    SD     = round(sd(x, na.rm = TRUE), 2))
}


## #########################################################################
## @S2.1 — TABLE S2.1. DESCRIPTIVE STATISTICS
## #########################################################################
## Trial-level mean, median and SD of samples per trial. Distributions are
## right skewed with substantial individual variation, so the published
## table carries a note that these are for description only.
## -------------------------------------------------------------------------

cat("\n\n@S2.1  Table S2.1 — Descriptive statistics\n")

## a) Samples by trait worry
cat("\na) Number of samples by trait worry\n")
s21a <- t(sapply(c(`High worry` = "highWorry", `Low worry` = "lowWorry"),
                 function(g) desc_block(ddR$openings[ddR$worryLevel == g])))
print(knitr::kable(s21a, format = "pipe"))

## b) Samples by trait worry x informational environment
cat("\nb) Number of samples by trait worry x informational environment\n")
s21b <- do.call(rbind, lapply(
  list(c("highWorry", "1 to 40"), c("highWorry", "1 to 5"),
       c("lowWorry",  "1 to 40"), c("lowWorry",  "1 to 5")),
  function(k) {
    x <- ddR$openings[ddR$worryLevel == k[1] & ddR$infoRatio == k[2]]
    data.frame(`Trait Worry` = ifelse(k[1] == "highWorry", "High", "Low"),
               Environment   = ifelse(k[2] == "1 to 40", "1:40 (poor)", "1:5 (rich)"),
               t(desc_block(x)), check.names = FALSE)
  }))
print(knitr::kable(s21b, row.names = FALSE, format = "pipe"))

## c) Samples by stimuli aversiveness
## The aversive condition is `angry`; the appetitive condition is
## `happyWithAngry`.
cat("\nc) Number of samples by stimuli aversiveness\n")
s21c <- t(sapply(c(High = "angry", Low = "happyWithAngry"),
                 function(g) desc_block(ddR$openings[ddR$emotion == g])))
print(knitr::kable(s21c, format = "pipe"))

## d) Samples by stimuli uncertainty
## Stimuli uncertainty is the variance of the initial image distribution, so
## high uncertainty is the LOW precision condition.
cat("\nd) Number of samples by stimuli uncertainty\n")
s21d <- t(sapply(c(High = "lowPrecision", Low = "highPrecision"),
                 function(g) desc_block(ddR$openings[ddR$precision == g])))
print(knitr::kable(s21d, format = "pipe"))


## #########################################################################
## @S2.2 — TABLE S2.2. STATISTICAL TESTS
## #########################################################################
## Rank-based group comparisons on sampling, the error-magnitude
## equivalence test, and the Bayes factor for information scarcity.
##
## All rank tests are run on PARTICIPANT MEANS, not individual trials, so
## each participant contributes one observation and the independence
## assumption holds.
## -------------------------------------------------------------------------

cat("\n\n@S2.2  Table S2.2 — Statistical tests\n")

## ---- Participant-level means used by the rank tests ---------------------
pt_means <- aggregate(openings ~ ptID + worryLevel, data = ddR, FUN = mean)

pt_means_env <- aggregate(openings ~ ptID + worryLevel + infoRatio,
                          data = ddR, FUN = mean)

pt_rich <- pt_means_env[pt_means_env$infoRatio == "1 to 5",  ]
pt_poor <- pt_means_env[pt_means_env$infoRatio == "1 to 40", ]

## Wide format for the paired tests: one row per participant, one column
## per information environment.
pt_wide <- reshape(pt_means_env[, c("ptID", "worryLevel", "infoRatio", "openings")],
                   idvar   = c("ptID", "worryLevel"),
                   timevar = "infoRatio",
                   direction = "wide")
col_rich <- "openings.1 to 5"
col_poor <- "openings.1 to 40"

pt_wide_high <- pt_wide[pt_wide$worryLevel == "highWorry", ]
pt_wide_low  <- pt_wide[pt_wide$worryLevel == "lowWorry", ]

## ---- S2.2a. Mann-Whitney U, high vs low worriers ------------------------
## Between-participant comparison, so unpaired. Effect size is the rank
## biserial correlation with its 95% CI.

mw_overall <- wilcox.test(openings ~ worryLevel, data = pt_means)
rb_overall <- rank_biserial(openings ~ worryLevel, data = pt_means)

mw_rich <- wilcox.test(openings ~ worryLevel, data = pt_rich)
rb_rich <- rank_biserial(openings ~ worryLevel, data = pt_rich)

mw_poor <- wilcox.test(openings ~ worryLevel, data = pt_poor)
rb_poor <- rank_biserial(openings ~ worryLevel, data = pt_poor)

s22a <- data.frame(
  Environment = c("Overall", "Information rich", "Information poor"),
  U   = round(c(mw_overall$statistic, mw_rich$statistic, mw_poor$statistic), 0),
  P   = sapply(c(mw_overall$p.value, mw_rich$p.value, mw_poor$p.value), fmt_p),
  rrb = round(c(rb_overall$r_rank_biserial, rb_rich$r_rank_biserial,
                rb_poor$r_rank_biserial), 2),
  `95% CI` = sprintf("(%.2f, %.2f)",
                     c(rb_overall$CI_low, rb_rich$CI_low, rb_poor$CI_low),
                     c(rb_overall$CI_high, rb_rich$CI_high, rb_poor$CI_high)),
  check.names = FALSE, stringsAsFactors = FALSE
)
cat("\na) Mann-Whitney U tests comparing samples per trial in high vs low worriers\n")
print(knitr::kable(s22a, row.names = FALSE, format = "pipe"))

## ---- S2.2b. Wilcoxon signed-rank, information poor vs rich --------------
## Within-participant comparison, so paired: every participant completed
## trials in both environments.

ws_all  <- wilcox.test(pt_wide[[col_poor]],      pt_wide[[col_rich]],      paired = TRUE)
ws_high <- wilcox.test(pt_wide_high[[col_poor]], pt_wide_high[[col_rich]], paired = TRUE)
ws_low  <- wilcox.test(pt_wide_low[[col_poor]],  pt_wide_low[[col_rich]],  paired = TRUE)

rb_all_ws  <- rank_biserial(pt_wide[[col_poor]],      pt_wide[[col_rich]],      paired = TRUE)
rb_high_ws <- rank_biserial(pt_wide_high[[col_poor]], pt_wide_high[[col_rich]], paired = TRUE)
rb_low_ws  <- rank_biserial(pt_wide_low[[col_poor]],  pt_wide_low[[col_rich]],  paired = TRUE)

s22b <- data.frame(
  Group = c("Overall", "High worriers", "Low worriers"),
  W   = round(c(ws_all$statistic, ws_high$statistic, ws_low$statistic), 0),
  P   = sapply(c(ws_all$p.value, ws_high$p.value, ws_low$p.value), fmt_p),
  rrb = round(c(rb_all_ws$r_rank_biserial, rb_high_ws$r_rank_biserial,
                rb_low_ws$r_rank_biserial), 2),
  `95% CI` = sprintf("(%.2f, %.2f)",
                     c(rb_all_ws$CI_low, rb_high_ws$CI_low, rb_low_ws$CI_low),
                     c(rb_all_ws$CI_high, rb_high_ws$CI_high, rb_low_ws$CI_high)),
  check.names = FALSE, stringsAsFactors = FALSE
)
cat("\nb) Wilcoxon signed-rank tests, information poor vs information rich\n")
print(knitr::kable(s22b, row.names = FALSE, format = "pipe"))

## Raw difference in samples per trial, quoted in the main text.
cat("\n   Raw difference in samples per trial (poor minus rich):\n")
for (g in c("highWorry", "lowWorry")) {
  w <- pt_wide[pt_wide$worryLevel == g, ]
  cat(sprintf("     %-10s poor %.2f, rich %.2f, difference %+.2f\n",
              g, mean(w[[col_poor]]), mean(w[[col_rich]]),
              mean(w[[col_poor]]) - mean(w[[col_rich]])))
}
cat(sprintf("     %-10s poor %.2f, rich %.2f, difference %+.2f\n", "overall",
            mean(pt_wide[[col_poor]]), mean(pt_wide[[col_rich]]),
            mean(pt_wide[[col_poor]]) - mean(pt_wide[[col_rich]])))

## ---- S2.2c. Error magnitude: high vs low worriers -----------------------
## Error magnitude is a continuous rating difference, not a count, so the
## frequentist test is a Gaussian linear mixed model with a participant
## random intercept. The equivalence claim is then made with a region of
## practical equivalence of +/- 0.1 trial-level SD, reporting the proportion
## of the 89% credible interval that falls inside that window.

ACC <- "accuracy_fullset"   ## measure reported in the manuscript

## Group means and SDs, at participant level to match Figure 2d.
acc_pt <- aggregate(as.formula(paste(ACC, "~ ptID + worryLevel")),
                    data = ddR, FUN = mean)
hi <- acc_pt$worryLevel == "highWorry"
lo <- acc_pt$worryLevel == "lowWorry"
mHi <- mean(acc_pt[[ACC]][hi]); sHi <- sd(acc_pt[[ACC]][hi])
mLo <- mean(acc_pt[[ACC]][lo]); sLo <- sd(acc_pt[[ACC]][lo])

## For this table only, trait worry is relevelled so that highWorry is the
## reference. The single coefficient then reads as low minus high, which is
## the direction reported in the published table.
ddR_acc <- ddR
ddR_acc$worryLevel <- relevel(ddR_acc$worryLevel, ref = "highWorry")

## Frequentist group difference, trial level.
f_acc    <- as.formula(paste(ACC, "~ worryLevel"))
acc_lme  <- lme(f_acc, random = ~ 1 | ptN, data = ddR_acc, na.action = "na.omit")
acc_co   <- summary(acc_lme)$tTable["worryLevellowWorry", ]
acc_beta <- acc_co["Value"]; acc_se <- acc_co["Std.Error"]; acc_p <- acc_co["p-value"]

## ROPE, scaled to the trial-level SD of the measure.
acc_sd  <- sd(ddR[[ACC]], na.rm = TRUE)
acc_win <- 0.1 * acc_sd
acc_bayes <- stan_glm(f_acc, data = ddR_acc, refresh = 0)
acc_rope  <- rope(acc_bayes, range = c(-acc_win, acc_win), ci = 0.89)
acc_pct   <- round(acc_rope$ROPE_Percentage[
                     acc_rope$Parameter == "worryLevellowWorry"] * 100, 1)

s22c <- data.frame(
  Statistic = c("Mean error, high worriers (SD)",
                "Mean error, low worriers (SD)",
                "Group difference (beta, low - high)",
                "Standard error",
                "p (linear mixed model)",
                "Standardised difference (SD units)",
                "ROPE (+/- 0.1 x SD)",
                "% of 89% CI inside ROPE"),
  Value = c(sprintf("%.2f (%.2f)", mHi, sHi),
            sprintf("%.2f (%.2f)", mLo, sLo),
            sprintf("%+.3f", acc_beta),
            sprintf("%.3f",  acc_se),
            sprintf("%.3f",  acc_p),
            sprintf("%+.3f", acc_beta / acc_sd),
            sprintf("+/-%.3f", acc_win),
            sprintf("%.1f%%", acc_pct)),
  stringsAsFactors = FALSE
)
cat("\nc) Comparison of error magnitude between high and low worriers\n")
print(knitr::kable(s22c, row.names = FALSE, format = "pipe"))

## ---- S2.2d. Bayes factor: information scarcity on sampling --------------
## Paired comparison of participant mean samples in the information-poor
## versus information-rich environment, split by worry group. JZS prior with
## the default scale r = 0.707. BF10 > 1 favours a difference; BF01 > 1
## favours no difference.

bf_high <- ttestBF(x = pt_wide_high[[col_poor]], y = pt_wide_high[[col_rich]],
                   paired = TRUE)
bf_low  <- ttestBF(x = pt_wide_low[[col_poor]],  y = pt_wide_low[[col_rich]],
                   paired = TRUE)
bf_all  <- ttestBF(x = pt_wide[[col_poor]],      y = pt_wide[[col_rich]],
                   paired = TRUE)

bf10 <- c(high = exp(bf_high@bayesFactor$bf),
          low  = exp(bf_low@bayesFactor$bf),
          all  = exp(bf_all@bayesFactor$bf))

s22d <- data.frame(
  Comparison = c("Info-poor vs Info-rich (high worriers)",
                 "Info-poor vs Info-rich (low worriers)"),
  BF10 = round(bf10[c("high", "low")], 2),
  BF01 = round(1 / bf10[c("high", "low")], 2),
  Interpretation = ifelse(bf10[c("high", "low")] > 10, "Strong evidence for effect",
                   ifelse(bf10[c("high", "low")] > 3,  "Moderate evidence for effect",
                   ifelse(bf10[c("high", "low")] > 1,  "Anecdotal / inconclusive",
                                                       "Evidence for no effect"))),
  row.names = NULL, stringsAsFactors = FALSE
)
cat("\nd) Bayes factor: effect of information scarcity on sampling\n")
cat("   Paired comparison of mean samples per trial, JZS prior r = 0.707\n")
print(knitr::kable(s22d, row.names = FALSE, format = "pipe"))
cat(sprintf("\n   All participants pooled: BF10 = %.2f\n", bf10[["all"]]))


## #########################################################################
## @S3 — TABLE S3. EFFECT OF TRAIT WORRY ON SAMPLING
## #########################################################################
## Negative binomial GLMM. worryLevel reference is lowWorry, so the trait
## worry coefficient is the high-worry group relative to low.
## -------------------------------------------------------------------------

cat("\n\n@S3  Table S3 — Trait worry on sampling\n")

model_S3 <- NB(openings ~ worryLevel + Age + Sex + SES, ddR)
print(summary(model_S3))
print(ci_table(model_S3, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S3"))

## Estimated marginal means give the group contrast on the response scale,
## quoted in the main text as an IRR with 95% CI.
## NOTE: the data must be passed to emmeans explicitly here. The NB() helper
## fits the model inside a function, so the fitted object only remembers its
## data as an argument called "data"; emmeans cannot find that object by
## itself and errors with "unable to reconstruct the data". Handing it
## data = ddR resolves it.
cat("\n   Estimated marginal means by worry group:\n")
print(emmeans(model_S3, ~ worryLevel, data = ddR, type = "response"))


## #########################################################################
## @S4 — TABLE S4. EFFECT OF STATE WORRY ON SAMPLING
## #########################################################################
## Within-participant: state worry is the pre-sampling rating on that trial,
## so this is a trial-by-trial association rather than a group difference.
## Rated 1-5 and centred at 1, so the intercept sits at the lowest value
## that occurs and each unit is one point on the scale.
## -------------------------------------------------------------------------

cat("\n\n@S4  Table S4 — State worry on sampling\n")

model_S4 <- NB(openings ~ worry_before_c1 + Age + Sex + SES, ddR)
print(summary(model_S4))

ci_S4 <- ci_table(model_S4, exponentiate = TRUE, ratio_label = "IRR",
                  table_ref = "Table S4")
print(ci_S4)

## Effect across the full range of the scale (4 units from 1 to 5), quoted
## in the main text alongside the per-unit IRR.
b_S4  <- summary(model_S4)$coef_table["worry_before_c1", "Estimate"]
se_S4 <- summary(model_S4)$coef_table["worry_before_c1", "Std.Err"]
cat(sprintf(
  "\n   Full scale range (4 units): IRR %.2f, 95%% CI (%.2f, %.2f)\n",
  exp(4 * b_S4), exp(4 * (b_S4 - 1.96 * se_S4)), exp(4 * (b_S4 + 1.96 * se_S4))))


## #########################################################################
## @S5 — TABLE S5. STIMULI AVERSIVENESS ON SAMPLING
## @S6 — TABLE S6. STIMULI AVERSIVENESS ON SELF-REPORTED WORRY
## #########################################################################
## Aversiveness is the experimental valence condition. The appetitive
## condition (`happyWithAngry`) is the reference, so the reported
## coefficient `emotionangry` is the aversive condition relative to it and
## increases with aversiveness.
## -------------------------------------------------------------------------

cat("\n\n@S5  Table S5 — Stimuli aversiveness on sampling\n")

model_S5 <- NB(openings ~ emotion + Age + Sex + SES, ddR)
print(summary(model_S5))
print(ci_table(model_S5, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S5"))

cat("\n\n@S6  Table S6 — Stimuli aversiveness on self-reported worry\n")
## Self-reported worry is a 1-5 rating, so Gaussian rather than count.

model_S6 <- lme(worry_before ~ emotion + Age + Sex + SES,
                random = ~ 1 | ptN, data = ddR, na.action = "na.omit")
print(summary(model_S6))
print(lme_ci(model_S6))


## #########################################################################
## @S7 — TABLE S7. STIMULI UNCERTAINTY ON SAMPLING
## @S8 — TABLE S8. STIMULI UNCERTAINTY ON WORRY
## #########################################################################
## Stimuli uncertainty is the variance of the initial image distribution,
## held in `precision` (inverse variance). The high-precision condition is
## the reference, so the reported coefficient `precisionlowPrecision` is the
## high-variance condition relative to the low-variance one and increases
## with stimulus uncertainty.
##
## The `precision` factor is the stimulus-uncertainty measure throughout the
## paper, in both these single-predictor tables and the combined models
## (Tables S11/S12), so the same operationalisation is used everywhere.
## -------------------------------------------------------------------------

cat("\n\n@S7  Table S7 — Stimuli uncertainty on sampling\n")

model_S7 <- NB(openings ~ precision + Age + Sex + SES, ddR)
print(summary(model_S7))
print(ci_table(model_S7, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S7"))

cat("\n\n@S8  Table S8 — Stimuli uncertainty on worry\n")

model_S8 <- lme(worry_before ~ precision + Age + Sex + SES,
                random = ~ 1 | ptN, data = ddR, na.action = "na.omit")
print(summary(model_S8))
print(lme_ci(model_S8))


## #########################################################################
## @S9  — TABLE S9.  SELF-REPORTED UNCERTAINTY ON SAMPLING
## @S10 — TABLE S10. SELF-REPORTED UNCERTAINTY ON WORRY
## #########################################################################
## Subjective uncertainty, rated before sampling. Fitted on `no_ones`:
## participants who never used the bottom of the scale, so the predictor has
## within-person variance for everyone in the model.
## -------------------------------------------------------------------------

cat("\n\n@S9  Table S9 — Self-reported uncertainty on sampling\n")

model_S9 <- NB(openings ~ uncertainty_before + Age + Sex + SES, no_ones)
print(summary(model_S9))
print(ci_table(model_S9, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S9"))

cat("\n\n@S10  Table S10 — Self-reported uncertainty on worry\n")

model_S10 <- lme(worry_before ~ uncertainty_before + Age + Sex + SES,
                 random = ~ 1 | ptN, data = no_ones, na.action = "na.omit")
print(summary(model_S10))
print(lme_ci(model_S10))


## #########################################################################
## @S11 — TABLE S11. BOTH UNCERTAINTY TYPES ON SAMPLING
## @S12 — TABLE S12. BOTH UNCERTAINTY TYPES ON WORRY
## #########################################################################
## Objective and subjective uncertainty entered together, to test whether
## the stimulus effect survives adjustment for the self-report.
##
## Stimulus uncertainty is the `precision` factor, the same measure used in
## Tables S7/S8. Note that the stimulus-uncertainty coefficient is positive
## when entered on its own (Table S7) but negative here, once self-reported
## uncertainty is also in the model; both are non-significant, so this is a
## redistribution of shared variance between two correlated predictors, not
## a reversal of any reliable effect.
## -------------------------------------------------------------------------

cat("\n\n@S11  Table S11 — Both uncertainty types on sampling\n")

model_S11 <- NB(openings ~ uncertainty_before + precision + Age + Sex + SES, no_ones)
print(summary(model_S11))
print(ci_table(model_S11, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S11"))

cat("\n\n@S12  Table S12 — Both uncertainty types on worry\n")

model_S12 <- lme(worry_before ~ uncertainty_before + precision + Age + Sex + SES,
                 random = ~ 1 | ptN, data = no_ones, na.action = "na.omit")
print(summary(model_S12))
print(lme_ci(model_S12))


## #########################################################################
## @S13 — TABLE S13. INFORMATION ENVIRONMENT, OVERALL GROUP
## @S14 — TABLE S14. INFORMATION ENVIRONMENT, HIGH WORRIERS
## @S15 — TABLE S15. INFORMATION ENVIRONMENT, LOW WORRIERS
## @S16 — TABLE S16. INFORMATION ENVIRONMENT x TRAIT WORRY
## #########################################################################
## infoRatio reference is "1 to 5" (information rich), so the coefficient
## reads as information scarcity relative to information wealth.
## -------------------------------------------------------------------------

cat("\n\n@S13  Table S13 — Information environment, overall group\n")

## "Overall" worry-group sample (high + low), matching the published table.
model_S13 <- NB(openings ~ infoRatio + Age + Sex + SES, ddR)
print(summary(model_S13))
print(ci_table(model_S13, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S13"))

cat("\n\n@S14  Table S14 — Information environment, high worriers\n")

model_S14 <- NB(openings ~ infoRatio + Age + Sex + SES, cleanedDataANOVA_highworry)
print(summary(model_S14))
print(ci_table(model_S14, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S14"))

cat("\n\n@S15  Table S15 — Information environment, low worriers\n")

model_S15 <- NB(openings ~ infoRatio + Age + Sex + SES, cleanedDataANOVA_lowworry)
print(summary(model_S15))
print(ci_table(model_S15, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S15"))

cat("\n\n@S16  Table S16 — Information environment x trait worry\n")
## The interaction term tests whether the information-scarcity effect
## differs between worry groups.

model_S16 <- NB(openings ~ infoRatio * worryLevel + Age + Sex + SES, ddR)
print(summary(model_S16))
print(ci_table(model_S16, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S16"))


## #########################################################################
## @S17 — TABLE S17. STOPPING SEARCH AND DECLARING
## #########################################################################
## Logistic GLMM on the search-level frame. Outcome is whether the
## participant stopped at this decision point (1) or continued (0), so a
## positive coefficient means a higher chance of stopping. Exponentiated
## coefficients are odds ratios.
## -------------------------------------------------------------------------

cat("\n\n@S17  Table S17 — Stopping search and declaring\n")

model_S17 <- mixed_model(
  fixed  = stopped ~ newInfo + infoRatio + search + Age + Sex + SES,
  random = ~ 1 | ptN,
  data   = df_long,
  family = binomial()
)
print(summary(model_S17))
print(ci_table(model_S17, exponentiate = TRUE, ratio_label = "OR",
               table_ref = "Table S17"))


## #########################################################################
## @S18 — TABLE S18. STOPPING, INTERACTION MODEL
## #########################################################################
## Adds a three-way interaction between informational novelty, state worry
## and trait worry. State worry is centred at 1 so that the main effects are
## read at the lowest value on the scale rather than at an off-scale zero;
## this leaves every interaction term unchanged.
## -------------------------------------------------------------------------

cat("\n\n@S18  Table S18 — Stopping, interaction model\n")

df_long$worry_before_c1 <- df_long$worry_before - 1

model_S18 <- mixed_model(
  fixed  = stopped ~ newInfo * worry_before_c1 * worryLevel +
                     infoRatio + search + Age + Sex + SES,
  random = ~ 1 | ptN,
  data   = df_long,
  family = binomial()
)
print(summary(model_S18))
print(ci_table(model_S18, exponentiate = TRUE, ratio_label = "OR",
               table_ref = "Table S18"))


## #########################################################################
## @S19 — TABLE S19. SUBGROUP ANALYSES
## #########################################################################
## Two sets of subgroup models, each fitted separately in high and low
## worriers:
##   (i)  the effect of informational novelty on stopping
##   (ii) the effect of self-reported uncertainty on sampling
## -------------------------------------------------------------------------

cat("\n\n@S19  Table S19 — Subgroup analyses\n")

## ---- (i) Novelty on stopping, by worry group ----------------------------
cat("\n  (i) Effect of novelty on stopping\n")

model_S19_stop_high <- mixed_model(
  fixed  = stopped ~ newInfo + worry_before + infoRatio + search + Age + Sex + SES,
  random = ~ 1 | ptN,
  data   = df_long_high,
  family = binomial(),
  iter_EM = 30
)
cat("\n  High worriers:\n")
print(summary(model_S19_stop_high))
print(ci_table(model_S19_stop_high, exponentiate = TRUE, ratio_label = "OR",
               table_ref = "Table S19"))

model_S19_stop_low <- mixed_model(
  fixed  = stopped ~ newInfo + worry_before + infoRatio + search + Age + Sex + SES,
  random = ~ 1 | ptN,
  data   = df_long_low,
  family = binomial(),
  iter_EM = 30
)
cat("\n  Low worriers:\n")
print(summary(model_S19_stop_low))
print(ci_table(model_S19_stop_low, exponentiate = TRUE, ratio_label = "OR",
               table_ref = "Table S19"))

## ---- (ii) Self-reported uncertainty on sampling, by worry group ---------
cat("\n  (ii) Effect of self-reported uncertainty on sampling\n")

no_ones_high <- no_ones[no_ones$worryLevel == "highWorry", ]
no_ones_low  <- no_ones[no_ones$worryLevel == "lowWorry", ]

model_S19_unc_high <- NB(openings ~ uncertainty_before + Age + Sex + SES, no_ones_high)
cat("\n  High worriers:\n")
print(summary(model_S19_unc_high))
print(ci_table(model_S19_unc_high, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S19"))

model_S19_unc_low <- NB(openings ~ uncertainty_before + Age + Sex + SES, no_ones_low)
cat("\n  Low worriers:\n")
print(summary(model_S19_unc_low))
print(ci_table(model_S19_unc_low, exponentiate = TRUE, ratio_label = "IRR",
               table_ref = "Table S19"))


## #########################################################################
## @ROPE — REGION OF PRACTICAL EQUIVALENCE (MAIN TEXT)
## #########################################################################
## The ROPE convention used throughout the paper is +/- 0.1 trial-level SD
## of the outcome, and the quantity reported is the proportion of the 89%
## credible interval that falls inside that window. The 89% interval, rather
## than the full posterior, is the quantity reported in every ROPE result.
##
## The information-poverty ROPE is fitted separately in each worry group on
## raw sample counts, with no covariates: with covariates in the model the
## intercept is no longer the group mean and the window would no longer be
## on the same scale as the reported effect.
##
## The window is symmetric about zero, so the proportion inside it does not
## depend on which information environment is the reference level.
## -------------------------------------------------------------------------

cat("\n\n@ROPE  Region of practical equivalence\n")

## Scaled to the trial-level SD of samples across both worry groups, so the
## same window applies to each group's model.
rope_window <- 0.1 * sd(ddR$openings, na.rm = TRUE)
cat(sprintf("  ROPE window for sampling: +/- %.3f samples (0.1 SD)\n", rope_window))

## Pull the information-environment coefficient out of a posterior draw
## matrix by name, so the code does not depend on column position.
env_draws <- function(fit) {
  d <- as.matrix(fit)
  j <- grep("^infoRatio", colnames(d))
  if (length(j) != 1)
    stop("Expected exactly one infoRatio coefficient in the posterior.")
  d[, j]
}

## ---- Information poverty, low worriers ----------------------------------
modelEnvLowWorry <- stan_glm(openings ~ infoRatio,
                             data = cleanedDataANOVA_lowworry, refresh = 0)
cat("\n  Low worriers:\n")
print(rope(env_draws(modelEnvLowWorry),
           range = c(-rope_window, rope_window),
           ci = 0.89))

## ---- Information poverty, high worriers ---------------------------------
modelEnvHighWorry <- stan_glm(openings ~ infoRatio,
                              data = cleanedDataANOVA_highworry, refresh = 0)
cat("\n  High worriers:\n")
print(rope(env_draws(modelEnvHighWorry),
           range = c(-rope_window, rope_window),
           ci = 0.89))

## ---- Error magnitude, high vs low worriers ------------------------------
## Same convention, scaled to the trial-level SD of error magnitude. This is
## the result tabulated in Table S2.2c and quoted in the main text.
## Refits the model if @S2.2 has not been run in this session.
if (!exists("acc_bayes")) {
  ACC     <- "accuracy_fullset"
  ddR_acc <- ddR
  ddR_acc$worryLevel <- relevel(ddR_acc$worryLevel, ref = "highWorry")
  acc_win   <- 0.1 * sd(ddR[[ACC]], na.rm = TRUE)
  acc_bayes <- stan_glm(as.formula(paste(ACC, "~ worryLevel")),
                        data = ddR_acc, refresh = 0)
}
cat("\n  Error magnitude (high vs low worriers):\n")
print(rope(acc_bayes, range = c(-acc_win, acc_win), ci = 0.89))

## ---- Convergence check --------------------------------------------------
cat("\n  Maximum Rhat across ROPE models:\n")
for (nm in c("modelEnvLowWorry", "modelEnvHighWorry", "acc_bayes")) {
  cat(sprintf("    %-20s %.3f\n", nm, max(summary(get(nm))[, "Rhat"], na.rm = TRUE)))
}


## #########################################################################
## @FIG2 — FIGURE 2. PANELS a-f
## #########################################################################
## Six panels, assembled with patchwork. Every panel is built from
## PARTICIPANT MEANS so that each density reflects one point per person.
##
##   a  samples by trait worry
##   b  samples by information environment and trait worry
##   c  trait worry x information scarcity interaction
##   d  error magnitude by trait worry
##   e  samples by stimuli aversiveness
##   f  samples by stimuli uncertainty
##
## Shared palette: teal for low worry, gold for high worry, purple pair for
## the stimulus manipulations.
## -------------------------------------------------------------------------

cat("\n\n@FIG2  Building Figure 2 panels\n")

col_low  <- "#00AFBB"   ## low worry
col_high <- "#E7B800"   ## high worry
col_stim <- c("#8a1a88", "#d69fd5")   ## stimulus manipulations

## ---- Panel a. Samples by trait worry ------------------------------------
openings_by_worryLevel <- aggregate(openings ~ ptID + worryLevel,
                                    data = ddR, FUN = mean)

SamplesByTW_v <- ggplot(openings_by_worryLevel,
                        aes(x = worryLevel, y = openings, fill = worryLevel)) +
  geom_violin(alpha = 0.4, aes(color = worryLevel), linewidth = 0.8) +
  stat_summary(fun = "mean", geom = "crossbar", width = 1, linewidth = 0.25,
               aes(color = worryLevel)) +
  scale_fill_manual(values = c(lowWorry = col_low, highWorry = col_high)) +
  scale_color_manual(values = c(lowWorry = col_low, highWorry = col_high)) +
  scale_x_discrete(labels = c(lowWorry = "Low", highWorry = "High"),
                   limits = c("lowWorry", "highWorry")) +
  labs(y = "Mean Samples", x = "Trait Worry") +
  theme_classic() +
  theme(legend.position = "none")

## ---- Panel b. Samples by information environment and trait worry --------
openings_by_infoRatio <- aggregate(openings ~ ptID + infoRatio + worryLevel,
                                   data = ddR, FUN = mean)
openings_by_infoRatio$worryGroup <- factor(
  ifelse(openings_by_infoRatio$worryLevel == "lowWorry", "Low", "High"),
  levels = c("Low", "High"))
openings_by_infoRatio$infoRatio <- factor(openings_by_infoRatio$infoRatio,
                                          levels = c("1 to 5", "1 to 40"))

InfoPlotViolin <- ggplot(openings_by_infoRatio,
                         aes(x = infoRatio, y = openings, fill = worryGroup)) +
  geom_violin(alpha = 0.4, aes(color = worryGroup), linewidth = 0.8,
              position = position_dodge(width = 0.8)) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.4, linewidth = 0.25,
               aes(color = NA), position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(Low = col_low, High = col_high),
                    labels = c(Low = "Low Worry", High = "High Worry"),
                    name = NULL) +
  scale_color_manual(values = c(Low = col_low, High = col_high), guide = "none") +
  scale_x_discrete(labels = c(`1 to 5` = "Rich", `1 to 40` = "Poor"),
                   limits = c("1 to 5", "1 to 40")) +
  labs(y = "Mean Samples", x = "Information") +
  guides(fill = guide_legend(override.aes = list(color = NA, linewidth = 0),
                             reverse = TRUE)) +
  theme_classic() +
  theme(legend.title = element_text(size = 10),
        legend.position = c(0.12, 1),
        legend.justification = c(0, 1),
        legend.key.width  = unit(1, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.background = element_rect(fill = "white", color = NA))

## ---- Panel c. Trait worry x information scarcity interaction ------------
## Estimated marginal means from a Gaussian mixed model on raw counts, used
## here for display: the panel shows the pattern of means with standard
## errors, while the inferential test is the negative binomial interaction
## model in @S16.

mod_emm <- lmer(openings ~ worryGroup * infoWealth + (1 | ptN), data = ddR)
emm_df  <- as.data.frame(emmeans(mod_emm, ~ infoWealth * worryGroup))

## Position the line-end labels at the rightmost x value.
last_x     <- tail(levels(emm_df$infoWealth), 1)
low_end_y  <- emm_df$emmean[emm_df$worryGroup == "Low"  & emm_df$infoWealth == last_x]
high_end_y <- emm_df$emmean[emm_df$worryGroup == "High" & emm_df$infoWealth == last_x]

intPlot <- ggplot(emm_df, aes(x = infoWealth, y = emmean,
                              colour = worryGroup, group = worryGroup)) +
  geom_point(position = position_dodge(width = 0.15), size = 2.5) +
  geom_line(position = position_dodge(width = 0.15), linewidth = 0.8) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.1, position = position_dodge(width = 0.15)) +
  scale_colour_manual(values = c(Low = col_low, High = col_high), guide = "none") +
  labs(x = "Information", y = "Estimated Mean Samples") +
  annotate("text", x = length(levels(emm_df$infoWealth)) + 0.2, y = low_end_y,
           label = "Low Worry",  angle = 90, size = 3, color = col_low,  hjust = 0.5) +
  annotate("text", x = length(levels(emm_df$infoWealth)) + 0.2, y = high_end_y,
           label = "High Worry", angle = 90, size = 3, color = col_high, hjust = 0.5) +
  theme_classic()

## ---- Panel d. Error magnitude by trait worry ----------------------------
## Built from the same participant means as Table S2.2c.
accuracy_by_worryLevel <- acc_pt
names(accuracy_by_worryLevel)[names(accuracy_by_worryLevel) == ACC] <- "err"

ErrorMagPlot_v <- ggplot(accuracy_by_worryLevel,
                         aes(x = worryLevel, y = err, fill = worryLevel)) +
  geom_violin(alpha = 0.4, aes(color = worryLevel), linewidth = 0.8) +
  stat_summary(fun = "mean", geom = "crossbar", width = 1, linewidth = 0.25,
               aes(color = worryLevel)) +
  scale_fill_manual(values = c(lowWorry = col_low, highWorry = col_high)) +
  scale_color_manual(values = c(lowWorry = col_low, highWorry = col_high)) +
  scale_x_discrete(labels = c(lowWorry = "Low", highWorry = "High"),
                   limits = c("lowWorry", "highWorry")) +
  labs(y = "Overall Error Magnitude", x = "Trait Worry") +
  theme_classic() +
  theme(legend.position = "none")

## ---- Panel e. Samples by stimuli aversiveness ---------------------------
openings_by_emotion <- aggregate(openings ~ ptID + emotion,
                                 data = ddR, FUN = mean)

SamplesByEmo_v <- ggplot(openings_by_emotion,
                         aes(x = emotion, y = openings, fill = emotion)) +
  geom_violin(alpha = 0.4, aes(color = emotion), linewidth = 0.8) +
  stat_summary(fun = "mean", geom = "crossbar", width = 1, linewidth = 0.25,
               aes(color = emotion)) +
  scale_fill_manual(values = col_stim) +
  scale_color_manual(values = col_stim) +
  scale_x_discrete(labels = c(happyWithAngry = "Low", angry = "High"),
                   limits = c("happyWithAngry", "angry")) +
  labs(y = "Mean Samples", x = "Stimuli Aversiveness") +
  theme_classic() +
  theme(legend.position = "none")

## ---- Panel f. Samples by stimuli uncertainty ----------------------------
## High stimuli uncertainty is the low-precision condition.
openings_by_precision <- aggregate(openings ~ ptID + precision,
                                   data = ddR, FUN = mean)

SamplesByVar_v <- ggplot(openings_by_precision,
                         aes(x = precision, y = openings, fill = precision)) +
  geom_violin(alpha = 0.4, aes(color = precision), linewidth = 0.8) +
  stat_summary(fun = "mean", geom = "crossbar", width = 1, linewidth = 0.25,
               aes(color = precision)) +
  scale_fill_manual(values = col_stim) +
  scale_color_manual(values = col_stim) +
  scale_x_discrete(labels = c(lowPrecision = "High", highPrecision = "Low"),
                   limits = c("lowPrecision", "highPrecision")) +
  labs(y = "Mean Samples", x = "Stimuli Variance") +
  theme_classic() +
  theme(legend.position = "none")

## ---- Assemble Figure 2 --------------------------------------------------
## Three columns by two rows, with panel b given double width. Panel tags
## are nudged above each plot so they clear the axis titles.

figure2 <- (SamplesByTW_v | InfoPlotViolin | intPlot |
            ErrorMagPlot_v | (SamplesByEmo_v | SamplesByVar_v)) +
  plot_layout(ncol = 3, nrow = 2, widths = c(1, 2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag.position = c(0, 1.08),
        plot.tag = element_text(hjust = 0, vjust = 0),
        plot.margin = margin(t = 25, r = 5, b = 5, l = 5))

print(figure2)

## ggsave("figure2.png", figure2, width = 11, height = 7, dpi = 300, bg = "white")


## #########################################################################
## END OF ANALYSIS
## #########################################################################

cat("\n\n## Analysis complete.\n")
cat("## Section tags: @SETUP @DATA @HELPERS @S2.1 @S2.2 @S3 @S4 @S5 @S6 @S7\n")
cat("##               @S8 @S9 @S10 @S11 @S12 @S13 @S14 @S15 @S16 @S17 @S18\n")
cat("##               @S19 @ROPE @FIG2\n")

# Analysis code consolidated, organised and commented with Claude
# (Opus 4.8, Anthropic). All models and outputs were checked by the
# authors to reproduce the results reported in the manuscript.
