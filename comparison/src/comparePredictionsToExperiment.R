####################Joseph Nasser
# October 16 2019
#
# Evaluate predictive model against experimental data
#
# Requirements: 
#       - R 3.4.0
#         - data.table (1.11.4)
#         - GenomicRanges (1.28.6)
#         - ROCR (1.0-7)
#         - ggplot2 (3.0.0)
#         - caTools (1.17.1
#
# Usage:
#       Rscript comparePredictionsToExperiment.R --predictions pred.txt --experimentalData expt.txt --plotConfig plot.config.txt --predConfig pred.config.txt
#
# This code is intended to evaluate a predictive model against a set of experimental data
# The code will find regions in predictions which overlap the regions listed in experimentalData (See below)
# The output will be a file (expt.pred.txt, containing one row per experimental region) which merges the experimental data with the predictions.
# PR curves and scatter plots will also be generated. 
# 
# The plotConfig file should contain one line per PR curve plot. Scatter plots will be generated for all prediction columns defined in predConfig
#
# Configuring Prediction Columns:
# In the ideal use case each experimental element should correspond to exactly one element-gene prediction. 
# However, there are instances where a single experimental element overlaps multiple predictions (eg a large deletion) or when an experimentally tested element-gene pair is not present in the predictions files.
# The predConfig file describes how to handle these cases. 
#
# Other:
# - Assumes prediction columns are monotonic increasing! (A hack is employed for distance)
# - Can't distinguish between a missing prediction and a non-prediction


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(caTools))

option.list <- list(
  make_option("--predictions", type="character", help="Predictions file (accepts comma delimited list)"),
  make_option("--experimentalData", type="character", help="File listing perturbational data (accepts comma delimited list)"),
  make_option("--experimentalPositiveColumn", type="character", default="Regulated", help="Column of experimentalData to consider an experimental positive"),
  make_option("--plotConfig", type="character", help="File describing which plots to make"),
  make_option("--predConfig", type="character", help="File describing how to aggregate/fill prediction columns"),
  make_option("--ignoreExptMissingPredictions", default=FALSE, action="store_true", help="Ignore EG pairs which do not have predictions. Do not fill based on predConfig"),
  make_option("--outDir", default = ".", type="character", help="Output directory"),
  make_option("--codeDir", default = ".", type="character", help="code directory")
)
opt <- parse_args(OptionParser(option_list=option.list))
write.table(t(as.data.frame(opt)), file.path(opt$outDir, "params.txt"), sep = "\t", quote = F, col.names = F, row.names = T)

source(file.path(opt$codeDir, "comparison.R"))

#For Testing at Broad
# basedir <- "/seq/lincRNA/RAP/ABC/180419_for_ABCPaper/190530_NG_Revision/190531_ABC/190619_public_code/191008_test_new_codebase/"
# opt$predictions <- file.path(basedir, "Predictions/K562/EnhancerPredictionsAllPutative.txt.gz")
# opt$experimentalData <- file.path(basedir, "known/K562.KnownEnhancers.txt")
# opt$plotConfig <- file.path(basedir, "known/plot.config.txt")
# opt$predConfig <- file.path(basedir, "known/pred.config.txt")
# opt$outDir <- paste0(basedir, "/known/")
# opt$codeDir <- "/seq/lincRNA/jnasser/EP_prediction/code/LanderLab-EP-Prediction/src/libs/Predictions/"

#Read input data
print("Reading input files")
pred <- loadFileString(opt$predictions)
expt <- loadFileString(opt$experimentalData)
expt <- subset(expt, IncludeInModel)

#Test boolean column
#pred$rand.bool <- sample.int(n = 2, size = nrow(pred), replace = T) - 1

#QC Input Files
qcExpt(expt, opt)

# Merge experimental data and predictions
print("Merging experiment and predictions")
print(opt$ignoreExptMissingPredictions)
merged <- combineExptPred(expt = expt, 
                          pred = pred,
                          config = opt$predConfig,
                          fill.missing = !opt$ignoreExptMissingPredictions)
write.table(merged, file.path(opt$outDir, "expt.pred.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

# Make plots
print("Making plots")
merged <- prepForPlotting(merged)

if ("class" %in% colnames(merged))
  if ("genic" %in% merged$class) 
    merged <- subset(merged, class != "promoter")

makePlots(merged, opt$plotConfig, opt$experimentalPositiveColumn, opt$outDir)
