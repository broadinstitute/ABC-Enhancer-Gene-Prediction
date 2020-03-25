qcExpt <- function(expt, opt) {
  #Check for duplicate experiments
  #dupe <- any(duplicated(expt[, c("CellType","TargetGene","chr","start","end")] ))
  dupe <- any(duplicated(expt[, c("TargetGene","chr","start","end")] ))
  if (dupe) {
    print("Error: The experimental data file contains duplicate experiments!")
    stop()
  }
  
  #check to make sure regulated column contains TRUE/FALSE
  reg.vals <- sort(unique(expt[, get(opt$experimentalPositiveColumn)]))
  if (!(all(reg.vals %in% c(FALSE, TRUE)) | all(reg.vals %in% c(0, 1)))) {
    print("Error: The experimental data column must contain TRUE/FALSE")
    stop()
  }
  if (length(reg.vals) == 1) {
    print("Note: all values are either positives or negatives. Plotting code will fail, but merged prediction/experimental table will be output.")
  }
}

combineExptPred <- function(expt, pred, config, fill.missing=TRUE) {
  config <- fread(config)
  pred.gr <- with(pred, GRanges(paste0(chr,":",TargetGene), IRanges(start, end)))
  expt.gr <- with(expt, GRanges(paste0(chr,":",TargetGene), IRanges(start, end)))
  #pred.gr <- with(pred, GRanges(paste0(CellType,":",chr,":",TargetGene), IRanges(start, end)))
  #expt.gr <- with(expt, GRanges(paste0(CellType,":",chr,":",TargetGene), IRanges(start, end)))
  ovl <- GenomicRanges::findOverlaps(expt.gr, pred.gr)
  
  #Rename duplicated columns
  col.bool <- colnames(pred) %in% colnames(expt)
  colnames(pred)[col.bool] <- paste0(colnames(pred)[col.bool], ".from.predictions")
   
  #Merge predictions with experimental data
  merged <- cbind(expt[queryHits(ovl)], pred[subjectHits(ovl)])
  
  #Sometimes a perturbed element will overlap multiple model elements (eg in the case of a large deletion)
  #In these cases need to summarize, Eg sum ABC.Score across model elements overlapping the deletion
  #This requires a config file describing how each prediction column should be aggregated
  #agg.cols <- c("chr","start","end","TargetGene","class","CellType","Regulated","Significant","PctChange")
  agg.cols <- c("chr","start","end","TargetGene","class","Regulated","Significant","PctChange")

  if (!all(agg.cols %in% colnames(merged))) {
      print(paste0("Ignoring the following columns that were missing from the file: ", setdiff(agg.cols, colnames(merged))))
      agg.cols <- intersect(agg.cols, colnames(merged))
  }

  merged <- collapseEnhancersOverlappingMultiplePredictions(merged, config, agg.cols)
  
  #Experimental data missing predictions
  #A tested enhancer element may not have a prediction
  #For ABC this is typically the case if the tested element does not overlap a DHS peak.
  #In this case we need to fill the predictions table
  #However, this code doesn't work if the predictor did not make a prediction for this gene.
  expt.missing.predictions <- expt[setdiff(seq(nrow(expt)), queryHits(ovl)),]
  print("The following experimental data is not present in predictions file: ")
  print(expt.missing.predictions[, setdiff(agg.cols, "class")])
  if (fill.missing) {
    expt.missing.predictions <- fillMissingPredictions(expt.missing.predictions, config, agg.cols)
    cols.we.want <- c(agg.cols, config$pred.col)
    merged <- rbind(merged, expt.missing.predictions[, ..cols.we.want])
    print("Experimental data missing predictions filled. Will be considered in PR curve!")
    print(expt.missing.predictions[, ..cols.we.want])
  } else {
    print("Experimental data missing predictions ignored. Will not be considered in PR curve!")
  }

  return(merged)
}

fillMissingPredictions <- function(df, config, agg.cols) {
  for (ii in seq(nrow(config))) {
    df[, config$pred.col[[ii]]] <- config$fill.val[ii]
  }
  
  unk.cols <- setdiff(agg.cols, unique(c(colnames(df), config$pred.cols)))
  df[, unk.cols] <- "Merge:UNKNOWN"
  return(df)
}

collapseEnhancersOverlappingMultiplePredictions <- function(df, config, agg.cols) {
  list.for.agg <- as.list(df[, ..agg.cols])
  
  all.list <- mapply(function(pred.col, agg.func) aggregate(df[, ..pred.col], by = list.for.agg, FUN = agg.func), 
                     config$pred.col, config$agg.func, SIMPLIFY=F)
  full.result <- Reduce(function(df1, df2) merge(df1, df2, by = agg.cols), all.list)
  
  return(full.result)
}

prepForPlotting <- function(df) {
  df$scatterplot.color <- with(df, ifelse(Regulated, "Activating", ifelse(Significant, "Repressive", "Not Significant")))
  return(df)
}

makePlots <- function(merged, config, pos.col, outdir) {
  config <- fread(config)  
  pred.cols <- unique(unlist(lapply(config$pred.cols, function(s) {strsplit(s,",")[[1]]})))
  
  #Make scatter plots
  lapply(pred.cols, function(s) {makeScatterPlot(merged, s, "PctChange", outdir)})
  
  #Compute PR curve
  #Hack for linear distance
  merged$distance <- -1*merged$distance
  pr <- sapply(pred.cols, 
               function(s) {performance(prediction(unlist(merged[, get(s)]), unlist(merged[, ..pos.col])), measure="prec", x.measure="rec")})
  merged$distance <- -1*merged$distance
  pr.df <- pr2df(pr)
  makeAUCTable(pr, outdir)
  write.table(pr.df, file.path(outdir, "pr.curve.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  #Make PR curve plots
  pct.pos <- sum(unlist(merged[, ..pos.col]))/nrow(merged)
  for (ii in seq(nrow(config))) {
    makePRCurvePlot(pr.df, config$plot.name[[ii]], config$pred.cols[[ii]], outdir, pct.pos = pct.pos)
  }
}

makeScatterPlot <- function(df, x.col, y.col, outdir) {
  g <- ggplot(df,
              aes(x = get(x.col),
                  y = get(y.col),
                  color = scatterplot.color)) + 
    geom_point() + 
    scale_color_manual(values=c("Activating" = "red", 
                                "Repressive" = "blue", 
                                "Not Significant" = "gray")) + 
    labs(x = x.col, y = y.col, color = "")
  
  ggsave(file.path(outdir, paste0(x.col, ".", y.col, ".scatter.pdf")), g, device = "pdf")
}

makePRCurvePlot <- function(pr.df, plot.name, col.list, outdir, pct.pos) {
  col.list <- strsplit(as.character(col.list), ",")[[1]]
  pr.df <- subset(pr.df, pred.col %in% col.list)
  
  #separate boolean predictors from continuous predictors
  pr.cutoff <- by(pr.df, pr.df$pred.col, function(df) unique(df$alpha))
  boolean.predictors <- names(pr.cutoff)[unlist(lapply(pr.cutoff, function(s) identical(s, c(Inf, 1, 0))), use.names=F)]
  
  cont.pred <- subset(pr.df, !(pred.col %in% boolean.predictors))
  bool.pred <- subset(pr.df, pred.col %in% boolean.predictors)
  bool.pred <- subset(bool.pred, alpha == 1)
  
  g <- ggplot(cont.pred,
           aes(x = recall,
               y = precision,
               color = pred.col)) + 
      geom_line() + 
      labs(title = plot.name, color = "") + 
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
      geom_hline(yintercept = pct.pos, linetype = 2, color = 'black')
    
  if (nrow(bool.pred) > 0) {
    g <- g + geom_point(data = bool.pred,
                        size = 3)
  }
      
  ggsave(file.path(outdir, paste0(plot.name, ".pr.pdf")), g, device = "pdf")
}

pr2df <- function(pr) {
  #Convert a list of ROCR performance objects into a dataframe
  
  doOne <- function(this.pr) {
    df <- as.data.frame(list(
      alpha = this.pr@alpha.values[[1]],
      precision = this.pr@y.values[[1]],
      recall = this.pr@x.values[[1]]
    ))
    return(df)
  }
  
  pr.list <- lapply(pr, doOne)
  
  for (ii in seq(length(pr.list))) {
    pr.list[[ii]]$pred.col <- names(pr.list)[ii]
  }
  
  return(rbindlist(pr.list))
}

makeAUCTable <- function(pr, outdir) {
  auc <- lapply(pr, function(s) computeAUC(s@x.values[[1]], s@y.values[[1]]))
  write.table(t(as.data.frame(auc)), file.path(outdir, "pr.summary.txt"), sep='\t', quote = F, row.names = T, col.names = T)
}

computeAUC <- function(x.vals, y.vals) {
  good.idx <- which(!is.na(x.vals) & !is.na(y.vals))
  return(trapz(x.vals[good.idx], y.vals[good.idx]))
}

fread_ignore_empty <- function(f, ...) {
  tryCatch({
    return(fread(f, fill = TRUE, ...))
  }, error = function(e){
    print(paste0("Could not open file: ", f))
    return()
  })
}

fread_gz_ignore_empty <- function(f, ...) {
  tryCatch({
    return(fread(paste0("gunzip -c ", f), ...))
  }, error = function(e){
    print(paste0("Could not open file: ", f))
    return()
  })
}

smart_fread <- function(f, ...) {
  if (summary(file(f))$class == "gzfile") {
    out <- fread_gz_ignore_empty(f, ...)
  } else {
    out <- fread_ignore_empty(f, ...)
  }
  
  #con <- file(f)
  #on.exit(close(con), add = TRUE)
  closeAllConnections()
  return(out)
}

loadDelimited <- function(file.list) {
  data.list <- lapply(file.list, smart_fread)
  return(rbindlist(data.list))
}

loadFileString <- function(file.str, delim = ",") {
  file.list <- strsplit(file.str, split = delim)[[1]]
  return(loadDelimited(file.list))
}
