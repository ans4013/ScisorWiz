#' ScisorWiz AllInfo File pipeline
#' @aliases ScisorWizAllInfo
#' @description This is the script to run the entire pipeline if the user inputs
#' an AllInfo File as their data. AllInfo file is an output of the scisorseqr
#' pipeline. Uses data from single-cell long-read RNA sequencing to visualize
#' the isoform expression patterns of a single gene across any number of cell
#' types. Data is clustered in one of 4 ways: by intron chain of the reads, by
#' TSS site, by PolyA site, or by all three clustering methods. User has option
#' to include a mismatch file by first running necessary files through the
#' MismatchFinder function. User can choose to create an interactive plot for
#' exploratory purposes by setting the "interactive" argument equal to "y".
#'
#' @param gencodeAnno Gencode annotation file from which data files are created
#' @param AllInfoInput File generated from the scisorseqr package that houses
#' all data in a uniquely organized fashion: readID, geneID, celltype, barcode,
#' intron chain, etc.
#' @param cellTypeFile Tab-separated file housing user specified cell types and
#' desired color of reads of those cell types for the output plot
#' @param gene User-specified gene of interest
#' @param cluster User-specified clustering method (1 = intron chain, 2 =
#' TSS site, 3 = PolyA site, 4 = intron chain, TSS site, and PolyA site).
#' Default is 1 (intron chain)
#' @param ci User-specified confidence interval for reads to be considered
#' alternate exons. Default is .05.
#' @param mismatchCutoff User-specified cutoff for SNV inclusion rate. Default
#' is .05
#' @param outputDir User-specified directory to store all output from the
#' pipeline
#' @param mismatchFile Output of MismatchFinder function if used. Default is
#' NULL.
# @param zoom Yes (y) or No (n) for zooming into a user-specified window on the
# plot. Default is No (n).
#' @param interactive Yes (y) or No (n) for creating an interactive html plot
#' file using a jpg produced from user input. Default is No (n).
#'
#' @return Plot visualizing isoform expression of a gene of interest among any
#'  number of user-specified cell types
#'
#' @usage ScisorWiz_AllInfo(gencodeAnno, AllInfoInput, cellTypeFile, gene,
#' cluster, ci, mismatchCutoff, outputDir, mismatchFile, interactive)
#'
#' @export

ScisorWiz_AllInfo <- function(gencodeAnno, AllInfoInput, cellTypeFile, gene,
                              cluster=1, ci=.05, mismatchCutoff=.05, outputDir,
                              mismatchFile=NULL, zoom="n", interactive = "n") {
  cat("================= Handling arguments =================\n")

  dir.create(outputDir, recursive = T)

  if(!file.exists(gencodeAnno)){
    warning('Gencode file path not valid. Please retry.')
  }

  if(!file.exists(AllInfoInput)){
    warning('AllInfo file path not valid. Please retry.')
  }

  if(!file.exists(cellTypeFile)){
    warning('Cell type file path not valid. Please retry.')
  }

  geneOutput <- paste0(outputDir, gene, "/")
  plotOutput <- paste0(outputDir,"Plots/")
  genePlotOutput <- paste0(plotOutput,gene,"/")

  if(!dir.exists(geneOutput)){
    dir.create(geneOutput)
  }

  if(!dir.exists(plotOutput)){
    dir.create(plotOutput)
  }

  if(!dir.exists(genePlotOutput)){
    dir.create(genePlotOutput)
  }

  if(!is.null(mismatchFile)){
    if(!file.exists(mismatchFile)){
    warning('Mismatch file path not valid. Please retry.')
    }
  }

  cat(paste("Output Directory:", geneOutput, "\n"))
  cat(paste("Plot Output Directory:", genePlotOutput, "\n"))
  cat(paste("Annotation File:", gencodeAnno, "\n"))
  cat(paste("Data File:", AllInfoInput, "\n"))
  cat(paste("Cell Type File:", cellTypeFile, "\n"))
  cat(paste("Gene of Interest:", gene, "\n"))
  if(!is.null(mismatchFile)){
    cat(paste("Mismatch File:", mismatchFile, "\n"))
  }

  cat("================= Running python script ==================\n")
  py_file <- system.file("python", "ClusterByIsoform_AllInfoInput.py",
                         package = "ScisorWiz")
  if(!is.null(mismatchFile)){
    runPy <- paste("python3", py_file, gencodeAnno, AllInfoInput, cellTypeFile,
                   gene, ci, outputDir, mismatchFile)
  }
  else{
    runPy <- paste("python3", py_file, gencodeAnno, AllInfoInput, cellTypeFile,
                   gene, ci, outputDir)
  }
  system(runPy)

  plotName <- paste0(gene, "_Isoform_Plot")
  annoRemap <- paste0(geneOutput, gene, ".anno_remap.gtf.gz")
  cellTypeFilewithFileNames <- paste0(geneOutput, gene,
                                      ".cellTypeFileWithFileNames.tab")
  orderFile <- paste0(geneOutput, gene, ".order.tab.gz")
  all5File <- paste0(geneOutput, gene,
                     ".all5DatasetsSpanningRegion.remap.gff.gz")
  altExonsFile <- paste0(geneOutput, gene, ".altExons.tab")
  projectionRemapFile <- paste0(geneOutput, gene, ".projection.tab.remap.gz")

  cat("==================== Running R script ====================\n")
  R_file <- system.file("RScript", "PlotIsoforms.r", package = "ScisorWiz")
  if(!is.null(mismatchFile)){
    SNVFile <- paste0(geneOutput, gene, ".SNVs.tab")
    insertionsFile <- paste0(geneOutput, gene, ".insertions.tab")
    deletionsFile <- paste0(geneOutput, gene, ".deletions.tab")
    if(zoom == "y"){
      cat("Please enter exon number for left side of zoom window:\n")
      windowStart <- scan(what = integer)
      cat("Please enter exon number for right side of zoom window:\n")
      windowEnd <- scan(what = integer)

      runR <- paste("Rscript", R_file, interactive, plotName, annoRemap,
                    cellTypeFilewithFileNames, orderFile, all5File, altExonsFile,
                    projectionRemapFile, gene, cluster, ci, mismatchCutoff,
                    plotOutput, SNVFile, insertionsFile, deletionsFile,
                    windowStart, windowEnd)
    }
    else {
      runR <- paste("Rscript", R_file, interactive, plotName, annoRemap,
                    cellTypeFilewithFileNames, orderFile, all5File, altExonsFile,
                    projectionRemapFile, gene, cluster, ci, mismatchCutoff,
                    plotOutput, SNVFile, insertionsFile, deletionsFile)
    }
  }
  else{
    if(zoom == "y"){
      cat("Please enter exon number for left side of zoom window:\n")
      windowStart <- scan(what = integer)
      cat("Please enter exon number for right side of zoom window:\n")
      windowEnd <- scan(what = integer)

      runR <- paste("Rscript", R_file, interactive, plotName, annoRemap,
                    cellTypeFilewithFileNames, orderFile, all5File, altExonsFile,
                    projectionRemapFile, gene, cluster, ci, mismatchCutoff,
                    plotOutput, windowStart, windowEnd)
    }
    else{
      runR <- paste("Rscript", R_file, interactive, plotName, annoRemap,
                    cellTypeFilewithFileNames, orderFile, all5File, altExonsFile,
                    projectionRemapFile, gene, cluster, ci, mismatchCutoff,
                    plotOutput)
    }
  }
  system(runR)

  if(interactive == "y"){
    cat("ENTERING INTERACTIVE PLOT\n")
    interactiveScript <- system.file("RScript", "interactivePlot.R", package = "ScisorWiz")
    plotPath <- paste0(genePlotOutput, plotName, ".jpg")
    htmlPath <- paste0(genePlotOutput, plotName, ".html")
    cat(paste(plotPath, "\n"))
    cat(paste(htmlPath, "\n"))
    runInteractive <- paste("Rscript", interactiveScript, plotPath, htmlPath)

    system(runInteractive)
  }

  rm(list=ls())
}
