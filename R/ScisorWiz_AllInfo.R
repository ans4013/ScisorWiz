#' ScisorWiz AllInfo File pipeline
#' @aliases ScisorWizAllInfo
#' @description This is the script to run the entire pipeline if the user enters
#' an AllInfo File as their data. AllInfo file is an output of the scisorseqr
#' pipeline. Uses the data from si  ngle-cell long-read RNA sequencing to
#' visualize the isoform expression patterns of a single gene across up to six
#' celltypes. The data is clustered in one of 4 ways: by intron chain of the
#' reads, by TSS site, by PolyA site, or by PolyA and TSS sites. User has option
#' to include a mismatch file by first running necessary files through the
#' MismatchFinder function.
#'
#' @param gencodeAnno Gencode annotation file from which data files are created
#' @param AllInfoInput File generated from the scisorseqr package that houses
#' all data in a uniquely organized fashion: readID, geneID, celltype, barcode,
#' intron chain, etc.
#' @param cellTypeFile File housing user specified celltypes and desired color
#' of reads of those celltypes for the output plot
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
#' @param interactive Yes (y) or No (n) for creating an interactive plot session
#' using the pdf produced from user input. Default is No (n).
#'
#' @return Plot visualizing isoform expression of the gene of interest among up
#' to 6 user-specified cell types
#'
#' @usage ScisorWiz_AllInfo(gencodeAnno, AllInfoInput, cellTypeFile, gene,
#' cluster, ci, mismatchCutoff, outputDir, mismatchFile, interactive)
#'
#' @export

ScisorWiz_AllInfo <- function(gencodeAnno, AllInfoInput, cellTypeFile, gene,
                              cluster=1, ci=.05, mismatchCutoff=.05, outputDir,
                              mismatchFile=NULL, zoom="n", interactive = "n") {
  print("================= Handling arguments =================")

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

  print(paste("Output Directory:", geneOutput))
  print(paste("Plot Output Directory:", genePlotOutput))
  print(paste("Annotation File:", gencodeAnno))
  print(paste("Data File:", AllInfoInput))
  print(paste("Cell Type File:", cellTypeFile))
  print(paste("Gene of Interest:", gene))
  if(!is.null(mismatchFile)){
    print(paste("Mismatch File:", mismatchFile))
  }

  print("=============== Running python script ================")
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

  print("================== Running R script ==================")
  R_file <- system.file("RScript", "PlotIsoforms.r", package = "ScisorWiz")
  if(!is.null(mismatchFile)){
    SNVFile <- paste0(geneOutput, gene, ".SNVs.tab")
    insertionsFile <- paste0(geneOutput, gene, ".insertions.tab")
    deletionsFile <- paste0(geneOutput, gene, ".deletions.tab")
    if(zoom == "y"){
      cat("Please enter exon number for left side of zoom window:")
      windowStart <- scan(what = integer)
      cat("Please enter exon number for right side of zoom window:")
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
      cat("Please enter exon number for left side of zoom window:")
      windowStart <- scan(what = integer)
      cat("Please enter exon number for right side of zoom window:")
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
    print("ENTERING INTERACTIVE PLOT")
    interactiveScript <- system.file("RScript", "interactivePlot.R", package = "ScisorWiz")
    plotPath <- paste0(genePlotOutput, plotName, ".jpg")
    htmlPath <- paste0(genePlotOutput, plotName, ".html")
    print(plotPath)
    print(htmlPath)
    runInteractive <- paste("Rscript", interactiveScript, plotPath, htmlPath)

    system(runInteractive)
  }

  rm(list=ls())
}
