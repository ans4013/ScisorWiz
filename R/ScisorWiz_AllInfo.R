#' ScisorWiz AllInfo File pipeline
#' @aliases ScisorWizAllInfo
#' @description This is the script to run the entire pipeline if the user enters
#' an AllInfo FIle as their data. AllInfo file is an output of the scisorseqr
#' pipeline. Uses the data from single-cell long-read RNA sequencing to
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
#' TSS site, 3 = PolyA site, 4 = intron chain, TSS site, and PolyA site)
#' @param ci User-specified confidence interval for reads to be considered
#' alternate exons. Is .05 by default.
#' @param outputDir User-specified directory to store all output from the
#' pipeline
#' @param mismatchFile Output of MismatchFinder function if used. Is NULL by
#' default.
#'
#' @return Plot visualizing isoform expression of the gene of interest among up
#' to 6 user-specified cell types
#'
#' @usage ScisorWiz_AllInfo(gencodeAnno, AllInfoInput, cellTypeFile, gene,
#' cluster, ci, outputDir, mismatchFile)
#'
#' @export

ScisorWiz_AllInfo <- function(gencodeAnno, AllInfoInput, cellTypeFile, gene,
                              cluster, ci=.05, outputDir, mismatchFile=NULL) {
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
    runR <- paste("Rscript", R_file, "pdf", plotName, 3, annoRemap,
                  cellTypeFilewithFileNames, orderFile, all5File, altExonsFile,
                  projectionRemapFile, gene, cluster, ci, plotOutput, SNVFile,
                  insertionsFile, deletionsFile)
  }
  else{
    runR <- paste("Rscript", R_file, "pdf", plotName, 3, annoRemap,
                cellTypeFilewithFileNames, orderFile, all5File, altExonsFile,
                projectionRemapFile, gene, cluster, ci, plotOutput)
  }
  system(runR)

  rm(list=ls())
}
