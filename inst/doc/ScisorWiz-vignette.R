## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ScisorWiz)

## ----ScisorWiz_AllInfo, eval=FALSE, echo=TRUE---------------------------------
#  
#  gencodeAnno <- system.file("extdata/", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  allInfoFile <- system.file("extdata/", "AllInfo.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata/", "userInput/celltypeFile",
#                           package = "ScisorWiz")
#  
#  ## Run command without plotting mismatches
#  ScisorWiz_AllInfo(gencodeAnno = gencodeAnno, AllInfoInput = allInfoFile,
#                    cellTypeFile = cTypeFile, gene = "Snap25", cluster = 1,
#                    ci = .05, outputDir = "extdata/outputDir/")
#  
#  ##___________________________________OR____________________________________##
#  
#  # If the optional mismatchFinder function hasn't been run, the mismatches file
#  # is located in extData/, rather than within the outputDir subdirectory
#  gencodeAnno <- system.file("extdata/", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  allInfoFile <- system.file("extdata/", "AllInfo.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata/", "userInput/celltypeFile",
#                           package = "ScisorWiz")
#  mismatches <- system.file("extdata/", "outputDir/Snap25.mismatches.txt.gz",
#                             package = "ScisorWiz")
#  
#  ## Run command with plotting mismatches
#  ScisorWiz_AllInfo(gencodeAnno = gencodeAnno, AllInfoInput = allInfoFile,
#                    cellTypeFile = cTypeFile, gene = "Snap25", cluster = 1,
#                    ci = .05, mismatchCutoff = .05,
#                    outputDir = "extdata/outputDir/", mismatchFile = mismatches)

## ----plot1, out.width = '60%', echo=F-----------------------------------------
knitr::include_graphics("../man/figures/Snap25_Isoform_Plot_noMis.pdf")

## ----ScisorWiz_2File, eval=FALSE, echo=TRUE-----------------------------------
#  
#  
#  ## Run command without plotting mismatches
#  ScisorWiz_2File(gencodeAnno = "gencodeAnnoFile.gz", gffInput = "CagePolyA.gff.gz",
#                  genesInput = "reads2genes.gz",cellTypeFile = "cellTypeFile_Snap25.tab",
#                  gene = "Snap25", ci = .05, outputDir = "extdata/outputDir/")
#  
#  ##___________________________________OR____________________________________##
#  
#  ## Run command with plotting mismatches
#  ScisorWiz_2File(gencodeAnno = "gencodeAnnoFile.gz",
#                  gffInput = "CagePolyA.gff.gz", genesInput = "reads2genes.gz",
#                  cellTypeFile = "cellTypeFile_Snap25.tab",gene = "Snap25",
#                  ci = .05, mismatchCutoff = .05, outputDir = "extdata/outputDir/",
#                  mismatchFile = "Snap25.mismatches.txt.gz")

## ----MismatchFinder, eval=FALSE, echo=TRUE------------------------------------
#  ## Run command
#  MismatchFinder(BAM = "sorted.bestperRead.mapping.bam", fasta = "mm10.fa",
#                 gencodeAnno = "gencode.vM21.annotation.gtf.gz", gene = "Snap25",
#                 outputDir = "extdata/outputDir/")

## ----plot2, out.width = '60%', echo=F-----------------------------------------
knitr::include_graphics("../man/figures/Snap25_Isoform_Plot.pdf")

## ----ScisorWiz_AllInfo_interactive, eval=FALSE, echo=TRUE---------------------
#  
#  gencodeAnno <- system.file("extdata/", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  allInfoFile <- system.file("extdata/", "AllInfo.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata/", "userInput/celltypeFile",
#                           package = "ScisorWiz")
#  
#  ## Run command without plotting mismatches
#  ScisorWiz_AllInfo(gencodeAnno = gencodeAnno, AllInfoInput = allInfoFile,
#                    cellTypeFile = cTypeFile, gene = "Snap25", cluster = 1,
#                    ci = .05, outputDir = "extdata/outputDir/", interactive = "y")

