## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----demonstrateMan, include = TRUE, eval = FALSE-----------------------------
#  ?ScisorWiz_AllInfo

## ----ViewAllInfoFile, eval=FALSE, echo=TRUE-----------------------------------
#  library(ScisorWiz)
#  allInfoExampleFile <- system.file("extdata", "userInput/AllInfo", package = "ScisorWiz")
#  system(paste("cat", allInfoExampleFile))

## ----ScisorWiz_AllInfo_NoMis, eval=FALSE, echo=TRUE---------------------------
#  library(ScisorWiz)
#  
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  allInfoFile <- system.file("extdata", "userInput/AllInfo.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata", "userInput/cellTypeFile_Snap25.tab",
#                           package = "ScisorWiz")
#  
#  ## Get static plot without plotting mismatches
#  ScisorWiz_AllInfo(gencodeAnno = gencodeAnnoFile, AllInfoInput = allInfoFile,
#                    cellTypeFile = cTypeFile, gene = "Snap25", cluster = 1,
#                    ci = .05, outputDir = "outputDir/")

## ----plot1, fig.align="center", out.width = '50%', echo=F---------------------
knitr::include_graphics("../man/figures/Snap25_Isoform_Plot_testdata_AllInfo_NoMis.pdf")

## ----ScisorWiz_AllInfo_Mismatches, eval=FALSE, echo=TRUE----------------------
#  ##___________________________________OR____________________________________##
#  
#  library(ScisorWiz)
#  
#  # In order to retain a smaller overall size for the ScisorWizpackage, the
#  # example data mismatches file is located in extData/, rather than within the
#  # outputDir subdirectory where normal output from the first step would go.
#  
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  allInfoFile <- system.file("extdata", "userInput/AllInfo.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata", "userInput/cellTypeFile_Snap25.tab",
#                           package = "ScisorWiz")
#  
#  # If you have not run the MismatchFinder function (see Section 2 below), use the
#  # included mismatches file by running the following line of code:
#  mismatchesFile <- system.file("extdata", "Snap25.mismatches.txt.gz",
#                             package = "ScisorWiz")
#  
#  # If you have run the MismatchFinder function, use the output created by the
#  # function by running the following line of code:
#  mismatchesFile <- "outputDir/Snap25/Snap25.mismatches.txt.gz"
#  
#  ## Get static plot with plotting mismatches
#  ScisorWiz_AllInfo(gencodeAnno = gencodeAnnoFile, AllInfoInput = allInfoFile,
#                    cellTypeFile = cTypeFile, gene = "Snap25", cluster = 1,
#                    ci = .05, mismatchCutoff = .05,
#                    outputDir = "outputDir/", mismatchFile = mismatchesFile)

## ----plot2, fig.align="center", out.width = '50%', echo=F---------------------
knitr::include_graphics("../man/figures/Snap25_Isoform_Plot_testdata_AllInfo.pdf")

## ----scisorseqr, eval=FALSE, echo=TRUE----------------------------------------
#  
#  library(scisorseqr)
#  
#  print("++++++++ Step 1: Getting barcodes and filtering those reads out")
#  
#  fqFolderPath <- system.file("extdata/", "userInput/fqFolder/", package = "ScisorWiz")
#  ctAssignments <- system.file("extdata/", "userInput/barcode_cellType_assignments",
#                               package = "ScisorWiz")
#  
#  GetBarcodes(fqFolder=fqFolderPath,
#          BCClustAssignFile=ctAssignments,
#          chemistry='v2', filterReads=TRUE, numProcesses=12,
#          concatenate=FALSE)
#  
#  print("++++++++ Step 2: Aligning with minimap2")
#  mmProgPath <- '~/minimap2/minimap2' ## Please provide path to minimap2 aligner
#  genomeFa <- '~/genomes/M.musculus/mm10.fa' ## Please provide path to a mouse genome
#  
#  MMalign(fqFolder=fqFolderPath,mmProgPath,
#          refGenome=genomeFa,
#          numThreads=12)
#  
#  print("++++++++ Step 3: Map and filter function")
#  
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                                 package = "ScisorWiz")
#  ## Please provide a path to a directory containing one fa.gz file per chromosome
#  ## below
#  chr_fa_dir <- '~/genomes/M.musculus/mm10/chromFa/'
#  
#  # Below two arguments are optional. Can instead set filterFullLength=FALSE in
#  # the function
#  polyABed_path <- 'atlas.clusters_chr.mm10.2-0.bed.gz'
#  cageBed_path <- 'mm10_fair+new_CAGE_peaks_phase1and2.bed.gz'
#  
#  MapAndFilter(numThreads=12, filterFullLength=TRUE,
#          polyABed=polyABed_path,
#          cageBed=cageBed_path,
#          annoGZ=gencodeAnnoFile,
#          seqDir=chr_fa_dir, genomeVersion='mm10')
#  
#  print("++++++++ Step 4: Getting All-Info files")
#  ## If default parameters for Step 1 are used, the the barcodeOutput file will be
#  ##autogenenerated in the "OutputFiltered" folder
#  
#  InfoPerLongRead(barcodeOutputFile='OutputFiltered/FilteredDeconvBC_P7HIPP_subset.csv',
#                  mapAndFilterOut='LRProcessingOutput/', minTimesIsoObserve=1,
#                  rmTmpFolder=FALSE)
#  
#  print("+++++++ All done!")

## ----ScisorWiz_2File_NoMis, eval=FALSE, echo=TRUE-----------------------------
#  
#  ## Run command without plotting mismatches
#  library(ScisorWiz)
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  gffFile <- system.file("extdata", "userInput/readinfo.gff.gz", package = "ScisorWiz")
#  genesFile <- system.file("extdata", "userInput/readstogenes.genes.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata", "userInput/cellTypeFile_Snap25.tab",
#                           package = "ScisorWiz")
#  
#  ScisorWiz_2File(gencodeAnno = gencodeAnnoFile, gffInput = gffFile,
#                  genesInput = genesFile, cellTypeFile = cTypeFile,
#                  gene = "Snap25", ci = .05, outputDir = "outputDir/")

## ----plot3, fig.align="center", out.width = '50%', echo=F---------------------
knitr::include_graphics("../man/figures/Snap25_Isoform_Plot_testdata_2File_NoMis.pdf")

## ----ScisorWiz_2File_Mismatches, eval=FALSE, echo=TRUE------------------------
#  ##___________________________________OR____________________________________##
#  
#  ## Run command with plotting mismatches
#  library(ScisorWiz)
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  gffFile <- system.file("extdata", "userInput/readinfo.gff.gz", package = "ScisorWiz")
#  genesFile <- system.file("extdata", "userInput/readstogenes.genes.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata", "userInput/cellTypeFile_Snap25.tab",
#                           package = "ScisorWiz")
#  
#  # If you have not run the MismatchFinder function (see Section 2 below), use the
#  # included mismatches file by running the following line of code:
#  mismatchesFile <- system.file("extdata", "Snap25.mismatches.txt.gz",
#                             package = "ScisorWiz")
#  
#  # If you have run the MismatchFinder function, use the output created by the
#  # function by running the following line of code:
#  mismatchesFile <- "outputDir/Snap25/Snap25.mismatches.txt.gz"
#  
#  ScisorWiz_2File(gencodeAnno = gencodeAnnoFile,
#                  gffInput = gffFile, genesInput = genesFile,
#                  cellTypeFile = cTypeFile, gene = "Snap25",
#                  ci = .05, mismatchCutoff = .05, outputDir = "outputDir/",
#                  mismatchFile = mismatchesFile)

## ----plot4, fig.align="center", out.width = '50%', echo=F---------------------
knitr::include_graphics("../man/figures/Snap25_Isoform_Plot_testdata_2File_Mis.pdf")

## ----MismatchFinder, eval=FALSE, echo=TRUE------------------------------------
#  ## Run command
#  library(ScisorWiz)
#  
#  # Optional: download the reference genome fasta
#  system(paste("wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
#               "&& gunzip mm10.fa.gz"))
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                                 package = "ScisorWiz")
#  bam <- system.file("extdata", "userInput/fqFolder/sorted.Snap25.bam",
#                     package = "ScisorWiz")
#  
#  # Find mismatches for Snap25 gene data
#  MismatchFinder(BAM = bam, fasta = "mm10.fa", gencodeAnno = gencodeAnnoFile,
#                 gene = "Snap25", outputDir = "outputDir/")

## ----ScisorWiz_AllInfo_interactive, eval=FALSE, echo=TRUE---------------------
#  library(ScisorWiz)
#  
#  gencodeAnnoFile <- system.file("extdata", "gencode.vM21.annotation.gtf.gz",
#                             package = "ScisorWiz")
#  allInfoFile <- system.file("extdata", "userInput/AllInfo.gz", package = "ScisorWiz")
#  cTypeFile <- system.file("extdata", "userInput/cellTypeFile_Snap25.tab",
#                           package = "ScisorWiz")
#  
#  ## Get interactive plot without plotting mismatches
#  ScisorWiz_AllInfo(gencodeAnno = gencodeAnnoFile, AllInfoInput = allInfoFile,
#                    cellTypeFile = cTypeFile, gene = "Snap25", cluster = 1,
#                    ci = .05, outputDir = "outputDir/", interactive = "y")

