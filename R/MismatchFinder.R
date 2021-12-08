#' ScisorWiz Mismatch Finder Function
#' @aliases MismatchFinder
#' @description Cross references data from .bam file against a reference genome
#' to determine mismatches in the genetic code. These mismatches include
#' Single Nucleotide Polymorphisms (SNPs), nucleotide insertions, and
#' nucleotide deletions. These mismatches are listed per read to be later
#' included in the ScisorWiz Pipeline function.
#'
#' @param BAM A sorted .bam file for the dataset of interest
#' @param fasta A reference genome .fasta file for cross referencing data in
#' .bam file to reference genome
#' @param gencodeAnno Gencode annotation file from which data files are created
#' @param gene User-specified gene of interest
#' @param outputDir User-specified directory to store output from the function
#'
#' @return Mismatch file with positions of SNPs, insertions, and deletions
#' organized by read
#'
#' @usage MismatchFinder(BAM, fasta, gencodeAnno, gene, outputDir)
#'
#' @export

MismatchFinder <- function(BAM, fasta, gencodeAnno, gene, outputDir) {
  print("================= Handling arguments =================")

  geneOutput <- paste0(outputDir, gene, "/")
  dir.create(geneOutput, recursive = T)

  if(!file.exists(BAM)){
    warning('.bam file path not valid. Please retry.')
  }

  if(!file.exists(fasta)){
    warning('Reference .fasta file path not valid. Please retry.')
  }

  if(!file.exists(gencodeAnno)){
    warning('Reference GENCODE annotation file path not valid. Please retry.')
  }

  print("===== Finding chromosome, start, and end of gene =====")

  startEnd_file <- system.file("python", "start_end_finder.py", package="ScisorWiz")
  runSEFinder <- paste("python3", startEnd_file, "--gencodeAnno", gencodeAnno,
                       "--gene", gene, "--outputDir", outputDir)
  #geneInfo <- vector(mode = "list", length = 3)
  system(runSEFinder)

  geneInfoFile <- paste0(geneOutput, gene, ".info.tab")
  geneInfo=read.table(geneInfoFile, sep = "\t")

  chrom = geneInfo$V1
  start = geneInfo$V2
  end = geneInfo$V3

  print(paste(".bam file =", BAM))
  print(paste(".fasta file =", fasta))
  print(paste("Chromosome =", chrom))
  print(paste("Gene start =", start))
  print(paste("Gene end =", end))
  print(paste("Output directory =", outputDir))

  print("=============== Determining Mismatches ===============")
  mismatch_file <- system.file("python", "mismatch_finder.py", package = "ScisorWiz")
  runMismatch <- paste("python3", mismatch_file, "--BAM", BAM, "--fasta", fasta, "--chrom",
                 chrom, "--start", start, "--end", end, "--gene", gene, "--outPrefix", geneOutput)
  system(runMismatch)

  rm(list=ls())
}
