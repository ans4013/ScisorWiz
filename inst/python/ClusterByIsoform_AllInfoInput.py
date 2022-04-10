#!/usr/bin/env python

###############################################################################
# Analyze data - Filter for gene, filter by cell type, remap gene locations,
# separate out mismatches
#   - Input: GENCODE annotation file, AllInfo file, cell type file, gene,
#            alt exon inclusion interval (ci), output directory, mismatch file
#   - Output: Many files used in the R script to visualize data
#
# Command: ./ClusterByIsoform_AllInfoInput.py GENCODE_Annotation.gtf AllInfo.gz
#                               cellTypeFile.tab, gene_name, ci, outputDir,
#                               (geneName).mismatchfile.txt.gz
#
# by Alexander N. Stein
################################################################################

import sys
import pandas as pd
import gzip
import os.path
from os import path
import numpy as np
import re


# Step 1 - Manage arguments and check for errors
print("Step 1 of 15 - Managing command line arguments")

annotationFile = sys.argv[1]
if not path.exists(sys.argv[1]):
    print("ERROR: ", sys.argv[1], " does not exist.")

allInfoFile = sys.argv[2]
if not path.exists(sys.argv[2]):
    print("ERROR: ", sys.argv[2], " does not exist.")

cellTypeFile = sys.argv[3]
if not path.exists(sys.argv[3]):
    print("ERROR: ", sys.argv[3], " does not exist.")

gene = sys.argv[4]

CI = sys.argv[5]

outputDir = sys.argv[6]
geneOutputDir = outputDir + gene

if len(sys.argv) > 7:
    mismatchFile = sys.argv[7]
    if not path.exists(sys.argv[7]):
        print("ERROR: ", sys.argv[7], " does not exist.")
else:
    mismatchFile = None


# Step 2 - Organize data
# Step 2a - Get annotation for the gene
print("Step 2 of 15 - Getting annotation for", gene, "-- May take some time depending on file size")
annoFile = gzip.open(annotationFile, "rb")

geneMatch = "\"" + gene + "\""
annoContents = []
for line in annoFile:
    # Decode data from gzipped version
    line = line.decode()
    columns = line.rstrip().split()
    for i in range(9, len(columns)):
	# First filter for match in column 3 for "exon" or "CDS", then filter
	#     the matched row for the gene of interest. If match, append row
	#     to array
        if "exon" in columns[2] or "CDS" in columns[2]:
            if "gene_name" in columns[i] and len(re.findall(geneMatch, columns[i+1])) > 0:
                annoContents.append(columns)

annoFile.close()

# Create output file path
outputAnnoGTFFile = gene + ".anno.gtf.gz"
outputPath1 = geneOutputDir + "/" + outputAnnoGTFFile
annoOutput = gzip.open(outputPath1, "w")

# Write to output file and encode for gzipping
for i in annoContents:
    for element in i:
        annoOutput.write(str(element).encode())
        annoOutput.write("\t".encode())
    annoOutput.write("\n".encode())

annoOutput.close()


# Get start, end, chromosome and strand for the gene
print("Step 3 of 15 - Getting start, end, chromosome and strand for", gene, "gene")
minimum = 1000000000
maximum = -1

# Find start (minimum) and end (maximum) of all reads for the gene of interest
for i in annoContents:
    if int(i[3]) < minimum:
        minimum = int(i[3])
    if int(i[4]) > maximum:
        maximum = int(i[4])

if len(annoContents[0][0]) == 0:
    raise Error("Gene name not found in dataset. Please check spelling and retry.")
else:
    chromosome = annoContents[0][0]

strand = annoContents[0][6]

ucscHeaderInfo = "browser position ", str(chromosome), ":", str(minimum), "-", str(maximum)

print(chromosome, " : ", minimum, "-", maximum, " : ", strand)


# Get all reads from all datasets and identify TSS and PolyA sites
print("Step 4 of 15 - Getting all reads from all datasets and identifying TSS/PolyA sites -- May take some time depending on file size")

## Grab geneID for the reads. This step assumes same geneID for all rows in
##     annoContents array
for i in range(len(annoContents[0])):
    if "gene_id" in annoContents[0][i]:
        geneID = annoContents[0][i+1].strip("\";")

# Open relevant files
if allInfoFile[-3:len(allInfoFile)] == ".gz":
    print("File is gzipped")
    allInfo = [read.strip('\n') for read in gzip.open(allInfoFile, 'rt', encoding = 'utf-8').readlines()]
else:
    print("File is not gzipped")
    allInfo = [read.strip('\n') for read in open(allInfoFile, 'rt', encoding = 'utf-8').readlines()]

# Get info of reads from the AllInfo file where geneID matches
IDs = []
for line in allInfo:
    columns = line.strip().split()
    if geneID in columns[1]:
        readID = columns[2] + ":" + columns[3] + ":" + columns[0]
        if len(columns) > 9:
            exonChain = columns[8].split(";%;")
        else:
            exonChain = columns[6].split(";%;")
        intronChain = columns[5].split("&&")
        if len(intronChain) > 2:
            tss = intronChain[0].split("_")
            polyA = intronChain[2].split("_")
        for i in exonChain:
            if len(i) > 0:
                read = []
                readInfo = i.strip().split("_")
                chrom = readInfo[0]
                start = readInfo[1]
                end = readInfo[2]
                strand = readInfo[3]
                
                read.append(chrom)
                read.append("ScisorWiz")
                read.append("cDNA_match")
                read.append(start)
                read.append(end)
                read.append(".")
                read.append(strand)
                read.append(".")
                read.append("read_id")
                read.append(readID)
                if len(intronChain) > 2:
                    read.append(tss[1])
                    read.append(tss[2])
                    read.append(polyA[1])
                    read.append(polyA[2])
                else:
                    for i in range(1,5):
                        read.append("None")
                read.append(columns[5])
                IDs.append(read)

# Create output file path
outputSpanningRegion = gene + ".all5DatasetsSpanningRegion.gff.gz"
outputPath2 = geneOutputDir + "/" + outputSpanningRegion
spanningRegionOutput = gzip.open(outputPath2, "w")

# Write to output file and encode for gzipping
geneIDMatch = []
for line in IDs:
    for element in line:
        spanningRegionOutput.write(str(element).encode())
        spanningRegionOutput.write("\t".encode())
    spanningRegionOutput.write("\n".encode())
    geneIDMatch.append(line)

spanningRegionOutput.close()


# Order by isoform
print("Step 5 of 15 - Ordering by isoform")
struc = []
intSt = []
info = ""

# Gather information for all reads with the same read ID which were sorted into
#     geneIDMatch in Step 4
for line in geneIDMatch:
    ID = line[9]
    
    # intSt will be empty to begin so first read will go to else function. If this
    #     ID already exists in intSt, then the info variable will add the read
    #     information associated with this ID
    if ID in intSt:
        end = int(line[3]) - 1
        info = info + ";%;" + line[0] + "_" + str(start) + "_" + str(end) + "_" + line[6]
        start = int(line[4]) + 1
    # if ID not found in intSt and intSt is not empty, append info to it, append
    #     that to struc then empty both intSt and info to restart for new read ID
    else:
        if len(intSt) != 0:
            intSt.append(info)
            struc.append(intSt)
            intSt = []
            info = ""
        intSt.append(ID) 
        start = int(line[4]) + 1
# Ensure last read ID is also appended to final array
intSt.append(info)
struc.append(intSt)

def myFunc(e):
    return e[1]
struc.sort(key=myFunc)

# Create output file path
outputOrderFile = gene + ".order.tab.gz"
outputPath3 = geneOutputDir + "/" + outputOrderFile
orderOutput = gzip.open(outputPath3, "w")

# Write struc to output file
for line in struc:
    for element in line:
        orderOutput.write(str(element).encode())
        orderOutput.write("\t".encode())
    orderOutput.write("\n".encode())

orderOutput.close()


# Make projection
print("Step 6 of 15 - Making projection")

# Find start (minimum) and end (maximum) for all reads sorted from Step 4
for i in geneIDMatch:
    if int(i[3]) < minimum:  
        minimum = int(i[3])
    if int(i[4]) > maximum:
        maximum = int(i[4])

mini = 1000000000
maxi = 0
p = [None] * (maximum + 1)
h = ["exonic", "intronic"]
proj = []
projection = []
element = 0
for i in annoContents:
    if "exon" in i[2] or "cDNA_match" in i[2]:
        mini = min(mini, int(i[3]))
        maxi = max(maxi, int(i[4]))
	# p is filled with 0's for all elements which fall within the start and
	#     end points of the exons
        for j in range(int(i[3]), int(i[4]) + 1):
            p[j] = 0

for i in geneIDMatch:
#    columns = i.strip().split()
    if "exon" in i[2] or "cDNA_match" in i[2]:
        mini = min(mini, int(i[3]))
        maxi = max(maxi, int(i[4]))
	# p is filled with 0's for all elements which fall within the start and
	#     end points of the exons
        for j in range(int(i[3]), int(i[4]) + 1):
            p[j] = 0

# p is filled with 1's for all elements which fall outside of the start and end
#     points of the exons
for k in range(mini, maxi + 1):
    if p[k] == None:
        p[k] = 1
    lastStart = mini

# proj builds the array which becomes an element of final projection array.
#     h simply tells you whether the read is intronic (1) or exonic (0)
for m in range(mini + 1, maxi + 1):
    if p[m] != p[m-1]:
        x = m - 1
        proj = [chromosome, lastStart, x, h[p[m-1]]]
        projection.append(proj)
        lastStart = m 

# Make sure that the last read is recorded in the final array
proj = [chromosome, lastStart, maxi, h[p[maxi-1]]]
projection.append(proj)

# Create output file path
outputProjFile = gene + ".projection.tab.gz"
outputPath4 = geneOutputDir + "/" + outputProjFile
projOutput = gzip.open(outputPath4, "w")

# Write to output file
for line in projection:
    for element in line:
        projOutput.write(str(element).encode())
        projOutput.write("\t".encode())
    projOutput.write("\n".encode())

projOutput.close()


# Remap projection to new coordinates
print("Step 7 of 15 - Remapping projection to new coordinates")

mini = 0
maxi = 0
remap = []

# Converts the old start and end coordinates of exons to relative coordinates
# in a smaller range
for line in projection:
    if "exonic" in line[3]:
        mini = maxi + 1
        maxi = mini + (int(line[2]) - int(line[1]) + 1) - 1
        line.append(str(mini))
        line.append(str(maxi))
        remap.append(line)

# Makes all introns 100 points long on new range
    if "intronic" in line[3]:
        mini = maxi + 1
        maxi = mini + 99
        line.append(str(mini))
        line.append(str(maxi))
        remap.append(line)

# Create output file path
remapFile = gene + ".projection.tab.remap.gz"
outputPath5 = geneOutputDir + "/" + remapFile
remapOutput = gzip.open(outputPath5, "w")

# Write to output file
for line in remap:
    for element in line:
        remapOutput.write(str(element).encode())
        remapOutput.write("\t".encode())
    remapOutput.write("\n".encode())

remapOutput.close()


# Remap reads to new coordinates for each cell type
print("Step 8 of 15 - Remapping reads to new coordinates for each cell type")

cellType = open(cellTypeFile, "rb")

cellTypeReads = []
# Runs for each specified cellType
for cType in cellType:
    readRemap = []
    cType = cType.decode()
    cType = cType.split()
    
    # Search for cellType matches throughout geneIDMatch array
    for line in geneIDMatch:
        ID = line[9]
        if cType[0] in ID:
            for i in remap:
		# Only interested in exonic sections
                if "exonic" in i[3]:
		    # If read does not fall within the start and end of an exonic section, move on
                    if int(line[3]) < int(i[1]) or int(line[4]) > int(i[2]):
                        continue
                    else:
			# If read falls within exonic section, recalc start and end point, and append to new readRemap
                        remapRead = []
                        remapRead.append(line[0])
                        remapRead.append(line[1])
                        remapRead.append(line[2])
                        x = int(line[3]) - int(i[1])
                        y = int(i[2]) - int(line[4])
                        staDis = int(i[4]) + x
                        endDis = int(i[5]) - y
                        remapRead.append(staDis)
                        remapRead.append(endDis)
                        remapRead.append(line[5])
                        remapRead.append(line[6])
                        remapRead.append(line[7])
                        remapRead.append(ID)
                        remapRead.append(line[10])
                        remapRead.append(line[11])
                        remapRead.append(line[12])
                        remapRead.append(line[13])
                        remapRead.append(line[14])
                        #remapRead.append(line[15])
                        readRemap.append(remapRead)
                        
    # Ensure naming convention works for a file. No "/" in cell type name
    if "/" in cType[0]:
        cType[0] = cType[0].replace("/", "_")
    # Create output path for cell type specific remap file
    remapReadFile = gene + "." + cType[0] + ".remap.gtf.gz"
    outputPath6 = geneOutputDir + "/" + remapReadFile
    remapReadOutput = gzip.open(outputPath6, "w")
    
    # Write to output file
    for line in readRemap:
        cellTypeReads.append(line)
        for element in line:
            remapReadOutput.write(str(element).encode())
            remapReadOutput.write("\t".encode())
        remapReadOutput.write("\n".encode())

    remapReadOutput.close()


# Find alternative exons
print("Step 9 of 15 - Finding alternative exons")

altExons = []
exIncl = 0
exTot = 0
exonCount = {}
internal = {}
readStart = {}
readEnd = {}
countKeys = exonCount.keys()

for line in cellTypeReads:
    if line[8] in countKeys:
        exonCount[line[8]].append(line[0] + "_" + str(line[3]) + "_" + str(line[4]) + "_" + line[6])
    else:
        exonCount[line[8]] = []
        exonCount[line[8]].append(line[0] + "_" + str(line[3]) + "_" + str(line[4]) + "_" + line[6])

for read in countKeys:
    startInfo = 1000000000
    endInfo = 0
    length = len(exonCount[read])
    for info in range(0,length):
        data = exonCount[read][info].split("_")
        if int(data[1]) < startInfo:
            startInfo = int(data[1])
        if int(data[2]) > endInfo:
            endInfo = int(data[2])
    readStart[read] = startInfo
    readEnd[read] = endInfo

    for i in range(1,length-1):
        internal[exonCount[read][i]] = []

intKeys = internal.keys()
for exon in intKeys:
    exIncl = 0
    exTot = 0
    data = exon.split("_")

    for read in countKeys:
        length = len(exonCount[read])
        if readStart[read] < int(data[1]) and readEnd[read] > int(data[2]):
            exTot = exTot + 1
        for i in range(1,length-1):
            if exonCount[read][i] in exon:
                exIncl = exIncl + 1

    if exTot != 0:
      psi = exIncl/exTot
    else:
      psi = 0
    internal[exon].append(exIncl)
    internal[exon].append(exTot)
    internal[exon].append(psi)

altExons = []
for exon in intKeys:
    if internal[exon][2] > float(CI) and internal[exon][2] < (1 - float(CI)):
        altExons.append(exon)

if len(altExons) < 1:
    altExons.append(0)

altExonFile = gene + ".altExons.tab"
outputPath7 = geneOutputDir + "/" + altExonFile
altExonOutput = open(outputPath7, "w")

for exon in altExons:
    altExonOutput.write(str(exon))
    altExonOutput.write("\n")

altExonOutput.close()


# Remap geneIDMatch to new coordinates
print("Step 10 of 15 - Remapping all5DatasetsSpanningRegion to the new coordinates.")
geneIDRemap = []
for line in geneIDMatch:
    ID = line[9]
    for i in remap:
        if "exonic" in i[3]:
            if int(line[3]) < int(i[1]) or int(line[4]) > int(i[2]):
                continue
            else:
                remapGeneID = []
                remapGeneID.append(line[0])
                remapGeneID.append(line[1])
                remapGeneID.append(line[2])
                x = int(line[3]) - int(i[1])
                y = int(i[2]) - int(line[4])
                staDis = int(i[4]) + x
                endDis = int(i[5]) - y
                remapGeneID.append(staDis)
                remapGeneID.append(endDis)
                remapGeneID.append(line[5])
                remapGeneID.append(line[6])
                remapGeneID.append(line[7])
                remapGeneID.append(line[8])
                remapGeneID.append(ID)
                geneIDRemap.append(remapGeneID)

geneIDRemapFile = gene + ".all5DatasetsSpanningRegion.remap.gff.gz"
outputPath8 = geneOutputDir + "/" + geneIDRemapFile
geneIDRemapOutput = gzip.open(outputPath8, "w")

for line in geneIDRemap:
    for element in line:
        geneIDRemapOutput.write(str(element).encode())
        geneIDRemapOutput.write("\t".encode())
    geneIDRemapOutput.write("\n".encode())

geneIDRemapOutput.close()


# Ordering remapped geneIDMatch by isoform
print("Step 11 of 15 - Ordering remapped reads by isoform.")
strucRemap = []
intStRemap = []
infoRemap = ""

# Gather information for all reads with the same read ID which were sorted into
#     geneIDRemap in Step 9
for line in geneIDRemap:
    ID = line[9]

    # intSt will be empty to begin so first read will go to else function. If this
    #     ID already exists in intSt, then the info variable will add the read
    #     information associated with this ID
    if ID in intStRemap:
        end = int(line[3]) - 1
        infoRemap = infoRemap + ";%;" + line[0] + "_" + str(start) + "_" + str(end) + "_" + line[6]
        start = int(line[4]) + 1
    # if ID not found in intSt and intSt is not empty, append info to it, append
    #     that to struc then empty both intSt and info to restart for new read ID
    else:
        if len(intStRemap) != 0:
            intStRemap.append(infoRemap)
            strucRemap.append(intStRemap)
            intStRemap = []
            infoRemap = ""
        intStRemap.append(ID)
        start = int(line[4]) + 1
# Ensure last read ID is also appended to final array
intStRemap.append(infoRemap)
strucRemap.append(intStRemap)

def myFunc(e):
    return e[1]
strucRemap.sort(key=myFunc)

# Create output file path
outputOrderRemapFile = gene + ".order.remap.tab.gz"
outputPath9 = geneOutputDir + "/" + outputOrderRemapFile
orderRemapOutput = gzip.open(outputPath9, "w")

# Write struc to output file
for line in strucRemap:
    for element in line:
        orderRemapOutput.write(str(element).encode())
        orderRemapOutput.write("\t".encode())
    orderRemapOutput.write("\n".encode())

orderRemapOutput.close()


# Remap annotation to new coordinates
print("Step 12 of 15 - Remapping annotation to the new coordinates")

annoRemap = []
for line in annoContents:
    line[9] = line[9].split("\"")
    ID = line[9][1]
    line[11] = line[11].split("\"")
    transID = line[11][1]
    for i in remap:
        # Only worried about reads that fall within an exonic section of the gene
        if "exonic" in i[3]:
	    # Do not worry about reads that do not fall within these exonic sections
            if int(line[3]) < int(i[1]) or int(line[4]) > int(i[2]):
                continue
	    # If they do, recalculate start and end of reads and append to annoRemap
            else:
                remapAnno = []
                remapAnno.append(line[0])
                remapAnno.append(line[1])
                remapAnno.append(line[2])
                x = int(line[3]) - int(i[1])
                y = int(i[2]) - int(line[4])
                staDis = int(i[4]) + x
                endDis = int(i[5]) - y
                remapAnno.append(staDis)
                remapAnno.append(endDis)
                remapAnno.append(line[5])
                remapAnno.append(line[6])
                remapAnno.append(line[7])
                remapAnno.append(line[8])
                remapAnno.append(ID)
                remapAnno.append(line[10])
                remapAnno.append(transID)
                remapAnno.append("gene_name")
                remapAnno.append(gene)
                annoRemap.append(remapAnno)

# Sort list of lists by the transcript ID to ensure the logic of the next step
# works correctly
annoRemap = sorted(annoRemap, key=lambda x:x[11])

# Now we look for all transcript IDs and assign remapped coordinates to the
#    first and last base associated with that transcript ID
trans = [] 
transRead = [None] * len(annoRemap[0])
minimum = 10000000
maximum = 0
for read in annoRemap:
    for i in range(0,len(read)-1):
        if str(read[i]) == "transcript_id":
            transID = read[i+1]
    # Array starts as a list of "None" values. Ensure that you are working with a transID
    if transRead[11] != None:
	# If you are working with a transID, but have completed going through the associated reads
	#    and are now seeing a read with a different transID, append what you have for the last
	#    transID and reset values of transRead back to "None"
        if transID not in transRead[11]:
            trans.append(transRead)
            transRead = [None] * len(annoRemap[0])
            minimum = 10000000
            maximum = 0
        transRead[0] = read[0]
        transRead[1] = read[1]
        transRead[2] = "transcript"
	# Assign minimum and maximum as you come across them for a specific transID
        minimum = min(minimum, int(read[3]))
        maximum = max(maximum, int(read[4]))
        transRead[3] = minimum
        transRead[4] = maximum
        transRead[5] = read[5]
        transRead[6] = read[6]
        transRead[7] = "."
        transRead[8] = read[8]
        transRead[9] = read[9]
        transRead[10] = read[10]
        transRead[11] = transID
        transRead[12] = read[12]
        transRead[13] = read[13]
    # This should only run on the first pass through the for loop
    else:
        transRead[0] = read[0]
        transRead[1] = read[1]
        transRead[2] = "transcript"
        minimum = min(minimum, int(read[3]))
        maximum = max(maximum, int(read[4]))
        transRead[3] = minimum
        transRead[4] = maximum
        transRead[5] = read[5]
        transRead[6] = read[6]
        transRead[7] = "."
        transRead[8] = read[8]
        transRead[9] = read[9]
        transRead[10] = read[10]
        transRead[11] = transID
        transRead[12] = read[12]
        transRead[13] = read[13]

trans.append(transRead)

# Put the transcript ID information at the top of the annoRemap array
for line in trans:
    annoRemap.insert(0, line)

# Create output path
remapAnnoFile = gene + ".anno_remap.gtf.gz"
outputPath10 = geneOutputDir + "/" + remapAnnoFile
remapAnnoOutput = gzip.open(outputPath10, "w")

# Write to output file
for line in annoRemap:
    for element in line:
        remapAnnoOutput.write(str(element).encode())
        remapAnnoOutput.write("\t".encode())
    remapAnnoOutput.write("\n".encode())

remapAnnoOutput.close()


# Set variables for the plot
print("Step 13 of 15 - Setting variables for the plot")

cellType = open(cellTypeFile, "rb")

# Format file for use in the R script
plotData = []
for line in cellType:
    data = []
    line = line.decode()
    line = line.split()
    data.append(line[0])
    # Make sure for the file paths that the names with "/" are changed to "_"
    if "/" in line[0]:
        line[0] = line[0].replace("/", "_")
    data.append(geneOutputDir + "/" + gene + "." + line[0] + ".remap.gtf.gz")
    data.append(line[0] + " isoforms")
    data.append(line[1])
    plotData.append(data)

# Create output path
plotNamesFile = gene + ".cellTypeFileWithFileNames.tab"
outputPath11 = geneOutputDir + "/" + plotNamesFile
plotNamesOutput = open(outputPath11, "w")

# Write to output file
for line in plotData:
    for element in line:
        plotNamesOutput.write(str(element))
        plotNamesOutput.write("\t")
    plotNamesOutput.write("\n")

plotNamesOutput.close()


# Identify and remap SNVs, insertions, and deletions if option is chosen by user
if mismatchFile != None:
    print("Step 14 of 15 - Identifying and remapping SNVs, insertions, and deletions")
    
    # Read in mismatches file
    mismatches = [read.strip('\n') for read in gzip.open(mismatchFile, 'rt', encoding = 'utf-8').readlines()]
    
    # Create dictionaries for each type of mismatch
    SNVDict = {}
    insertDict = {} 
    delDict = {}
    # Separate out data by mismatch type
    for read in range(1,len(mismatches)):
        infoSNV = mismatches[read].split("\t")
        SNVDict[infoSNV[1]] = infoSNV[2]
        insertDict[infoSNV[1]] = infoSNV[3]
        delDict[infoSNV[1]] = infoSNV[4]
    
    # Organize SNVs by readID in dictionary
    SNVDictList = {}
    SNVKeys = SNVDict.keys()
    for readID in SNVKeys:
        # If SNV entry exists for the readID, then record each instance in the
        # dictionary
        if SNVDict[readID] != "NA":
            allSNVs = []
            readSNV = SNVDict[readID].split(';')
            for i in readSNV:
                SNV = i.split("_")
                allSNVs.append(SNV[0])
            SNVDictList[readID] = allSNVs
    
    # Organize insertions by readID in dictionary        
    insertDictList = {}
    insertKeys = insertDict.keys()
    for insertID in SNVKeys:
        # If SNV entry exists for the readID, then record each instance in the
        # dictionary
        if insertDict[insertID] != "NA":
            allInserts = []
            readInsert = insertDict[insertID].split(';')
            for j in readInsert:
                insert = j.split("_")
                allInserts.append(insert[0])
            insertDictList[insertID] = allInserts
    
    # Organize deletions by readID in dictionary            
    delDictList = {}
    delKeys = delDict.keys()
    for delID in delKeys:
        # If SNV entry exists for the readID, then record each instance in the
        # dictionary
        if delDict[delID] != "NA":
            allDels = []
            readDel = delDict[delID].split(';')
            for k in readDel:
                deletion = k.split("_")
                allDels.append(deletion[0])
            delDictList[delID] = allDels
    
    # Remap SNVs          
    SNVList = {}
    keysSNV = SNVDictList.keys()
    for key in keysSNV:
        snvRemap = []
        for i in range(0, len(SNVDictList[key])):
            # Comparing to the remapped projections table
            for line in remap:
                if "exonic" in line[3]:
                    # Use exonic start and end to remap SNV positions
                    if int(SNVDictList[key][i]) >= int(line[1]) and int(SNVDictList[key][i]) <= int(line[2]):
                        diff = int(line[2]) - int(SNVDictList[key][i])
                        remapSNV = int(line[5]) - diff
                        snvRemap.append(remapSNV)
        if len(snvRemap) > 0:
            SNVList[key] = snvRemap
    
    # Remap insertions
    insertList = {}
    keysInsert = insertDictList.keys()
    for key in keysInsert:
        insertRemap = []
        for i in range(0, len(insertDictList[key])):
            # Comparing to the remapped projections table
            for line in remap:
                if "exonic" in line[3]:
                    # Use exonic start and end to remap insertion positions
                    if int(insertDictList[key][i]) >= int(line[1]) and int(insertDictList[key][i]) <= int(line[2]):
                        diff = int(line[2]) - int(insertDictList[key][i])
                        remapInsert = int(line[5]) - diff
                        insertRemap.append(remapInsert)
        if len(insertRemap) > 0:
            insertList[key] = insertRemap
    
    #Remap deletions
    delList = {}
    keysDelete = delDictList.keys()
    for key in keysDelete:
        deleteRemap = []
        for i in range(0, len(delDictList[key])):
            # # Comparing to the remapped projections table
            for line in remap:
                if "exonic" in line[3]:
                    # Use exonic start and end to remap deletion positions
                    if int(delDictList[key][i]) >= int(line[1]) and int(delDictList[key][i]) <= int(line[2]):
                        diff = int(line[2]) - int(delDictList[key][i])
                        remapDelete = int(line[5]) - diff
                        deleteRemap.append(remapDelete)
        if len(deleteRemap) > 0:
            delList[key] = deleteRemap
    
    # Create SNVs file path        
    snvFile = gene + ".SNVs.tab"
    outputPath12 = geneOutputDir + "/" + snvFile
    snvOutput = open(outputPath12, "w")
    
    # Write SNVs to file by readID in dictionary
    SNVk = SNVList.keys()
    for key in SNVk:
        snvOutput.write(str(key))
        snvOutput.write("\t")
        for i in range(0, len(SNVList[key])):
            if i > 0:
                snvOutput.write(",")
            snvOutput.write(str(SNVList[key][i]))
        snvOutput.write("\n")
    snvOutput.close()
    
    # Create insertions file path        
    insertFile = gene + ".insertions.tab"
    outputPath13 = geneOutputDir + "/" + insertFile
    insertOutput = open(outputPath13, "w")
    
    # Write insertions to file by readID in dictionary
    insertK = insertList.keys()
    for key in insertK:
        insertOutput.write(str(key))
        insertOutput.write("\t")
        for i in range(0, len(insertList[key])):
            if i > 0:
                insertOutput.write(",")
            insertOutput.write(str(insertList[key][i]))
        insertOutput.write("\n")
    insertOutput.close()
    
    # Create deletions file path        
    deleteFile = gene + ".deletions.tab"
    outputPath14 = geneOutputDir + "/" + deleteFile
    deleteOutput = open(outputPath14, "w")
    
    # Write deletions to file by readID in dictionary
    delK = delList.keys()
    for key in delK:
        deleteOutput.write(str(key))
        deleteOutput.write("\t")
        for i in range(0, len(delList[key])):
            if i > 0:
                deleteOutput.write(",")
            deleteOutput.write(str(delList[key][i]))
        deleteOutput.write("\n")
    deleteOutput.close()


# Create UCSC upload file
print("Step 15 of 15 - Creating UCSC reference file for specified cell types.")

cellType = open(cellTypeFile, "rb")

ucscRefFile = gene + ".ucscReference.gtf.gz"
outputPath15 = geneOutputDir + "/" + ucscRefFile
ucscRefOutput = gzip.open(outputPath15, "w")

for cType in cellType:
    cType = cType.decode()
    cType = cType.split()
    
    readsPerCType = []
    
    color = list(np.random.choice(range(256), size=3))

    
    ucscTrackInfo = "track name=" + "\'" + str(cType[0]) + "\' " + "description=" + "\'" + str(cType[0]) + "\'" + " color=" + str(color[0]) + "," + str(color[1]) + "," + str(color[2]) + "\n"
    ucscRefOutput.write(str(ucscTrackInfo).encode())

    # Search for cellType matches throughout geneIDMatch array
    for line in geneIDMatch:
        ID = line[9]
        if cType[0] in ID:
            for i in remap:
		        # Only interested in exonic sections
                if "exonic" in i[3]:
		            # If read does not fall within the start and end of an exonic section, move on
                    if int(line[3]) < int(i[1]) or int(line[4]) > int(i[2]):
                        continue
                    else:
			            # If read falls within exonic section, recalc start and end point, and append to new readRemap
                        ucscRead = []
                        ucscRead.append(line[0])
                        ucscRead.append(line[1])
                        ucscRead.append(line[2])
                        ucscRead.append(line[3])
                        ucscRead.append(line[4])
                        ucscRead.append(line[5])
                        ucscRead.append(line[6])
                        ucscRead.append(line[7])
                        ucscRead.append("read_id " + "\"" + ID + "\";")
                        #ucscRead.append(ID)
                        readsPerCType.append(ucscRead)
    
    # Write to output file
    for line in readsPerCType:
        line = '\t'.join(line) + '\n'
        ucscRefOutput.write(line.encode())
        #for element in line:
        #    ucscRefOutput.write(str(element).encode())
        #    ucscRefOutput.write("\t".encode())
        #ucscRefOutput.write("\n".encode())

cellType.close()
ucscRefOutput.close()
