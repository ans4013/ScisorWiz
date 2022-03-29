#!/usr/bin/env python

###############################################################################
# Visualize data read-by-read; option for different cluster methods, option to
# plot genetic mismatches
#   - Input: interactive option, output file name, remapped anno file,
#            cellTypeFileWithFileNames.tab, order.tab.gz, all5 file,
#            altExons file, projection remap file, gene name,
#            cluster method, altExon inclusion cutoff (ci),
#            mismatch inclusion cutoff (optional), plot output path,
#            SNV file (optional), insertion file (optional),
#            deletion file (optional)
#   - Output: Plot of all isoforms of gene across all celltype choices
#
# Command: ./PlotIsoforms.r interactive option, plot name, anno remap file,
#           cellTypeFilewithFileNames, order.tab.gz, all5 File, altExons file,
#           projection remap file, gene name, cluster method, ci,
#           mismatch inclusion cutoff, plotOutput, SNVFile, insertionsFile,
#           deletionsFile
#
# by Hagen U. Tilgner (2012); modified and updated by Alexander N. Stein (2021)
################################################################################

if(!require(hash)) install.packages('hash')
if(!require(stats)) install.packages('stats')
if(!require(dplyr)) install.packages('dplyr')

# Determine clustering methods to use
getClustering<-function(a, m, cluster, num) {
    
    # If clustering method is specified do this
    if(length(a) > 9 & cluster != "Intron") {
        if(cluster == "TSS") {
            a <- a[order(a[,10]),]
        }
        if(cluster == "PolyA") {
            a <- a[order(a[,13]),]
        }
        if(cluster == "all") {
            a <- a[order(a[,10], a[,13], a[,14]),]
        }
        
        # Get only reads that match readID. If number of reads is >75, randomly
        # choose 75 reads to plot
        readOrder=rownames(m)[which(rownames(m) %in% a[,9])]
        if(length(readOrder) > 75) {
            readSample = sample(readOrder, 75)
            readOrder = readSample[order(match(readSample,a[,9]))]
        }
        else {
            readOrder = readOrder[order(match(readOrder, a[,9]))]
        }
    }
    
    # If unspecified clustering method or "intron chain" clustering method do
    # this.
    else {
        readOrder=rownames(m)[which(rownames(m) %in% a[,9])]
        if(length(readOrder) > 75) {
            readSample = sample(readOrder, 75)
            readOrder = readSample[order(match(readSample,readOrder))]
        }
    }
    printOut <- paste0("length(order", num ,"): ", length(readOrder))
    cat(printOut, "\n")
    return(readOrder)
}

# Function for getting mismatches specific to the reads of each cell type
getMismatches<-function(alignment, SNVs, insertions, deletions){
    reads2Mismatches = hash()
    
    # Enter removePath function to remove ".path1" if it exists in the readID
    alignment <- removePath(get(alignment))
    
    # Get all mismatches for each read for the cell type
    for(i in 1:length(alignment[,9])){
        readID = as.list(strsplit(as.character(alignment[i,9]), ":"))[[1]]
        x = dplyr::filter(SNVs, grepl(as.character(readID[3]), V1))
        y = dplyr::filter(insertions, grepl(as.character(readID[3]), V1))
        z = dplyr::filter(deletions, grepl(as.character(readID[3]), V1))
        mismatches = list(as.character(x[,2]), as.character(y[,2]), as.character(z[,2]))
        reads2Mismatches[as.character(alignment[i,9])] = mismatches
    }
    return(reads2Mismatches)
}

# Function to remove ".path1" from the readID strings
removePath<-function(alignments){
    alignments[,9]=as.vector(alignments[,9])
    
    # Iterate through all readIDs in table
    for(i in 1:length(alignments[,9])){
        readAlign <- unlist(strsplit(as.character(alignments[i,9]), ".path1"))
        alignments[i,9] = as.vector(readAlign)
    }
    return(alignments)
}

# Count number of mismatches for SNVs, insertion, and deletions to determine
# whether or not this mismatch is statistically relevant. This function removes
# random sequencing errors or outliers
countMismatches<-function(mismatch, totalReads, mismatchCutoff){
    
    # Separate out SNVs, insertions, and deletions for individual analysis
    SNVtmp <- as.list(values(mismatch)[1,])
    insTmp <- as.list(values(mismatch)[2,])
    delTmp <- as.list(values(mismatch)[3,])
    mismatchKeys <- as.list(keys(mismatch))

    # Split strings with SNV locations and add them to a list
    SNVList = list()
    for(i in 1:length(SNVtmp)){
        snvs <- as.list(strsplit(as.character(SNVtmp[i]), ","))
        for(i in 1:length(snvs[[1]])){
            suppressWarnings(SNVList <- append(SNVList, as.numeric(snvs[[1]][i])))
        }
    }

    # Count instances of a specific SNV location and determine if it is
    # statistically relevant. Tabulate gives tally for every position 0 to the
    # largest number in the list
    snvCount <- tabulate(as.numeric(unlist(SNVList)))
    snvIncl = list()
    for(i in 1:length(snvCount)){
        
        # Only enter if one or more instances of SNV at this position
        if(snvCount[[i]][1] != 0){
            
            # Determine inclusion rate for each position
            incl <- snvCount[[i]][1]/totalReads
            
            # Determine statistical significance of SNV and add it to list of
            # statistically significant SNV positions
            if(snvCount[[i]][1] != 0 & incl > mismatchCutoff & incl < (1-mismatchCutoff)) {
                snvIncl <- append(snvIncl, i)
            }
        }
    }

    # Split strings with insertion locations and add them to a list
    insList = list()
    for(i in 1:length(insTmp)){
        ins <- as.list(strsplit(as.character(insTmp[i]), ","))
        for(i in 1:length(ins[[1]])){
            suppressWarnings(insList <- append(insList, as.numeric(ins[[1]][i])))
        }
    }

    # Count instances of a specific insertion location and determine if it is
    # statistically relevant. Tabulate gives tally for every position 0 to the
    # largest number in the list
    insCount <- tabulate(as.numeric(unlist(insList)))
    insIncl = list()
    for(i in 1:length(insCount)){
        
        # Only enter if one or more instances of insertion at this position
        if(insCount[[i]][1] != 0){
            
            # Determine inclusion rate for each position
            incl <- insCount[[i]][1]/totalReads
            
            # Determine statistical significance of insertion and add it to list
            # of statistically significant SNV positions
            if(insCount[[i]][1] != 0 & incl > mismatchCutoff & incl < (1-mismatchCutoff)) {
                insIncl <- append(insIncl, i)
            }
        }
    }

    # Split strings with deletion locations and add them to a list
    delList = list()
    for(i in 1:length(delTmp)){
        del <- as.list(strsplit(as.character(delTmp[i]), ","))
        for(i in 1:length(del[[1]])){
            suppressWarnings(delList <- append(delList, as.numeric(del[[1]][i])))
        }
    }

    # Count instances of a specific deletion location and determine if it is
    # statistically relevant. Tabulate gives tally for every position 0 to the
    # largest number in the list
    delCount <- tabulate(as.numeric(unlist(delList)))
    delIncl = list()
    for(i in 1:length(delCount)){
        
        # Only enter if one or more instances of deletion at this position
        if(delCount[[i]][1] != 0){
            
            # Determine inclusion rate for each position
            incl <- delCount[[i]][1]/totalReads
            
            # Determine statistical significance of deletion and add it to list
            # of statistically significant SNV positions
            if(delCount[[i]][1] != 0 & incl > mismatchCutoff & incl < (1-mismatchCutoff)) {
                delIncl <- append(delIncl, i)
            }
        }
    }

    # Add all statistically significant positions for each type of mismatch to
    # a returnable list of lists with each type as its own list element
    includedMis <- list()
    if(length(snvIncl) > 0){
        includedMis <- append(includedMis, list(snvIncl))
    }
    
    # If no statistically relevant mismatches are present, NONE should take its
    # place in the list of lists
    else{
        includedMis <- append(includedMis, list("NONE"))
    }
    if(length(insIncl) > 0){
        includedMis <- append(includedMis, list(insIncl))
    }
    else{
        includedMis <- append(includedMis, list("NONE"))
    }
    if(length(delIncl) > 0){
        includedMis <- append(includedMis, list(delIncl))
    }
    else{
        includedMis <- append(includedMis, list("NONE"))
    }
    return(includedMis)
}

#getWindowStart<-function(windowStart, ){
#
#}

# Function to plot each read from a specific cell type
plotReads<-function(alignments,startLineNumber,numAlignedReads,alignmentHeader,
                    alignedReads,from,to,center,cexT,colReads,localExonRadius,
                    mismatchPos,tar1,tar2,totalCount,CI,altExons,specialColor){
    
    # Create header for cell type
    alignmentHeader = gsub("_", " ", alignmentHeader)
    
    # Add space before header
    lineNumber=startLineNumber-15
    text(center,lineNumber,alignmentHeader,col=colReads,cex=1*cexT)
    
    # Add space after header
    lineNumber=lineNumber-15;

    # If readID ends with ".path1", remove this substring
    if(grepl(".path1", as.character(alignments[1,9]), fixed=TRUE)){
        alignments = removePath(alignments)
    }

    # If user wants mismatches plotted, enter. Otherwise skip this code.
    if(!(is.empty(mismatchPos))){
        
        # Create lists
        incl <- list()
        snvIncl <- list()
        insIncl <- list()
        delIncl <- list()

        # Get the count of all mismatches. Only statistically significant
        # mismatches will be included in incl
        incl <- countMismatches(mismatchPos, totalCount, as.double(mismatchCutoff))
        
        # Separate out each type of mismatch
        snvIncl <- incl[[1]]
        insIncl <- incl[[2]]
        delIncl <- incl[[3]]
    }

    # Plot each read alignment
    for(readAlignment in alignedReads){
        
        # If read alignment has ".path1", remove it. Otherwise, do nothing and
        # set readAlign to current string.
        if(grepl(".path1", as.character(readAlignment), fixed=TRUE)){
            readAlign <- unlist(strsplit(as.character(readAlignment), ".path1"))
        }
        else{
            readAlign <- readAlignment
        }
        
        # Filter for only alignments which match this readID
        exons=alignments[which(alignments[,9]==readAlign),];
        
        # Get start and end of read
        start=min(exons[,4]);
        end=max(exons[,5]);

        # If alternate exons exist, specify and plot them
        for(i in 1:length(exons[,1])){
            
            # Test all exons in the read being plotted for alternate status
            t=paste(exons[i,1],exons[i,4],exons[i,5],exons[i,7],sep="_");
            
            # If match found, plot this exon as orange
            if(t %in% altExons$V1) {
                chosenColor=specialColor;
            }
            
            # If no match, plot as user chosen color
            else{
                chosenColor=colReads;
            }
            
            # Plot the exon as colored rectangle
            rect(exons[i,4],lineNumber-localExonRadius,exons[i,5],lineNumber+localExonRadius,col=chosenColor,border=chosenColor)
        }

        # If user wants mismatches plotted, enter this step.
        if(!(is.empty(mismatchPos))){
            
            # If this read alignment has mismatches, enter this step.
            if(length(mismatchPos[[readAlign]]) > 0){
                
                # Get positions of mismatches for this read alignment
                mismatches=values(mismatchPos, readAlign)
                
                # Separate out SNVs
                snvPos = as.list(strsplit(as.character(mismatches[[1]][1]), ","))[[1]]
                
                # Separate out insertions
                insertPos = as.list(strsplit(as.character(mismatches[[2]][1]), ","))[[1]]
                
                # Separate out deletions
                delPos = as.list(strsplit(as.character(mismatches[[3]][1]), ","))[[1]]
                
                # For each SNV, check to see if it is statistically significant
                # by comparing it to the list of statistically significant SNV
                # positions
                for(i in 1:length(snvPos)){
                    if(snvPos[i] %in% snvIncl){ # || as.integer(snvPos[i])-1 %in% snvIncl || snvPos[i]+1 %in% snvIncl){
                
                        # Remove any sequencing errors that can occur within the
                        # first 20 or last 20 bases and plot SNV as a cyan dot
                        if((as.integer(snvPos[i]) > start + 20) & (as.integer(snvPos[i]) < end - 20)) {
                            chosenCol="cyan";
                            rect(as.integer(snvPos[i])-10,lineNumber-localExonRadius,as.integer(snvPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                        }
                    }
                }
            
                # For each insertion, check to see if it is statistically
                # significant by comparing it to the list of statistically
                # significant insertion positions
                for(i in 1:length(insertPos)){
                    if(insertPos[i] %in% insIncl){
                        
                        # Remove any sequencing errors that can occur within the
                        # first 20 or last 20 bases and plot insertion as a
                        # green dot
                        if(insertPos[i] > start + 20 & insertPos[i] < end - 20){
                            chosenCol="chartreuse"
                            rect(as.integer(insertPos[i])-10,lineNumber-localExonRadius,as.integer(insertPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                        }
                    }
                }
                
                # For each deletion, check to see if it is statistically
                # significant by comparing it to the list of statistically
                # significant deletion positions
                for(i in 1:length(delPos)){
                    if(delPos[i] %in% delIncl){ # || delPos[i]-1 %in% delIncl || delPos[i]+1 %in% delIncl){
                        
                        # Remove any sequencing errors that can occur within the
                        # first 20 or last 20 bases and plot deletion as a red
                        # dot
                        if(delPos[i] > start + 20 & delPos[i] < end - 20){
                            chosenCol="red"
                            rect(as.integer(delPos[i])-10,lineNumber-localExonRadius,as.integer(delPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                        }
                    }
                }
            }
        }
        
        # Move line number to begin next read with down 2 for spacing purposes
        lineNumber=lineNumber-2;
    }
    return(lineNumber);
}

# Main plotting function
plotGenes<-function(annoGTF, aNames, annotationHeader, headerNames, colAnno, colReads,
                    cexT=1, SNVFile, insertFile, deleteFile, to, orderNames, windowStart,
                    windowEnd, drawAxis, altExons, projections, readLengths, CI,
                    mismatchCutoff, specialColor, numCT){
    
    cat("## starting plotGenes\n")
    anno=read.table(annoGTF)

    cat("## reading alignments\n")
    emptyDF=data.frame(factor(),factor(),factor(),integer(),integer(),factor(),factor(),factor(),factor())

    #for(i in 1:numCT){
    #    alignmentNames <- paste("alignments", i, sep='')
    #}
    #for(i in 1:numCT){
    #    assign(alignmentNames[i], tryCatch({res=read.table(get(alignedGFFs[i]), sep="\t",
    #                                                       colClasses = c("character", "character", "character", "integer", "integer", "character", "character", "character", "character", "integer", "integer", "integer", "integer", "character"))},
    #                                  warning = function(w) {cat("a warning was raised\n");
    #                                      return(emptyDF); },
    #                                  error = function(e) {return(emptyDF); }))
    #}
    #print(paste(alignmentNames[7], get(alignmentNames[7])))

    if(is.na(windowStart)){
        from = 100000
        plotEnd = 0
        
        # Find outer bounds of the plot
        for(i in 1:numCT ) {
            alignedMin = min(get(aNames[i])$V4)
            alignedMax = max(get(aNames[i])$V5)
            from = min(from, alignedMin)
            plotEnd = max(plotEnd, alignedMax)
        }
    }
    else{
        for(i in 1:nrow(projections)){
            if(i == windowStart){
                from = projections[i,5]
            }
            else if(i == windowEnd){
                plotEnd = projections[i,6]
            }
        }
    }

    # Add 100 to end of plot for spacing purposes
    plotEnd = plotEnd + 100
    
    # If the plot ends before the predetermined end, then that is the new plot
    # end for spacing purposes
    if(plotEnd < as.integer(to)) {
        to = plotEnd
    }

    # Extend start 100 bases left, if that number is less than 1, make it equal
    # to 1
    from = from - 100
    if(from < 1) {
        from = 1
    }
    cat("## the general plot\n");
    annoTransIndexes=which(anno[,3]=="transcript");
    numAnnoTrans=length(annoTransIndexes);
    y=anno;
    for(i in 1:length(annoTransIndexes)){
        x=anno[which(anno[,3]=="transcript" & anno[,12]==anno[i,12]),]
        if(x[1,4]>to | x[1,5]<from){
            y=y[which(y[,12]!=anno[i,12]),]
        }
    }
    anno=y
    annoTransIndexes=which(anno[,3]=="transcript");
    numAnnoTrans=length(annoTransIndexes);
    cat("numAnnoTrans =",numAnnoTrans,"\n");

    # Set width and center of plot
    showLength=to-from+1;
    center=(from+to)/2

    r2MNames <- paste("reads2Mismatches", 1:numCT, sep='')
    for(i in 1:numCT){
        assign(r2MNames[i], hash())
    }

    # If user chose to plot mismatches, enter this step.
    if(length(args) > 14){
        cat("## getting alignment IDs and numbers\n")
        
        # Create tables from the mismatch files
        SNVTable = read.table(SNVFile)
        insertTable = read.table(insertFile)
        deleteTable = read.table(deleteFile)

        # Get mismatches for each cell type
        cat("Getting SNVs, insertions, and deletions\n")
        for(i in 1:numCT){
            assign(r2MNames[i], getMismatches(aNames[i], SNVTable, insertTable, deleteTable))
        }
    }

    numAlignedReads <- paste("numAlignedReads", 1:numCT, sep='')
    for(i in 1:numCT){
        assign(numAlignedReads[i], length(get(orderNames[i])))
    }

    # Start plotting from top down
    cat("## the plot\n");
    ytop = 120 + (12*numAnnoTrans)
    for(i in 1:numCT){
        ytop = ytop + (2*get(numAlignedReads[i])) + 25
    }
    
    plot(c(from,to),c(1,ytop),col="white",axes=FALSE,xlab="",ylab="");

    xPositions=c((from-100):(to+100));
    xPositions=xPositions[which(xPositions %% 100 ==0)];

    # Plot lines
    for(i in 1:length(xPositions)){
        lines(c(xPositions[i],xPositions[i]),c(-1,ytop+2),col="lightcyan1")
    }

    axisPoints=c(from:to);
    axisPoints=axisPoints[which(axisPoints %% 2000 ==0)];
    # Plotting differently for different outputs. Default (hardcoded) for pdf
    if(interactive =="n"){
        axisOffset=50;
        localExonRadius=0.9;
    }
    else{
        axisOffset=110
        localExonRadius=0.15;
    }

    if(drawAxis==TRUE){
        for(i in c(1:length(axisPoints))){
            text(axisPoints[i]+axisOffset,ytop-3,ucscNumberVector(axisPoints[i]),cex=0.67*cexT,pos=2);
            lines(c(axisPoints[i],axisPoints[i]),c(ytop-6,ytop))
        }
        ytop=ytop-10;
    }
    else{
        ytop=ytop-5;
    }

    # Begin plotting cell types, adding space to prevent overlap of cell type
    # names and reads
    for(i in 1:numCT){
        cat("# the aligned reads",i, "\n");
        cat("enter with ",ytop,"; using ", get(numAlignedReads[i]),"reads\n")
        newLineNumber = plotReads(get(aNames[i]), ytop, get(numAlignedReads[i]),
                                  get(headerNames[i]), get(orderNames[i]), from, to,
                                  center, cexT, get(colReads[i]), localExonRadius,
                                  get(r2MNames[i]), target1, target2,
                                  readLengths[[i]], CI, altExons, specialColor);
        ytop=ytop - (2*get(numAlignedReads[i])) - 40
    }

    cat("# the annotated transcripts\n")
    cat("enter with ",ytop,";\n")

    lineNumber=ytop
    lineNumber=lineNumber-10-5
    if(interactive == "y"){
        extraShift=500;
    }
    else{
        extraShift=0;
    }

    # Plot annotations - starting with title of annotation section of plot
    text(center+extraShift,lineNumber,paste(annotationHeader,geneName[1],sep=" "),col=colAnno,cex=1*cexT);
    lineNumber=lineNumber-20

    # Plot transcript annotations, including only those exons which are labeled
    # as "exon" or "CDS"
    for(trIndex in annoTransIndexes){
        tr=anno[trIndex,12]
        gene=anno[trIndex,14]
        lines(c(anno[trIndex,4],anno[trIndex,5]),c(lineNumber,lineNumber),col=colAnno,lwd=0.2)

        exons=anno[which(anno[,3]=="exon" & anno[,12]==tr),];
        for(i in 1:length(exons[,1])){
            rect(exons[i,4],lineNumber-2,exons[i,5],lineNumber+2,col=colAnno,border=colAnno)
        }
        cds=anno[which(anno[,3]=="CDS" & anno[,12]==tr),];
        if(length(cds[,1])>0){
            for(i in 1:length(cds[,1])){
                rect(cds[i,4],lineNumber-4.7,cds[i,5],lineNumber+4.7,col=colAnno,border=colAnno);
            }
        }
        lineNumber=lineNumber-12;
    }

    # Plot exon numbers
    exonOrder = data.frame(start=integer(0), end=integer(0))
    
    # Ensure the numbers are located within exonic sections of gene
    for(proj in 1:nrow(projections)) {
        if(projections[proj,4] == "exonic") {
            exonOrder[nrow(exonOrder)+1,] = c(projections[proj,5], projections[proj,6])
        }
    }

    # Write in the actual numbers based on the center of the exon which with
    # they align
    for(i in 1:nrow(exonOrder)){
        loc = list(exonOrder$start[i], exonOrder$end[i])
        if(loc[[2]] > to) {
            textLoc = (loc[[1]] + to)/2
        }
        else {
            textLoc = (loc[[1]] + loc[[2]])/2
        }
        text(textLoc, lineNumber, i, col=colAnno, cex=.5*cexT)
    }
    lineNumber = lineNumber-12
}


# 0. reading parameters and data
args<-commandArgs(trailingOnly=TRUE);

interactive=args[1];
cat("interactive =", interactive,"\n");

outGraphFileName=args[2];
cat("outGraphFileName =", outGraphFileName,"\n");

exampleAnnotationGTF=args[3];
cat("exampleAnnotationGTF =", exampleAnnotationGTF,"\n");

cellTypeFileWithFileNames=args[4];
cat("cellTypeFileWithFileNames =", cellTypeFileWithFileNames,"\n");

cellTypeTableWithFileNames=read.table(cellTypeFileWithFileNames,sep="\t")

# Number of cell types listed
numCT <- nrow(cellTypeTableWithFileNames)

orderMatrixFile=args[5];
cat("orderMatrixFile =", orderMatrixFile,"\n");

alignedReadsFile=args[6];
cat("alignedReadsFile =", alignedReadsFile, "\n");

altExonsFile=args[7]
cat("altExonsFile =", altExonsFile, "\n")

projFile=args[8]
cat("projectionFile =", projFile, "\n")

geneName=args[9];
cat("geneName =", geneName,"\n");

cluster=args[10]
if(cluster == 1){
    cluster = "Intron"
}
if(cluster == 2){
    cluster = "TSS"
}
if(cluster == 3){
    cluster = "PolyA"
}
if(cluster == 4){
    cluster = "all"
}
cat("Cluster method =", cluster, "\n")

CI=args[11];
CI = format(round(as.integer(CI), 2), nsmall=2)
upper = 1 - as.double(CI)
upper = format(round(upper, 2), nsmall=2)
cat("Confidence Interval = ", as.double(CI), "-", as.double(upper),"\n");

mismatchCutoff=args[12]
cat("Mismatch Cutoff = ", mismatchCutoff, "\n")

output=args[13]
outputDir = paste0(output,geneName)
cat("outputDir = ", outputDir, "\n");

# If user wants to zoom in on a specific window, but doesn't want to plot
# mismatches, the length of arguments will be 16, enter this step
if(length(args) == 15){
    windowStart=args[14]
    windowEnd=args[15]
} else if(length(args) == 18){
    windowStart=args[17]
    windowEnd=args[18]
} else {
    windowStart=NA
    windowEnd=NA
}

# If user wants mismatches plotted, length of the arguments will be >=17, if
# there are more, then enter this step
if(length(args) >= 16){
    SNVFile=args[14]
    cat("SNVFile =", SNVFile, "\n")

    insertFile=args[15]
    cat("insertFile =", insertFile, "\n")

    deleteFile=args[16]
    cat("deleteFile =", deleteFile, "\n")

    # If user wants to zoom in on a specific window and want mismatches plotted,
    # the length of arguments will be > 17, enter this step
    #if(length(args) > 17){
    #    windowStart=args[18]
    #    windowEnd=args[19]
    #}
    #else{
    #    windowStart=NA
    #    windowEnd=NA
    #}
}

#
########################################
# 1. parameters

legendscale=1.0;

# 1h. names for the subfigures:
subfigureNames=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o");



########################################
# 2. the name of the plot

cat("# 2. the name of the output-file\n");

if(interactive == "y"){
    outGraphFileNameJpg=paste(outputDir,paste(outGraphFileName,"jpg",sep="."),sep="/");
    jpeg(file=outGraphFileNameJpg,width=8,height=numCT*2,unit="in",res=1200)
    cexMain=0.9;
    cexLab=1;
    legendscale=0.9
    marVector=c(2, 2, 2, 1) +0.1;
    mypadj=-3;
    myadj=2.5;
    omaVector=c(2,2,2,2)-1;
    mysubfigureLabelCex=1.2;
    mgpVector=c(1.2,0.5,0);
    mycex=2.5;
    cexText=0.65;
    mylwd=1.25;
    exampleY=0.49;
    box1X=c(-0.1,0.6);
    box1Y=c(0.0,9);
    box2X=c(8.3,9.7);
    box2Y=c(0.0,9);
    padjOffsetForFigure4=9.5
    xMtextLine=0.5;
} else {
    pdfWidth=11;
    pdfHeight=numCT*2;
    outGraphFileNamePdf=paste(outputDir,paste(outGraphFileName,"pdf",sep="."),sep = "/");
    pdf(file=outGraphFileNamePdf,width=pdfWidth,height=pdfHeight,title="all exons");
    cexMain=1.8;
    cexLab=2;
    legendscale=1.8;
    marVector=c(4, 4, 4, 2) +0.1;
    mypadj=-5.3;
    mycex=1.6;
    cexText=1.7;
    myadj=2;
    omaVector=c(1,1,1,3);
    mysubfigureLabelCex=2.5;
    mgpVector=c(2.4,1,0);
    mylwd=2.5;
    exampleY=0.495;
    box1X=c(0.4,1.1);
    box1Y=c(0.8,10);
    box2X=c(8.25,9.5);
    box2Y=c(0.8,10);
    xMtextLine=2.3;
}

layout(matrix(c(1,1,3,3,3,3,3,3,
                1,1,3,3,3,3,3,3,
                1,1,3,3,3,3,3,3,
                2,2,3,3,3,3,3,3,
                2,2,3,3,3,3,3,3,
                2,2,3,3,3,3,3,3),nrow=6,ncol=8,byrow=TRUE));


par(mar = marVector);
par(oma = omaVector);

########################################
# 3. parsing data and plotting:
cat("# 3. parsing data and plotting:\n");

plot(c(1:10),c(1:10),col="white",xlab="",ylab="",main="",axes=FALSE)
plot(c(1:10),c(1:10),col="white",xlab="",ylab="",main="",axes=FALSE)

cat("number of cell types:", numCT, "\n")

emptyDF=data.frame(factor(),factor(),factor(),integer(),integer(),factor(),factor(),factor(),factor())

# Get the alignments using names stored in the cellTypeTableWithFileNames
# file for each cell type.
if(interactive == "y"){
    par(mar = marVector-c(1.5,2,2,1))
} else{
    par(mar = marVector-c(3.5,4,4,2))
}

aNames <- paste('a',1:numCT,sep='')
for(i in 1:numCT){
    cat("reading", as.character(cellTypeTableWithFileNames[i,1]), "data\n")
    assign(aNames[i], tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[i,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }))
}

cat("reading in data\n")
orderMatrix=read.table(orderMatrixFile);
orderMatrixSub=orderMatrix
m=as.data.frame(orderMatrixSub[,2])
rownames(m)=orderMatrixSub[,1]
alignedReads=read.table(alignedReadsFile, sep = "\t");
altExons=read.table(altExonsFile)
projections=read.table(projFile)

# Get only sections which are exonic and find position of last base of exons
# in data
projected=projections[which(projections[,4]=="exonic"),]
rightPlotEnd = max(projected[,6])

# Cluster the data if user specified via AllInfo path. Default is 1 (intron
# chain).
orderNames <- paste('order', 1:numCT, sep='')
for(i in 1:numCT){
    temp <- getClustering(get(aNames[i]), m, cluster, i)
    assign(orderNames[i], temp)
}

# Get number of reads for total read count used in countMismatches function
readLengths <- list()
readOrderNames <- paste("readOrder", 1:numCT, sep='')
for(i in 1:numCT){
    assign(readOrderNames[i], rownames(m)[which(rownames(m) %in% get(aNames[i])[,9])])
    readLengths <- append(readLengths, length(get(readOrderNames[i])))
}

headerNames <- paste("annotationHeader", 1:numCT, sep='')
for(i in 1:numCT){
    assign(headerNames[i], as.character(cellTypeTableWithFileNames[i,3]))
}

colReadNames <- paste("colReads", 1:numCT, sep='')
for(i in 1:numCT){
    assign(colReadNames[i], as.character(cellTypeTableWithFileNames[i,4]))
}

# Call to plotting function
plotGenes(annoGTF=exampleAnnotationGTF, aNames, annotationHeader="Gencode annotation", headerNames,
          colAnno="black", colReadNames, cexT=mycex, SNVFile, insertFile, deleteFile, rightPlotEnd,
          orderNames, windowStart, windowEnd, drawAxis=FALSE, altExons=altExons, projections=projections, readLengths=readLengths,
          CI=CI, mismatchCutoff, specialColor="darkorange", numCT);
warnings()


dev.off();
