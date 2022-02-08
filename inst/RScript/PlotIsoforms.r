#!/usr/bin/env python

###############################################################################
# Visualize data read-by-read; option for different cluster methods, option to
# plot genetic mismatches
#   - Input: output file type, output file name, last figure,
#            remapped anno file, cellTypeFileWithFileNames.tab, order.tab.gz,
#            all5 file, altExons file, projection remap file, gene name,
#            cluster method, altExon inclusion cutoff (ci),
#            mismatch inclusion cutoff (optional), plot output path,
#            SNV file (optional), insertion file (optional),
#            deletion file (optional)
#   - Output: Plot of all isoforms of gene across all celltype choices
#
# Command: ./PlotIsoforms.r pdf, plot name, last figure, anno remap file,
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

# Not sure if this is really necessary
ucscNumber<-function(x){

    options(scipen=10)
    if(x<1000){
	    x=as.character(x);
	    return(x);
    }
    if(x<1000000){
	    x=as.character(x);
	    y=paste(substr(x,1,nchar(x)-3), substr(x,nchar(x)-2,nchar(x)),sep=",");
	    return(y);
    }

    if(x<1000000000){
	    x=as.character(x);
	    y=paste(substr(x,1,nchar(x)-6), substr(x,nchar(x)-5,nchar(x)-3),  substr(x,nchar(x)-2,nchar(x)),sep=",");
	    return(y);
    }
    stop("cant deal with numbers >=1 billion")
}

ucscNumberVector<-function(y){
    y=sapply(X=y,FUN=ucscNumber)
    return(y)
}

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
    alignment <- removePath(alignment)
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
            SNVList <- append(SNVList, as.numeric(snvs[[1]][i]))
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
                #print(paste("snv at location:", i, "value:", snvCount[[i]], "- inclusion:", incl))
                snvIncl <- append(snvIncl, i)
            }
        }
    }

    # Split strings with insertion locations and add them to a list
    insList = list()
    for(i in 1:length(insTmp)){
        ins <- as.list(strsplit(as.character(insTmp[i]), ","))
        for(i in 1:length(ins[[1]])){
            insList <- append(insList, as.numeric(ins[[1]][i]))
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
                #print(paste("insertion at location:", i, "value:", insCount[[i]], "- inclusion:", incl))
                insIncl <- append(insIncl, i)
            }
        }
    }

    # Split strings with deletion locations and add them to a list
    delList = list()
    for(i in 1:length(delTmp)){
        del <- as.list(strsplit(as.character(delTmp[i]), ","))
        for(i in 1:length(del[[1]])){
            delList <- append(delList, as.numeric(del[[1]][i]))
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
            #print(paste(delCount[[i]][1], totalReads, incl))
            #print(paste(delCount[[i]], totalReads, incl))
            # Determine statistical significance of deletion and add it to list
            # of statistically significant SNV positions
            if(delCount[[i]][1] != 0 & incl > mismatchCutoff & incl < (1-mismatchCutoff)) {
                #print(paste("deletion at location:", i, "value:", delCount[[i]], "- inclusion:", incl))
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
    #print("includedMis")
    #print(includedMis)
    return(includedMis)
}

# Function to plot each read from a specific cell type
plotReads<-function(alignments,startLineNumber,numAlignedReads,alignmentHeader,alignedReads,from,to,center,cexT,colReads,localExonRadius,mismatchPos,tar1,tar2,totalCount,CI,altExons,specialColor){
    # Create background
    #arrowDis=(max(alignments[,5])-min(alignments[,4]))/40;     #testing what would happen if this was removed
    #arrowPos=min(alignments[,4])+c(1:40)*arrowDis;     #testing what would happen if this was removed

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
        #print("Post countMismaches")
        # Separate out each type of mismatch
        snvIncl <- incl[[1]]
        insIncl <- incl[[2]]
        delIncl <- incl[[3]]

        #print(paste("snv:", snvIncl))
        #print(paste("ins:", insIncl))
        #print(paste("del:", delIncl))
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
        #arrowPosLoc=c();       #testing what would happen if this was removed

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
                #print("snvPos")
                #print(snvPos)
                # Separate out insertions
                insertPos = as.list(strsplit(as.character(mismatches[[2]][1]), ","))[[1]]
                #print("insertPos")
                #print(insertPos)
                # Separate out deletions
                delPos = as.list(strsplit(as.character(mismatches[[3]][1]), ","))[[1]]
                #print("delPos")
                #print(delPos)

                # For each SNV, check to see if it is statistically significant
                # by comparing it to the list of statistically significant SNV
                # positions
                for(i in 1:length(snvPos)){
                    if(snvPos[i] %in% snvIncl){ # || as.integer(snvPos[i])-1 %in% snvIncl || snvPos[i]+1 %in% snvIncl){
                        # Remove any sequencing errors that can occur within the
                        # first 20 or last 20 bases and plot SNV as a cyan dot
                        if((as.integer(snvPos[i]) > start + 20) & (as.integer(snvPos[i]) < end - 20)){
                            #print("valid snv found")
                            #print(snvPos[i])
                            chosenCol="cyan";
                            rect(as.integer(snvPos[i])-10,lineNumber-localExonRadius,as.integer(snvPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                        }
                    }
                }
                # For each insertion, check to see if it is statistically
                # significant by comparing it to the list of statistically
                # significant insertion positions
                for(i in 1:length(insertPos)){
                    if(insertPos[i] %in% insIncl){ # || insertPos[i]-1 %in% insIncl || insertPos[i]+1 %in% insIncl){
                        # Remove any sequencing errors that can occur within the
                        # first 20 or last 20 bases and plot insertion as a
                        # green dot
                        if(insertPos[i] > start + 20 & insertPos[i] < end - 20){
                            #print("valid ins found")
                            #print(insertPos[i])
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
                            #print("valid del found")
                            #print(delPos[i])
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
plotGenes<-function(annoGTF,alignedGFF1,alignedGFF2,alignedGFF3,alignedGFF4,alignedGFF5,alignedGFF6,annotationHeader,
                    alignmentHeader1,alignmentHeader2,alignmentHeader3,alignmentHeader4,alignmentHeader5,alignmentHeader6,
                    colAnno="darkgray",colReads1="steelblue4",colReads2="darkorange",colReads3="black",colReads4="black",
                    colReads5="black",colReads6="black",cexT=1,SNVFile,insertFile,deleteFile,to,alignedReads1,alignedReads2,alignedReads3,alignedReads4,
                    alignedReads5,alignedReads6,target1,target2,drawAxis,altExons,projections,readLengths,CI,mismatchCutoff,specialColor){
    cat("## starting plotGenes\n")
    anno=read.table(annoGTF)

    cat("## reading alignments\n")
    emptyDF=data.frame(factor(),factor(),factor(),integer(),integer(),factor(),factor(),factor(),factor())
    alignments1=tryCatch({res=read.table(alignedGFF1, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    alignments2=tryCatch({res=read.table(alignedGFF2, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    alignments3=tryCatch({res=read.table(alignedGFF3, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    alignments4=tryCatch({res=read.table(alignedGFF4, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    alignments5=tryCatch({res=read.table(alignedGFF5, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    alignments6=tryCatch({res=read.table(alignedGFF6, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })

    from = 100000
    plotEnd = 0
    # Find outer bounds of the plot
    for(i in 1:6) {
        alignment = paste0("alignments", i)
        alignments = eval(parse(text = alignment))
        alignedMin = min(alignments$V4)
        alignedMax = max(alignments$V5)
        from = min(from, alignedMin)
        plotEnd = max(plotEnd, alignedMax)
    }

    # Add 100 to end of plot for space purposes
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
    cat("numAnnoTrans=",numAnnoTrans,"\n");
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
    cat("numAnnoTrans=",numAnnoTrans,"\n");

    # Set width and center of plot
    showLength=to-from+1;
    center=(from+to)/2

    reads2Mismatches1 = hash()
    reads2Mismatches2 = hash()
    reads2Mismatches3 = hash()
    reads2Mismatches4 = hash()
    reads2Mismatches5 = hash()
    reads2Mismatches6 = hash()

    # If user chose to plot mismatches, enter this step.
    if(length(args) > 14){
        cat("## getting alignment IDs and numbers\n")
        # Create tables from the mismatch files
        SNVTable = read.table(SNVFile)
        insertTable = read.table(insertFile)
        deleteTable = read.table(deleteFile)

        # Get mismatches for each cell type
        cat("Getting SNVs, insertions, and deletions")
        reads2Mismatches1 = getMismatches(alignments1, SNVTable, insertTable, deleteTable)
        reads2Mismatches2 = getMismatches(alignments2, SNVTable, insertTable, deleteTable)
        reads2Mismatches3 = getMismatches(alignments3, SNVTable, insertTable, deleteTable)
        reads2Mismatches4 = getMismatches(alignments4, SNVTable, insertTable, deleteTable)
        reads2Mismatches5 = getMismatches(alignments5, SNVTable, insertTable, deleteTable)
        reads2Mismatches6 = getMismatches(alignments6, SNVTable, insertTable, deleteTable)
    }

    numAlignedReads1=length(alignedReads1);
    numAlignedReads2=length(alignedReads2);
    numAlignedReads3=length(alignedReads3);
    numAlignedReads4=length(alignedReads4);
    numAlignedReads5=length(alignedReads5);
    numAlignedReads6=length(alignedReads6);

    # Start plotting from top down
    cat("## the plot\n");
    ytop= 60 + 2*numAlignedReads1+25 + 2*numAlignedReads2+25 + 2*numAlignedReads3+25 + 2*numAlignedReads4+25 + 2*numAlignedReads5+25 + 2*numAlignedReads6+25  +  (12*numAnnoTrans) + 60
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
    if(eps=="eps" | eps=="pdf"){
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
    cat("# the aligned reads 1\n");
    cat("enter with ",ytop,"; using ", numAlignedReads1,"reads\n")
    newLineNumber=plotReads(alignments1,ytop,numAlignedReads1,alignmentHeader1,alignedReads1,from,to,center,cexT,colReads1,
                            localExonRadius,reads2Mismatches1,target1,target2,readLengths[[1]],CI,altExons,specialColor);

    cat("# the aligned reads 2\n");
    ytop=ytop- (2*numAlignedReads1)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads2,"reads\n")
    newLineNumber=plotReads(alignments2,ytop,numAlignedReads2,alignmentHeader2,alignedReads2,from,to,center,cexT,colReads2,
                            localExonRadius,reads2Mismatches2,target1,target2,readLengths[[2]],CI,altExons,specialColor);

    cat("# the aligned reads 3\n");
    ytop=ytop- (2*numAlignedReads2)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads3,"reads\n")
    newLineNumber=plotReads(alignments3,ytop,numAlignedReads3,alignmentHeader3,alignedReads3,from,to,center,cexT,colReads3,
                            localExonRadius,reads2Mismatches3,target1,target2,readLengths[[3]],CI,altExons,specialColor);

    cat("# the aligned reads 4\n");
    ytop=ytop- (2*numAlignedReads3)-40;
    cat("enter with ",ytop ,"; using ", numAlignedReads4,"reads\n")
    newLineNumber=plotReads(alignments4,ytop,numAlignedReads4,alignmentHeader4,alignedReads4,from,to,center,cexT,colReads4,
                            localExonRadius,reads2Mismatches4,target1,target2,readLengths[[4]],CI,altExons,specialColor);

    cat("# the aligned reads 5\n");
    ytop=ytop- (2*numAlignedReads4)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads5,"reads\n")
    newLineNumber=plotReads(alignments5,ytop,numAlignedReads5,alignmentHeader5,alignedReads5,from,to,center,cexT,colReads5,
                            localExonRadius,reads2Mismatches5,target1,target2,readLengths[[5]],CI,altExons,specialColor);

    cat("# the aligned reads 6\n");
    ytop=ytop- (2*numAlignedReads5)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads6,"reads\n")
    newLineNumber=plotReads(alignments6,ytop,numAlignedReads6,alignmentHeader6,alignedReads6,from,to,center,cexT,colReads6,
                            localExonRadius,reads2Mismatches6,target1,target2,readLengths[[6]],CI,altExons,specialColor);

    cat("# the annotated transcripts\n")
    ytop=ytop- (2*numAlignedReads6)-40;
    cat("enter with ",ytop,";\n")

    lineNumber=ytop
    lineNumber=lineNumber-10-5
    if(eps=="jpg"){
        extraShift=500;
    }
    else{
        extraShift=0;
    }

    # Plot annotations - starting with title of annotation section of plot
    text(center+extraShift,lineNumber,paste(annotationHeader,geneName[1],sep=" "),col=colAnno,cex=1*cexT);
    lineNumber=lineNumber-20
    #arrowDis=(max(anno[,5])-min(anno[,4]))/40;     #testing what would happen if this was removed
    #arrowPos=min(anno[,4])+c(1:40)*arrowDis;   #testing what would happen if this was removed

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

eps=args[1];
cat("eps =", eps,"\n");

outGraphFileName=args[2];
cat("outGraphFileName =", outGraphFileName,"\n");

lastFigure=as.numeric(args[3]);
cat("lastFigure =", lastFigure,"\n");

exampleAnnotationGTF=args[4];
cat("exampleAnnotationGTF =", exampleAnnotationGTF,"\n");

cellTypeFileWithFileNames=args[5];
cat("cellTypeFileWithFileNames =", cellTypeFileWithFileNames,"\n");

orderMatrixFile=args[6];
cat("orderMatrixFile =", orderMatrixFile,"\n");

alignedReadsFile=args[7];
cat("alignedReadsFile =", alignedReadsFile, "\n");

altExonsFile=args[8]
cat("altExonsFile =", altExonsFile, "\n")

projFile=args[9]
cat("projectionFile =", projFile, "\n")

geneName=args[10];
cat("geneName =", geneName,"\n");

cluster=args[11]
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

CI=args[12];
CI = format(round(as.integer(CI), 2), nsmall=2)
upper = 1 - as.double(CI)
upper = format(round(upper, 2), nsmall=2)
cat("Confidence Interval = ", as.double(CI), "-", as.double(upper),"\n");

mismatchCutoff=args[13]
cat("Mismatch Cutoff = ", mismatchCutoff, "\n")

output=args[14]
outputDir = paste0(output,geneName)
cat("outputDir = ", outputDir, "\n");

# If user wants mismatches plotted, length of the arguments will be <=14, if
# there are more, then enter this step
if(length(args) > 14){
    SNVFile=args[15]
    cat("SNVFile =", SNVFile, "\n")

    insertFile=args[16]
    cat("insertFile =", insertFile, "\n")

    deleteFile=args[17]
    cat("deleteFile =", deleteFile, "\n")
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

if(eps=="pdf"){
    pdfWidth=11;
    pdfHeight=8;
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


if(eps=="eps"){
    epsWidth=11.4;
    epsHeight=10;
    outGraphFileNameEps=paste(outputDir,paste(outGraphFileName,"eps",sep="."),sep="/");
    postscript(file=outGraphFileNameEps,width=epsWidth,height=epsHeight,title="",horizontal=T);
    cexMain=2.1;
    cexLab=2.1;
    legendscale=2.5;
    marVector=c(4, 4, 4, 2) +0.1;
    mypadj=-3.5;
    mycex=2;
    cexText=1.4;
    myadj=2;
    omaVector=c(1,1,1,3);
    mysubfigureLabelCex=2.5;
    mgpVector=c(2.6,0.9,0);
    mylwd=3.5;
    exampleY=0.495;
    box1X=c(0.4,1.1);
    box1Y=c(0.8,10);
    box2X=c(8.25,9.5);
    box2Y=c(0.8,10);
    xMtextLine=1.4;
    padjOffsetForFigure4=8.5;
}

if(eps=="png"){
    pngWidth=1080;
    pngHeight=820;
    outGraphFileNamePng=paste(outputDir,paste(outGraphFileName,"png",sep="."),sep="/");
    png(file=outGraphFileNamePng,width=pngWidth,height=pngHeight,units = "px");
    cexMain=2.2;
    cexLab=2.5;
    legendscale=2.5
    marVector=c(6, 6, 4, 2) +0.1;
    mypadj=-12.0;
    mycex=2.2;
    cexText=1.6;
    myadj=1.5;
    omaVector=c(3,2,3,3);
    mysubfigureLabelCex=3.2;
    mgpVector=c(3,1,0)
    mylwd=4.5;
    exampleY=0.495;
    box1X=c(0.7,1.4);
    box1Y=c(0.8,8.2);
    box2X=c(8.1,9.5);
    box2Y=c(0.8,8.2);
    padjOffsetForFigure4=9.5;
    xMtextLine=0.6;
}

if(eps=="jpg"){
    outGraphFileNameJpg=paste(outputDir,paste(outGraphFileName,"jpg",sep="."),sep="/");
    jpeg(file=outGraphFileNameJpg,width=5,height=4,unit="in",res=1200)
    cexMain=0.8;
    cexLab=0.8;
    legendscale=0.8
    marVector=c(2,2,2,1)+0.1;
    mypadj=-3;
    myadj=2.5;
    omaVector=c(2,2,2,2)-1;
    mysubfigureLabelCex=1.2;
    mgpVector=c(0.95,0.3,0);
    mycex=2.5;
    cexText=0.65;
    mylwd=1.6;
    exampleY=0.49;
    box1X=c(-0.1,0.6);
    box1Y=c(0.0,9);
    box2X=c(8.3,9.7);
    box2Y=c(0.0,9);
    padjOffsetForFigure4=9.5
    xMtextLine=0.5;
}

if(eps!="eps" && eps!="png" && eps!="jpg" && eps!="pdf"){
    stop(eps, " is not a suppoted file format\n")
}


layout(matrix(c(1,1,3,3,3,3,3,3,
                1,1,3,3,3,3,3,3,
                1,1,3,3,3,3,3,3,
                2,2,3,3,3,3,3,3,
                2,2,3,3,3,3,3,3,
                2,2,3,3,3,3,3,3,
                4,4,5,5,6,6,7,7),nrow=7,ncol=8,byrow=TRUE));



par(mar=marVector)
par(oma= omaVector);
########################################
# 3. parsing data and plotting:
cat("# 3. parsing data and plotting:\n");



if(lastFigure>=1){
    plot(c(1:10),c(1:10),col="white",xlab="",ylab="",main="",axes=FALSE)
    plot(c(1:10),c(1:10),col="white",xlab="",ylab="",main="",axes=FALSE)
}

if(lastFigure>=3){
    #cat("############################################## example\n");

    cellTypeTableWithFileNames=read.table(cellTypeFileWithFileNames,sep="\t")

    emptyDF=data.frame(factor(),factor(),factor(),integer(),integer(),factor(),factor(),factor(),factor())

    # Get the alignments using names stored in the cellTypeTableWithFileNames
    # file for each cell type.
    par(mar=marVector-c(3.5,4,4,2))
    cat("reading", as.character(cellTypeTableWithFileNames[1,1]), "data\n")
    a1=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[1,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    #print(dim(a1))

     cat("reading", as.character(cellTypeTableWithFileNames[2,1]), "data\n")
    a2=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[2,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    #print(dim(a2))

    cat("reading", as.character(cellTypeTableWithFileNames[3,1]), "data\n")
    a3=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[3,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    #print(dim(a3))

    cat("reading", as.character(cellTypeTableWithFileNames[4,1]), "data\n")
    a4=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[4,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    #print(dim(a4))

    cat("reading", as.character(cellTypeTableWithFileNames[5,1]), "data\n")
    a5=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[5,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); })
    #print(dim(a5))

    cat("reading", as.character(cellTypeTableWithFileNames[6,1]), "data\n")
    a6=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[6,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }) #, finally = {cat("try-catch done\n")})
    #print(dim(a6))

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
    order1 <- getClustering(a1, m, cluster, 1)
    order2 <- getClustering(a2, m, cluster, 2)
    order3 <- getClustering(a3, m, cluster, 3)
    order4 <- getClustering(a4, m, cluster, 4)
    order5 <- getClustering(a5, m, cluster, 5)
    order6 <- getClustering(a6, m, cluster, 6)

    # Get number of reads for total read count used in countMismatches function
    readLengths <- list()
    readOrder1=rownames(m)[which(rownames(m) %in% a1[,9])]
    readLengths <- append(readLengths, length(readOrder1))
    readOrder2=rownames(m)[which(rownames(m) %in% a2[,9])]
    readLengths <- append(readLengths, length(readOrder2))
    readOrder3=rownames(m)[which(rownames(m) %in% a3[,9])]
    readLengths <- append(readLengths, length(readOrder3))
    readOrder4=rownames(m)[which(rownames(m) %in% a4[,9])]
    readLengths <- append(readLengths, length(readOrder4))
    readOrder5=rownames(m)[which(rownames(m) %in% a5[,9])]
    readLengths <- append(readLengths, length(readOrder5))
    readOrder6=rownames(m)[which(rownames(m) %in% a6[,9])]
    readLengths <- append(readLengths, length(readOrder6))

    # Call to plotting function
    plotGenes(annoGTF=exampleAnnotationGTF,alignedGFF1=as.character(cellTypeTableWithFileNames[1,2]),
              alignedGFF2=as.character(cellTypeTableWithFileNames[2,2]),alignedGFF3=as.character(cellTypeTableWithFileNames[3,2]),
              alignedGFF4=as.character(cellTypeTableWithFileNames[4,2]),alignedGFF5=as.character(cellTypeTableWithFileNames[5,2]),
              alignedGFF6=as.character(cellTypeTableWithFileNames[6,2]),annotationHeader="Gencode annotation",
              alignmentHeader1=as.character(cellTypeTableWithFileNames[1,3]), alignmentHeader2=as.character(cellTypeTableWithFileNames[2,3]),
              alignmentHeader3=as.character(cellTypeTableWithFileNames[3,3]), alignmentHeader4=as.character(cellTypeTableWithFileNames[4,3]),
              alignmentHeader5=as.character(cellTypeTableWithFileNames[5,3]),alignmentHeader6=as.character(cellTypeTableWithFileNames[6,3]),
              colAnno="black",colReads1=as.character(cellTypeTableWithFileNames[1,4]),colReads2=as.character(cellTypeTableWithFileNames[2,4]),
              colReads3=as.character(cellTypeTableWithFileNames[3,4]),colReads4=as.character(cellTypeTableWithFileNames[4,4]),
              colReads5=as.character(cellTypeTableWithFileNames[5,4]),colReads6=as.character(cellTypeTableWithFileNames[6,4]),cexT=mycex,
              SNVFile,insertFile,deleteFile,rightPlotEnd,alignedReads1=order1,alignedReads2=order2,alignedReads3=order3,alignedReads4=order4,alignedReads5=order5,
              alignedReads6=order6,56889953,56890775,drawAxis=FALSE,altExons=altExons,projections=projections,readLengths=readLengths,
              CI=CI,mismatchCutoff,specialColor="darkorange");
    warnings()
}


dev.off();