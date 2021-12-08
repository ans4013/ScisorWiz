# by Hagen, January 1st, 2012
if(!require(hash)) install.packages('hash')
if(!require(stats)) install.packages('stats')
if(!require(dplyr)) install.packages('dplyr')


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
        readOrder=rownames(m)[which(rownames(m) %in% a[,9])]
        if(length(readOrder) > 75) {
            readSample = sample(readOrder, 75)
            readOrder = readSample[order(match(readSample,a[,9]))]
        }
        else {
            readOrder = readOrder[order(match(readOrder, a[,9]))]
        }
    }
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

getMismatches<-function(alignment, SNVs, insertions, deletions){
    reads2Mismatches = hash()
    print("getting SNVs, insertions, and deletions")
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

removePath<-function(alignments){
    for(i in length(alignments[,9])){
        readAlign <- as.list(strsplit(as.character(alignments[i,9]), ".path"))[[1]]
        print(readAlign)
        alignments[i,9] = readAlign[1]
        print(alignments[i,9])
    }
    return(alignments)
}

plotReads<-function(alignments,startLineNumber,numAlignedReads,alignmentHeader,alignedReads,center,cexT,colReads,localExonRadius,mismatchPos,tar1,tar2,totalCount,CI,altExons,specialColor){
    arrowDis=(max(alignments[,5])-min(alignments[,4]))/40;
    arrowPos=min(alignments[,4])+c(1:40)*arrowDis;

    alignmentHeader = gsub("_", " ", alignmentHeader)
    lineNumber=startLineNumber-15
    text(center,lineNumber,alignmentHeader,col=colReads,cex=1*cexT)
    lineNumber=lineNumber-15;

    for(readAlignment in alignedReads){
        exons=alignments[which(alignments[,9]==readAlignment),];
        start=min(exons[,4]);
        end=max(exons[,5]);
        arrowPosLoc=c();

        for(i in 1:length(exons[,1])){
            t=paste(exons[i,1],exons[i,4],exons[i,5],exons[i,7],sep="_");
            if(t %in% altExons$V1) {
                chosenColor=specialColor;
            }
            else{
                chosenColor=colReads;
            }
            rect(exons[i,4],lineNumber-localExonRadius,exons[i,5],lineNumber+localExonRadius,col=chosenColor,border=chosenColor)
        }

        # showing mismatches
        if(!(is.empty(mismatchPos))){
            #print(alignments[1,9])
            if(grepl(".path1", as.character(alignments[1,9]), fixed=TRUE)){
                print(alignments[1,9])
                removePath(alignments)
            }
            #print(alignments[1,9])
            exons=alignments[which(alignments[,9] == readAlignment),]
            if(length(mismatchPos[[readAlignment]]) > 0){
                #print(readAlignment)
                #print(mismatchPos[[readAlignment]])
                mismatches=values(mismatchPos, readAlignment)
                #print("HI")
                #print(mismatches)
                snvPos = as.list(strsplit(as.character(mismatches[1,]), ","))[[1]]
                insertPos = as.list(strsplit(as.character(mismatches[2,]), ","))[[1]]
                delPos = as.list(strsplit(as.character(mismatches[3,]), ","))[[1]]
                for(i in 1:length(snvPos)){
                    chosenCol="cyan";
                    rect(as.integer(snvPos[i])-10,lineNumber-localExonRadius,as.integer(snvPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                }
                for(i in 1:length(insertPos)){
                    chosenCol="chartreuse"
                    rect(as.integer(insertPos[i])-10,lineNumber-localExonRadius,as.integer(insertPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                }
                for(i in 1:length(delPos)){
                    chosenCol="red"
                    rect(as.integer(delPos[i])-10,lineNumber-localExonRadius,as.integer(delPos[i])+10,lineNumber+localExonRadius,col=chosenCol,border=chosenCol);
                }
            }
        }
        lineNumber=lineNumber-2;
    }
    return(lineNumber);
}


plotGenes<-function(annoGTF,alignedGFF1,alignedGFF2,alignedGFF3,alignedGFF4,alignedGFF5,alignedGFF6,annotationHeader,
                    alignmentHeader1,alignmentHeader2,alignmentHeader3,alignmentHeader4,alignmentHeader5,alignmentHeader6,
                    colAnno="darkgray",colReads1="steelblue4",colReads2="darkorange",colReads3="black",colReads4="black",
                    colReads5="black",colReads6="black",cexT=1,SNVFile,insertFile,deleteFile,to,alignedReads1,alignedReads2,alignedReads3,alignedReads4,
                    alignedReads5,alignedReads6,target1,target2,drawAxis,totalCount,altExons,projections,CI,specialColor){
    cat("## starting plotGenes\n")
    anno=read.table(annoGTF)

    totalCount <- nrow(orderMatrix)

    cat("## reading alignments\n")
    emptyDF=data.frame(factor(),factor(),factor(),integer(),integer(),factor(),factor(),factor(),factor())
    alignments1=tryCatch({res=read.table(alignedGFF1, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    alignments2=tryCatch({res=read.table(alignedGFF2, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    alignments3=tryCatch({res=read.table(alignedGFF3, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    alignments4=tryCatch({res=read.table(alignedGFF4, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    alignments5=tryCatch({res=read.table(alignedGFF5, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    alignments6=tryCatch({res=read.table(alignedGFF6, sep="\t")}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})

    from = 100000
    plotEnd = 0
    allAligned = data.frame()
    for(i in 1:6) {
        alignment = paste0("alignments", i)
        alignments = eval(parse(text = alignment))
        alignedMin = min(alignments$V4)
        alignedMax = max(alignments$V5)
        from = min(from, alignedMin)
        plotEnd = max(plotEnd, alignedMax)
    }

    plotEnd = plotEnd + 100
    if(plotEnd < as.integer(to)) {
        to = plotEnd
    }

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

    showLength=to-from+1;
    center=(from+to)/2

    reads2Mismatches1 = hash()
    reads2Mismatches2 = hash()
    reads2Mismatches3 = hash()
    reads2Mismatches4 = hash()
    reads2Mismatches5 = hash()
    reads2Mismatches6 = hash()

    print(length(args))
    if(length(args) > 13){
        cat("## getting alignment IDs and numbers\n")
        SNVTable = read.table(SNVFile)
        #print(SNVTable)
        insertTable = read.table(insertFile)
        deleteTable = read.table(deleteFile)

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

    cat("## the plot\n");
    ytop= 60 + 2*numAlignedReads1+25 + 2*numAlignedReads2+25 + 2*numAlignedReads3+25 + 2*numAlignedReads4+25 + 2*numAlignedReads5+25 + 2*numAlignedReads6+25  +  (12*numAnnoTrans) + 60
    plot(c(from,to),c(1,ytop),col="white",axes=FALSE,xlab="",ylab="");

    xPositions=c((from-100):(to+100));
    xPositions=xPositions[which(xPositions %% 100 ==0)];
    for(i in 1:length(xPositions)){
        lines(c(xPositions[i],xPositions[i]),c(-1,ytop+2),col="lightcyan1")
    }


    axisPoints=c(from:to);
    axisPoints=axisPoints[which(axisPoints %% 2000 ==0)];
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

    cat("# the aligned reads 1\n");
    cat("enter with ",ytop,"; using ", numAlignedReads1,"reads\n")
    newLineNumber=plotReads(alignments1,ytop,numAlignedReads1,alignmentHeader1,alignedReads1,center,cexT,colReads1,
                            localExonRadius,reads2Mismatches1,target1,target2,totalCount,CI,altExons,specialColor);

    cat("# the aligned reads 2\n");
    ytop=ytop- (2*numAlignedReads1)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads2,"reads\n")
    newLineNumber=plotReads(alignments2,ytop,numAlignedReads2,alignmentHeader2,alignedReads2,center,cexT,colReads2,
                            localExonRadius,reads2Mismatches2,target1,target2,totalCount,CI,altExons,specialColor);

    cat("# the aligned reads 3\n");
    ytop=ytop- (2*numAlignedReads2)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads3,"reads\n")
    newLineNumber=plotReads(alignments3,ytop,numAlignedReads3,alignmentHeader3,alignedReads3,center,cexT,colReads3,
                            localExonRadius,reads2Mismatches3,target1,target2,totalCount,CI,altExons,specialColor);

    cat("# the aligned reads 4\n");
    ytop=ytop- (2*numAlignedReads3)-40;
    cat("enter with ",ytop ,"; using ", numAlignedReads4,"reads\n")
    newLineNumber=plotReads(alignments4,ytop,numAlignedReads4,alignmentHeader4,alignedReads4,center,cexT,colReads4,
                            localExonRadius,reads2Mismatches4,target1,target2,totalCount,CI,altExons,specialColor);

    cat("# the aligned reads 5\n");
    ytop=ytop- (2*numAlignedReads4)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads5,"reads\n")
    newLineNumber=plotReads(alignments5,ytop,numAlignedReads5,alignmentHeader5,alignedReads5,center,cexT,colReads5,
                            localExonRadius,reads2Mismatches5,target1,target2,totalCount,CI,altExons,specialColor);

    cat("# the aligned reads 6\n");
    ytop=ytop- (2*numAlignedReads5)-40;
    cat("enter with ",ytop,"; using ", numAlignedReads6,"reads\n")
    newLineNumber=plotReads(alignments6,ytop,numAlignedReads6,alignmentHeader6,alignedReads6,center,cexT,colReads6,
                            localExonRadius,reads2Mismatches6,target1,target2,totalCount,CI,altExons,specialColor);

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

    text(center+extraShift,lineNumber,paste(annotationHeader,geneName[1],sep=" "),col=colAnno,cex=1*cexT);
    lineNumber=lineNumber-20
    arrowDis=(max(anno[,5])-min(anno[,4]))/40;
    arrowPos=min(anno[,4])+c(1:40)*arrowDis;
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

    exonOrder = data.frame(start=integer(0), end=integer(0))
    for(proj in 1:nrow(projections)) {
        if(projections[proj,4] == "exonic") {
            exonOrder[nrow(exonOrder)+1,] = c(projections[proj,5], projections[proj,6])
        }
    }
    #startEnd = data.frame(start=integer(0), end=integer(0))
    #for(trIndex in annoTransIndexes){
    #    tr=anno[trIndex,12]
    #    allExons=anno[which(anno[,3]=="CDS" | anno[,3] == "exon" & anno[,12]==tr),]
    #    print("allExons") #test
    #    print(length(allExons[1,])) #test
    #    for(i in 1:length(allExons)){
    #        len = nrow(startEnd)
    #        startEnd[nrow(startEnd)+1,] = c(allExons[i,4], allExons[i,5])
    #    }
    #}
    #print(startEnd)
    #print(startEnd[order(startEnd$start),])
    #startEnd <- unique(startEnd)
    #startEnd <- startEnd[order(startEnd$start),]

    #print(startEnd)
    #exonOrder = data.frame(start=integer(0), end=integer(0))
    #count = 1
    #for(exon in 1:nrow(startEnd)){
    #    if(exon == count){
    #        readStart = startEnd$start[exon]
    #        readEnd = startEnd$end[exon]
    #        data <- startEnd %>% filter(start == readStart | end == readEnd)
    #        exonOrder[nrow(exonOrder)+1,] = c(min(data$start), min(data$end))
    #        count = count + nrow(data)
    #    }
    #}
    #print(exonOrder)

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
upper = 1 - as.integer(CI)
upper = format(round(upper, 2), nsmall=2)
cat("Confidence Interval = ", CI, "-", upper,"\n");

output=args[13]
outputDir = paste0(output,geneName)
cat("outputDir = ", outputDir, "\n");

if(length(args) > 13){
    SNVFile=args[14]
    cat("SNVFile =", SNVFile, "\n")

    insertFile=args[15]
    cat("insertFile =", insertFile, "\n")

    deleteFile=args[16]
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
    outGraphFileNamePdf=paste(outputDir,paste(paste(outGraphFileName,lastFigure,sep="_"),"pdf",sep="."),sep = "/");
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
    outGraphFileNameEps=paste(outputDir,paste(paste(outGraphFileName,lastFigure,sep="_"),"eps",sep="."),sep="/");
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
    outGraphFileNamePng=paste(outputDir,paste(paste(outGraphFileName,lastFigure,sep="_"),"png",sep="."),sep="/");
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
    outGraphFileNameJpg=paste(outputDir,paste(paste(outGraphFileName,lastFigure,sep="_"),"jpg",sep="."),sep="/");
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
    cat("############################################## example\n");

    cellTypeTableWithFileNames=read.table(cellTypeFileWithFileNames,sep="\t")

    emptyDF=data.frame(factor(),factor(),factor(),integer(),integer(),factor(),factor(),factor(),factor())

    par(mar=marVector-c(3.5,4,4,2))
    cat("reading", as.character(cellTypeTableWithFileNames[1,1]), "data\n")
    a1=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[1,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    print(dim(a1))

     cat("reading", as.character(cellTypeTableWithFileNames[2,1]), "data\n")
    a2=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[2,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    print(dim(a2))

    cat("reading", as.character(cellTypeTableWithFileNames[3,1]), "data\n")
    a3=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[3,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    print(dim(a3))

    cat("reading", as.character(cellTypeTableWithFileNames[4,1]), "data\n")
    a4=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[4,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    print(dim(a4))

    cat("reading", as.character(cellTypeTableWithFileNames[5,1]), "data\n")
    a5=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[5,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    print(dim(a5))

    cat("reading", as.character(cellTypeTableWithFileNames[6,1]), "data\n")
    a6=tryCatch({res=read.table(as.character(cellTypeTableWithFileNames[6,2]))}, warning = function(w) {cat("a warning was raised\n"); return(emptyDF); }, error = function(e) {return(emptyDF); }, finally = {cat("try-catch done\n")})
    print(dim(a6))

    cat("reading order matrix\n")
    orderMatrix=read.table(orderMatrixFile);
    orderMatrixSub=orderMatrix
    m=as.data.frame(orderMatrixSub[,2])
    rownames(m)=orderMatrixSub[,1]
    alignedReads=read.table(alignedReadsFile, sep = "\t");
    altExons=read.table(altExonsFile)
    projections=read.table(projFile)

    projected=projections[which(projections[,4]=="exonic"),]
    #print(projected)
    rightPlotEnd = max(projected[,6])

    order1 <- getClustering(a1, m, cluster, 1)
    order2 <- getClustering(a2, m, cluster, 2)
    order3 <- getClustering(a3, m, cluster, 3)
    order4 <- getClustering(a4, m, cluster, 4)
    order5 <- getClustering(a5, m, cluster, 5)
    order6 <- getClustering(a6, m, cluster, 6)

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
              alignedReads6=order6,56889953,56890775,drawAxis=FALSE,altExons=altExons,projections=projections,
              CI=CI,specialColor="darkorange");
    warnings()
}


dev.off();
