import gzip, argparse

def findInfo(annoFile, gene):
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
                if "gene_name" in columns[i] and gene in columns[i+1]:#len(re.findall(geneMatch, columns[i+1])) > 0:
                    annoContents.append(columns)

    minimum = 1000000000
    maximum = -1
  
    # Find start (minimum) and end (maximum) of all reads for the gene of interest
    for i in annoContents:
        if int(i[3]) < minimum:
            minimum = int(i[3])
        if int(i[4]) > maximum:
            maximum = int(i[4])
  
    chromosome = annoContents[0][0]
  
    return [chromosome, minimum, maximum]

def main(args):
    annotationFile = args['gencodeAnno'][0]
    
    gene=args['gene'][0]
    
    outputDir=args['outputDir'][0]
      
    annoFile = gzip.open(annotationFile, "rb")

    info = list(findInfo(annoFile, gene))

    outputInfo = gene + ".info.tab"
    geneOutputDir = outputDir + gene
    outputPath1 = geneOutputDir + "/" + outputInfo
    infoOutput = open(outputPath1, "w")

    for i in info:
        infoOutput.write(str(i))
        infoOutput.write("\t")
    infoOutput.write("\n")
        
    infoOutput.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Finding gene info')
    parser.add_argument('--gencodeAnno', required=True, nargs=1)
    parser.add_argument('--gene', required=True, nargs=1)
    parser.add_argument('--outputDir', required=True, nargs=1)
    args = vars(parser.parse_args())
    main(args)
