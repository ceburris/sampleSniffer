#' Finds the Closest Genetic Match for Each Sample
#'
#' Reads in all VCF files in input directory, calculates genetic distance between each
#' pair of samples, and finds the closest genetic match for each sample. Either snps
#' file or chromosome number must be specified, but only one can be used at a time.
#' Valid values for chromosome option are numbers 1-23 and X. If output filename is
#' given, value is defaulted to sampleSniffer.pdf.
#'
#' @param input Input directory containing VCF files
#' @param snps Name of file containing the SNP IDs to be used
#' @param chromsome Chromosome number to be used
#' @param output Name of output file
#' @export
sampleSniffer <- function(input=NULL,snps=NULL,chromosome=NULL,output="sampleSniffer.pdf"){
    if (is.null(snps) && is.null(chromosome)){
        stop("Either snps or chromosome option must be supplied")
    } else if (!is.null(snps) && !is.null(chromosome)){
        stop("Only one of either snps or chromosome options can be used")
    } else if (is.null(input)){
        print("Input directory is defaulted to current directory")
        input = "."
    }

    if (!is.null(snps)){
        snps.table = read.table(snps)
        snps = snps.table$V1
        calls = getSnps(snps=snps,input=input)
    } else {
        calls = getSnps(subsetByChromosome=paste("chr",chromosome,sep=""),input=input)
    }

    calls = calls[,-1]
    distance = calculateDistance(calls)

    m <- apply(distance, 1, getClosest)
    print(cbind(names(distance), names(distance)[m]))

    pheatmap::pheatmap(distance,display_numbers=TRUE)
}
