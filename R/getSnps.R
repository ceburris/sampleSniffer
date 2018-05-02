#' Get Table of SNPs by Sample
#'
#' This funtion reads a set of VCF files and creates a table of SNPs with SNP IDs as rows and
#' sample names (one per file) as columns
#'
#' @param subsetByChromosome Chromosome number to use for comparison, left NULL if using SNP list
#' @param snps List of SNP IDs to use for comparison
#' @param input Input Directory Containing annotated VCF files
#' @return Matrix containing SNP calls for all samples
getSnps <- function(subsetByChromosome=NULL,snps=snps,input=input){

options(stringsAsFactors = FALSE)

files = list.files(input,pattern="*.vcf", full.names=TRUE)
filenames = list.files(input,pattern="*.vcf")

  for (i in 1:length(files)){
    file <- files[i]
    vcf <- vcfR::read.vcfR(file)
    vcf <- vcf[which(!is.na(vcfR::getID(vcf))),]
    vcf <- vcf[which(!duplicated(vcfR::getID(vcf))),]

    if (is.null(subsetByChromosome)) {
      ids <- vcfR::getID(vcf)
      snp.sub <- vcf[ids %in% snps,]
    } else {
      chromosomes <- vcfR::getCHROM(vcf)
      snp.sub <- vcf[chromosomes == subsetByChromosome,]
    }

    snp.gt <- vcfR::extract.gt(snp.sub)
    snp.ids <- vcfR::getID(snp.sub)

    sub.final <- data.frame(ID=as.character(snp.ids),calls=as.character(snp.gt))
    sub.final <- sub.final[!is.na(sub.final$ID) & nchar(sub.final$calls == 3), ]
    name = substring(as.character(filenames[i]),1,10)
    colnames(sub.final) = c("ID",name)

    if (i == 1){
      all.calls = sub.final
    } else {
      all.calls = dplyr::full_join(all.calls,sub.final)
    }
  }
  return(all.calls)
}
