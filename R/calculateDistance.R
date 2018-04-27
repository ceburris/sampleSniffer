#' Calculate the Genetic Distance Between Samples
#'
#' Calculates the percentage of SNP matches. The number of complete and matching
#' SNP loci is the numerator and the number of SNPs non-missing for both samples
#' is the denominator
#'
#' @param snp_matrix Matrix containing SNP calls
#' Matrix containing genetic distance between all samples
calculateDistance <- function(snp_matrix) {
  result <- matrix(data=NA, nrow=ncol(snp_matrix), ncol=ncol(snp_matrix))
  for(i in 1:ncol(snp_matrix)) {
    for(j in i:ncol(snp_matrix)) {
      if (i == j) {
        next
      }
      mat <- cbind(snp_matrix[, i], snp_matrix[, j])
      mat <- mat[complete.cases(mat), ]
      result[i,j] <- sum(mat[, 1] == mat[, 2]) / nrow(mat)
      result[j,i] <- result[i,j]
    }
  }
  result <- data.frame(result)
  rownames(result) <- colnames(snp_matrix)
  colnames(result) <- colnames(snp_matrix)
  result
}
