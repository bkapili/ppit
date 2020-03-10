#' Calculate pairwise percent identity.
#'
#' Wrapper function that uses the \code{dist.gene} function in \code{ape} to calculate
#' pairwise percent identities between query and reference sequences. Loci in input alignment that only contain
#' gaps are removed prior to calculation.
#'
#' @param type Type of query sequence (i.e., "partial" if amplicon or "full" if full-length). For partial sequences,
#' pairwise percent identity calculated based on loci between first and last non-gap amplicon positions in the alignment.
#' @param query Name of the query sequence to compare to references.
#' @param group Vector of reference sequence names to compare to query sequence.
#' @param alignment Nucleotide alignment used for constructing phylogenetic tree of query and reference sequences (e.g., SEPP output alignment).
#' @return An object of class \code{data frame}
#' @importMethodsFrom Biostrings as.matrix
#' @export

percent.identity <- function(type, query, group, alignment_as_matrix) {
  # Keep only the query and the closest reference
  alignment_as_matrix <- alignment_as_matrix[match(c(query,group),rownames(alignment_as_matrix)),]

  # If the query is a partial sequence, include a step where
  # we only keep the positions in the alignment that are between the first and last base of the query
  if(type == 'partial') {
    # Only keep the positions in the alignment that are between the first and last base of the query
    subset <- alignment_as_matrix[,which(alignment_as_matrix[1,] != "-")[2]:sort(which(alignment_as_matrix[1,] != "-"),decreasing = TRUE)[2]]

    # For each position in the alignment:
    delete_column <- numeric()
    for(n in 1:dim(subset)[2]) {
      # Check if both sequences have a gap at that position, and if it does, record the index of the column (the empty position)
      if(length(unique(subset[,n])) == 1 & unique(subset[,n])[1] == '-') {
        delete_column <- append(n, delete_column)
      }
    }

    # Delete the columns that are only gaps
    if(length(delete_column) > 0) {
      subset <- subset[,-delete_column]
    }

    # Calculate percent identity between these two sequences, and record the result in the corresponding cell in the dataframe
    pid <- 1-ape::dist.gene(x = subset, method = 'percentage', pairwise.deletion = FALSE)
  } else {

    # For each position in the alignment:
    delete_column <- numeric()
    for(n in 1:dim(alignment_as_matrix)[2]) {
      # Check if both sequences have a gap at that position, and if it does, record the index of the column (the empty position)
      if(length(unique(alignment_as_matrix[,n])) == 1 & unique(alignment_as_matrix[,n])[1] == '-') {
        delete_column <- append(n, delete_column)
      }
    }

    # Delete the columns that are only gaps
    if(length(delete_column) > 0) {
      subset <- alignment_as_matrix[,-delete_column]
    }

    # Calculate percent identity between these two sequences, and record the result in the corresponding cell in the dataframe
    pid <- 1-ape::dist.gene(x = subset, method = 'percentage', pairwise.deletion = FALSE)
  }
  return(pid)
}
