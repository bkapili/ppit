#' Calculate patristic distances between query and reference sequences.
#'
#' Wrapper function that reformats the patristic distance output from the \code{cophenetic}
#' function in \code{ape}. Distance matrix output
#' reformatted to have query sequences as rows and reference sequences as columns.
#'
#' @param query Vector of query sequence names.
#' @param tree Phylogenetic tree in Newick format containing query sequences placed onto a reference tree.
#' @return An object of class \code{"dist"}

patristic.distance <- function(query, tree) {
  # Calculate patristic distances (i.e., sum of branch lengths) for each pair of tips on the tree
  df <- as.data.frame(ape::cophenetic.phylo(tree))

  # Create dataframe that contains patristic distances between just queries (as rows) and references (as columns)
  df <- df[match(query, rownames(df)),]
  df <- df[,-which(colnames(df) %in% query)]

  return(df)
}
