#' Determine optimal phylogenetic neighborhood threshold.
#'
#' Determines patristic distance threshold that maximizes the number of host inferences
#' returned at the phylum rank.
#'
#' @param query Vector of query sequence names.
#' @param tree Phylogenetic tree in Newick format containing query sequences placed onto a reference tree.
#' @param alignment Nucleotide alignment used for constructing phylogenetic tree of query and reference sequences (e.g., SEPP output alignment).
#' @param taxonomy Data frame containing taxonomic information and genomic location for each reference sequence.
#' @param cutoffs Data frame containing patristic distance and pairwise percent identity cutoffs for each taxonomic rank.
#' @return An object of class \code{numeric}
#' @export

tree.partition <- function(type, query, tree, alignment, taxonomy, cutoffs) {
  print('Patristic distance optimization: round 1 of 4')
  df_num <- as.data.frame(matrix(nrow = 0, ncol = 2))
  colnames(df_num) <- c('opt_thresh', 'num_phyl')

  df_temp <- as.data.frame(matrix(nrow = 1, ncol = 2))
  colnames(df_temp) <- c('opt_thresh', 'num_phyl')

  c <- 1
  k_vals <- as.numeric(c('0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'))

  eval <- ppit(type=type, query=query, tree=tree, alignment=alignment, taxonomy=taxonomy, opt_thresh = k_vals[c], cutoffs = cutoffs)
  num_phyl <- length(which(eval$Phylum != ''))
  df_num[c,] <- c(k_vals[c], num_phyl)
  c <- c+1

  eval <- ppit(type, query, tree, alignment, taxonomy, k_vals[c], cutoffs)
  num_phyl <- length(which(eval$Phylum != ''))
  df_num[c,] <- c(k_vals[c], num_phyl)

  while(num_phyl >= df_num$num_phyl[c-1]) {
    c <- c+1

    eval <- ppit(type, query, tree, alignment, taxonomy, k_vals[c], cutoffs)
    num_phyl <- length(which(eval$Phylum != ''))

    df_num[c,] <- c(k_vals[c], num_phyl)
  }

  print('Patristic distance optimization: round 2 of 4')
  max <- df_num$opt_thresh[which(df_num$num_phyl == max(df_num$num_phyl))]
  k_vals <- unique(c(max-0.05, max+0.05))

  for(x in k_vals) {
    eval <- ppit(type, query, tree, alignment, taxonomy, x, cutoffs)
    num_phyl <- length(which(eval$Phylum != ''))
    df_temp[1,] <- c(x, num_phyl)
    df_num <- rbind(df_num, df_temp)
  }

  print('Patristic distance optimization: round 3 of 4')
  max <- df_num$opt_thresh[which(df_num$num_phyl == max(df_num$num_phyl))]
  k_vals <- unique(c(max-0.03, max+0.03))


  for(x in k_vals) {
    eval <- ppit(type, query, tree, alignment, taxonomy, x, cutoffs)
    num_phyl <- length(which(eval$Phylum != ''))
    df_temp[1,] <- c(x, num_phyl)
    df_num <- rbind(df_num, df_temp)
  }

  print('Patristic distance optimization: round 4 of 4')
  max <- df_num$opt_thresh[which(df_num$num_phyl == max(df_num$num_phyl))]
  k_vals <- c(max-0.02, max-0.01, max+0.01, max+0.02)

  for(x in k_vals) {
    eval <- ppit(type, query, tree, alignment, taxonomy, x, cutoffs)
    num_phyl <- length(which(eval$Phylum != ''))
    df_temp[1,] <- c(x, num_phyl)
    df_num <- rbind(df_num, df_temp)
  }

  # Return the patristic distance that maximizes the number of phylum inferences. If there are ties,
  # choose the one with the largest value. That value uses the most data to draw inferences.
  opt <- df_num[order(-df_num[,2],-df_num[,1]),][1,1]

  # If the optimal patristic distance is smaller than the genus threshold, print that the genus threshold will be used as the cutoff.
  if(opt < cutoffs[1,5]) {
    print(paste0("Optimal threshold (", opt, ") < genus cutoff. Using genus patristic distance cutoff as phylogenetic neighborhood."))
    opt <- cutoffs[1,5]
  }
  return(opt)
}
