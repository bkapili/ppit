#' Infer host identities of query sequences.
#'
#' Function that infers taxonomic identity of marker gene source organisms based on phylogenetic placement and sequence identity. For a tutorial, type \code{Vignette('ppit')}.
#'
#' @param type Either "partial" or "full".
#' @param query Vector of query sequence names.
#' @param tree Phylogenetic tree in Newick format containing query sequences placed onto a reference tree.
#' @param alignment Nucleotide alignment used for constructing phylogenetic tree of query and reference sequences (e.g., SEPP output alignment).
#' @param opt_thresh Phylogenetic neighborhood threshold.
#' @param cutoffs Data frame containing patristic distance and pairwise percent identity cutoffs for each taxonomic rank.
#' @return An object of class \code{data frame}
#' @export

ppit <- function(type, query, tree, alignment, taxonomy, opt_thresh, cutoffs) {
  # Define dataframe to hold taxonomic inferences
  df_inferences <- as.data.frame(matrix(nrow = length(query), ncol = 14)) # Define dataframe to hold taxonomic inferences
  df_inferences[,] <- '' # Make all cells empty
  rownames(df_inferences) <- query # Name the rows using the amplicon names
  colnames(df_inferences) <- c('Suspected_homolog', 'Percent_id_fail', 'Pat_dist_fail', 'Potential_HGT','Percent_id', 'Pat_dist', 'Domain', 'Phylum' ,'Class', 'Order', 'Family', 'Genus', 'Plasmid', 'Num_taxa') # Name the columns according to the information we want to obtain

  # Create vectors containing suspected homologs and sequences located on plasmids
  homologs <- taxonomy$Tip_label[which(taxonomy$Suspected_NifH_homolog == 'X')] # List of suspected homologs
  plasmids <- taxonomy$Tip_label[which(taxonomy$Gene_location == 'plasmid')] # List of sequences on plasmids

  # Calculate patristic distances
  p_dist <- patristic.distance(query, tree) # Calculate patristic distances

  # Coerce DNAMultipleAlignment object to matrix, so we can calculate % identity using ape::dist.gene
  alignment_as_matrix <- as.matrix(alignment)

  # Here, we will loop through each amplicon and first check if the closest reference
  # sequence is a suspected homolog. If it is, we will mark that amplicon as a suspected homolog
  # and proceed to the next amplicon. If not, we will gather the references that are
  # within the patristic distance cutoff for each rank. At each rank, we will
  # compare the taxonomic identity at that rank of all the gathered references
  # and record the identity only if they're all the same. Once there are no references
  # below the cutoff for a rank OR there is an inconsistency in taxonomic identity,
  # the loop is stopped and we move to the next amplicon.

  for(x in query) {
    print(paste("Starting inference", match(x, query), "of", length(query)))

    failed <- as.character() # Reset the failed ranks tracker
    row_ind <- which(rownames(p_dist) == x) # Record the row index in p_dist that corresponds to the amplicon

    min_pdist <- min(p_dist[row_ind,]) # Record the minimum patristic distance
    df_inferences$Pat_dist[row_ind] <- min_pdist

    min_ind <- which(p_dist[row_ind,] == min_pdist) # Record the column index that corresponds to the closest reference
    closest_ref <- colnames(p_dist)[min_ind][1] # Record the name of the closest reference. If there's two, take the first one

    est_max_id <- percent.identity(type, x, closest_ref, alignment_as_matrix) # Record pairwise percent identity with closest neighbor

    df_inferences$Percent_id[row_ind] <- est_max_id

    diverged_nifH_ref <- "Gaiellales_bacterium--42930-43757-QZKG01000006_1-RJQ47724_1" # This is the name of the most diverged nifH sequence, as determined by position on tree

    if(closest_ref %in% homologs | (closest_ref == diverged_nifH_ref & est_max_id < cutoffs[2,1])) { # If the closest reference is in the list of suspected homologs OR the closest ref is the most diverged nifH reference but is below phylum %id threshold,
      df_inferences$Suspected_homolog[row_ind] <- 'X' # Mark this amplicon as a suspected homolog and move to next amplicon
      #print('Done.')
    } else if(est_max_id < cutoffs[2,1]){ # If the estimated maximum percent id is lower than the phlyum threshold
      df_inferences$Percent_id_fail[row_ind] <- 'X' # Mark that this amplicon was below minimum percent identity cutoff
      #print('Done.')
    } else if(min_pdist > opt_thresh) { # If the minimum patristic distance is higher than the optimum threshold
      df_inferences$Pat_dist_fail[row_ind] <- 'X' # Mark this this amplicon was above the minimum patristic distance cutoff
      #print('Done.')
    } else { # Otherwise, proceed with the source organism inferencing
      # Record the taxonomic ranks at which we can infer identity based on the lowest patristic distance and estimated max id
      possible_ranks <- colnames(cutoffs[which(min_pdist < cutoffs[1,] & est_max_id > cutoffs[2,])])

      opt_set <- p_dist[row_ind, which(p_dist[row_ind,] < opt_thresh)] # Subset patristic distances to those within the optimal patristic distance threshold

      # At the lowest taxonomic rank possible, check if all the rank names are consistent
      lowest_pdist_rank <- possible_ranks[length(possible_ranks)] # Record the lowest possible rank we can infer at based on the smallest patristic distance
      lowest_ids <- colnames(p_dist)[which(p_dist[row_ind,] < cutoffs[1,match(lowest_pdist_rank, colnames(cutoffs))] & p_dist[row_ind,] < opt_thresh)] # Record the sequence names within the patristic distance cutoff for that rank
      ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences

      # If any are Betaproteobacteria, replace class with "Gammaproteobacteria" since they now appear to be an order
      # within the Gammaproteobacteria (Parks et al., 2018)
      if("Betaproteobacteria" %in% ref_tax$Class) {
        ref_tax$Class[which(ref_tax$Class == "Betaproteobacteria")] <- "Gammaproteobacteria"
      }

      # Determine all the ranks at which those taxon names are consistent, neglecting blank records
      df_consistent <- as.data.frame(matrix(nrow = 1, ncol = length(possible_ranks)))
      colnames(df_consistent) <- possible_ranks

      consistency.check <- function(rank) {
        a <- length(which(unique(ref_tax[,match(rank, colnames(ref_tax))]) != ''))

        return(a)
      }
      df_consistent[1,] <- lapply(possible_ranks,consistency.check)

      # From the set of sequences at the lowest rank at which we can hope to infer, determine the lowest rank at which all the names are consistent
      consistent_ranks <- colnames(df_consistent[which(df_consistent == 1)])
      lowest_consistent_rank <- consistent_ranks[length(consistent_ranks)]


      if("Phylum" %in% consistent_ranks == FALSE) { # If there are no ranks at phylum or lower that are consistent
        df_inferences$Potential_HGT[row_ind] <- 'X' # Mark that this sequence is in a region of the tree with potential HGT and move on to the next sequence
        #print('Done.')

      } else if(length(lowest_consistent_rank) > 0) { # But if there is at least one rank that is taxonomically consistent below phylum
        if(lowest_pdist_rank == lowest_consistent_rank) { # And if the consistent rank is the same threshold rank we used for creating the sequence subset
          ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Reset reference taxonomies in case any were Betaproteobacteria
          df_inferences[row_ind,match('Domain', colnames(df_inferences)):match(lowest_consistent_rank, colnames(df_inferences))] <- as.matrix(ref_tax[1,1:match(lowest_consistent_rank, colnames(ref_tax))]) # Record the taxonomic profile of those in-group as the taxonomic inference

          ref_species <- taxonomy$Species[which(taxonomy$Tip_label %in% lowest_ids)] # Record references' species names
          df_inferences$Num_taxa[row_ind] <- length(unique(ref_species)) # Record how many unique species were used to infer the identity

          location <- taxonomy$Gene_location[match(lowest_ids, taxonomy$Tip_label)] # Record where these sequences were located (i.e., chromosome, plasmid, undetermined)

          df_inferences$Percent_id[row_ind] <- est_max_id
          df_inferences$Pat_dist[row_ind] <- min_pdist

          # If one of the reference sequences is located on a plasmid,
          if("plasmid" %in% location) {
            df_inferences$Plasmid[row_ind] <- 'X' # Mark that taxonomic inference in the 'Plasmid' column with an 'X'
          }
          #print('Done.')
        } else { # But if the consistent rank is not the same threshold rank we used for creating the sequence subset OR if the closest sequence is not above the required percent identity threshold
          # Create a new subset using the threshold at the consistent rank
          # At the lowest taxonomic rank possible, check if all the rank names are consistent
          lowest_ids <- colnames(p_dist)[which(p_dist[row_ind,] < cutoffs[1,match(lowest_consistent_rank, colnames(cutoffs))] & p_dist[row_ind,] < opt_thresh)] # Record the sequence names within the patristic distance cutoff for that rank
          ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences

          # If any are Betaproteobacteria, replace class with "Gammaproteobacteria" since they now appear to be an order
          # within the Gammaproteobacteria (Parks et al., 2018)
          if("Betaproteobacteria" %in% ref_tax$Class) {
            ref_tax$Class[which(ref_tax$Class == "Betaproteobacteria")] <- "Gammaproteobacteria"
          }

          # Determine all the ranks at which those taxon names are consistent, neglecting blank records
          df_consistent <- as.data.frame(matrix(nrow = 1, ncol = length(possible_ranks)))
          colnames(df_consistent) <- possible_ranks

          df_consistent[1,] <- lapply(possible_ranks,consistency.check)

          # From the set of sequences at the lowest rank at which we can hope to infer, determine the lowest rank at which all the names are consistent
          consistent_ranks <- colnames(df_consistent[which(df_consistent == 1)])

          if("Phylum" %in% consistent_ranks == FALSE) { # If there are no ranks at phylum or lower that are consistent, do not infer taxonomy and move to the next amplicon.
            df_inferences$Potential_HGT[row_ind] <- 'X' # Mark that this sequence is in a region of the tree with potential HGT and move on to the next sequence
            #print('Done.')

          } else if(length(consistent_ranks[length(consistent_ranks)]) > 0) { # But if there is at least one rank that is taxonomically consistent below phylum
            if(lowest_consistent_rank == consistent_ranks[length(consistent_ranks)]) { # And if the consistent rank is the same threshold rank we used for creating the sequence subset AND the closest sequence is above the required percent identity threshold
              ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Reset reference taxonomies in case any were Betaproteobacteria
              df_inferences[row_ind,match('Domain', colnames(df_inferences)):match(lowest_consistent_rank, colnames(df_inferences))] <- as.matrix(ref_tax[1,1:match(lowest_consistent_rank, colnames(ref_tax))]) # Record the taxonomic profile of those in-group as the taxonomic inference

              ref_species <- taxonomy$Species[which(taxonomy$Tip_label %in% lowest_ids)] # Record references' species names
              df_inferences$Num_taxa[row_ind] <- length(unique(ref_species)) # Record how many unique species were used to infer the identity

              location <- taxonomy$Gene_location[match(lowest_ids, taxonomy$Tip_label)] # Record where these sequences were located (i.e., chromosome, plasmid, undetermined)

              df_inferences$Percent_id[row_ind] <- est_max_id
              df_inferences$Pat_dist[row_ind] <- min_pdist
              # If one of the reference sequences is located on a plasmid,
              if("plasmid" %in% location) {
                df_inferences$Plasmid[row_ind] <- 'X' # Mark that taxonomic inference in the 'Plasmid' column with an 'X'
              }
              #print('Done.')

            } else { # But if the consistent rank is not the same threshold rank we used for creating the sequence subset OR if the closest sequence is not above the required percent identity threshold
              # Create a new subset using the threshold at the consistent rank
              # At the lowest taxonomic rank possible, check if all the rank names are consistent
              lowest_consistent_rank <- consistent_ranks[length(consistent_ranks)]
              lowest_ids <- colnames(p_dist)[which(p_dist[row_ind,] < cutoffs[1,match(lowest_consistent_rank, colnames(cutoffs))] & p_dist[row_ind,] < opt_thresh)] # Record the sequence names within the patristic distance cutoff for that rank
              ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences

              # If any are Betaproteobacteria, replace class with "Gammaproteobacteria" since they now appear to be an order
              # within the Gammaproteobacteria (Parks et al., 2018)
              if("Betaproteobacteria" %in% ref_tax$Class) {
                ref_tax$Class[which(ref_tax$Class == "Betaproteobacteria")] <- "Gammaproteobacteria"
              }

              # Determine all the ranks at which those taxon names are consistent, neglecting blank records
              df_consistent <- as.data.frame(matrix(nrow = 1, ncol = length(possible_ranks)))
              colnames(df_consistent) <- possible_ranks

              df_consistent[1,] <- lapply(possible_ranks,consistency.check)

              # From the set of sequences at the lowest rank at which we can hope to infer, determine the lowest rank at which all the names are consistent
              consistent_ranks <- colnames(df_consistent[which(df_consistent == 1)])

              if("Phylum" %in% consistent_ranks == FALSE) { # If there are no ranks at phylum or lower that are consistent, do not infer taxonomy and move to the next amplicon.
                df_inferences$Potential_HGT[row_ind] <- 'X' # Mark that this sequence is in a region of the tree with potential HGT and move on to the next sequence
                #print('Done.')

              } else if(length(lowest_consistent_rank) > 0) { # But if there is at least one rank that is taxonomically consistent below phylum
                if(lowest_consistent_rank == consistent_ranks[length(consistent_ranks)]) { # And if the consistent rank is the same threshold rank we used for creating the sequence subset AND the closest sequence is above the required percent identity threshold
                  ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences
                  df_inferences[row_ind,match('Domain', colnames(df_inferences)):match(lowest_consistent_rank, colnames(df_inferences))] <- as.matrix(ref_tax[1,1:match(lowest_consistent_rank, colnames(ref_tax))]) # Record the taxonomic profile of those in-group as the taxonomic inference

                  ref_species <- taxonomy$Species[which(taxonomy$Tip_label %in% lowest_ids)] # Record references' species names
                  df_inferences$Num_taxa[row_ind] <- length(unique(ref_species)) # Record how many unique species were used to infer the identity

                  location <- taxonomy$Gene_location[match(lowest_ids, taxonomy$Tip_label)] # Record where these sequences were located (i.e., chromosome, plasmid, undetermined)

                  df_inferences$Percent_id[row_ind] <- est_max_id
                  df_inferences$Pat_dist[row_ind] <- min_pdist
                  # If one of the reference sequences is located on a plasmid,
                  if("plasmid" %in% location) {
                    df_inferences$Plasmid[row_ind] <- 'X' # Mark that taxonomic inference in the 'Plasmid' column with an 'X'
                  }
                  #print('Done.')
                }
              } else { # But if the consistent rank is not the same threshold rank we used for creating the sequence subset OR if the closest sequence is not above the required percent identity threshold
                # Create a new subset using the threshold at the consistent rank
                # At the lowest taxonomic rank possible, check if all the rank names are consistent
                lowest_consistent_rank <- consistent_ranks[length(consistent_ranks)]
                lowest_ids <- colnames(p_dist)[which(p_dist[row_ind,] < cutoffs[1,match(lowest_consistent_rank, colnames(cutoffs))] & p_dist[row_ind,] < opt_thresh)] # Record the sequence names within the patristic distance cutoff for that rank
                ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences

                # If any are Betaproteobacteria, replace class with "Gammaproteobacteria" since they now appear to be an order
                # within the Gammaproteobacteria (Parks et al., 2018)
                if("Betaproteobacteria" %in% ref_tax$Class) {
                  ref_tax$Class[which(ref_tax$Class == "Betaproteobacteria")] <- "Gammaproteobacteria"
                }

                # Determine all the ranks at which those taxon names are consistent, neglecting blank records
                df_consistent <- as.data.frame(matrix(nrow = 1, ncol = length(possible_ranks)))
                colnames(df_consistent) <- possible_ranks

                df_consistent[1,] <- lapply(possible_ranks,consistency.check)

                # From the set of sequences at the lowest rank at which we can hope to infer, determine the lowest rank at which all the names are consistent
                consistent_ranks <- colnames(df_consistent[which(df_consistent == 1)])

                if("Phylum" %in% consistent_ranks == FALSE) { # If there are no ranks at phylum or lower that are consistent, do not infer taxonomy and move to the next amplicon.
                  df_inferences$Potential_HGT[row_ind] <- 'X' # Mark that this sequence is in a region of the tree with potential HGT and move on to the next sequence
                  #print('Done.')

                } else if(length(lowest_consistent_rank) > 0) { # But if there is at least one rank that is taxonomically consistent below phylum
                  if(lowest_consistent_rank == consistent_ranks[length(consistent_ranks)]) { # And if the consistent rank is the same threshold rank we used for creating the sequence subset AND the closest sequence is above the required percent identity threshold
                    ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences
                    df_inferences[row_ind,match('Domain', colnames(df_inferences)):match(lowest_consistent_rank, colnames(df_inferences))] <- as.matrix(ref_tax[1,1:match(lowest_consistent_rank, colnames(ref_tax))]) # Record the taxonomic profile of those in-group as the taxonomic inference

                    ref_species <- taxonomy$Species[which(taxonomy$Tip_label %in% lowest_ids)] # Record references' species names
                    df_inferences$Num_taxa[row_ind] <- length(unique(ref_species)) # Record how many unique species were used to infer the identity

                    location <- taxonomy$Gene_location[match(lowest_ids, taxonomy$Tip_label)] # Record where these sequences were located (i.e., chromosome, plasmid, undetermined)

                    df_inferences$Percent_id[row_ind] <- est_max_id
                    df_inferences$Pat_dist[row_ind] <- min_pdist
                    # If one of the reference sequences is located on a plasmid,
                    if("plasmid" %in% location) {
                      df_inferences$Plasmid[row_ind] <- 'X' # Mark that taxonomic inference in the 'Plasmid' column with an 'X'
                    }
                    #print('Done.')
                  }
                } else { # But if the consistent rank is not the same threshold rank we used for creating the sequence subset OR if the closest sequence is not above the required percent identity threshold
                  # Create a new subset using the threshold at the consistent rank
                  # At the lowest taxonomic rank possible, check if all the rank names are consistent
                  lowest_consistent_rank <- consistent_ranks[length(consistent_ranks)]
                  lowest_ids <- colnames(p_dist)[which(p_dist[row_ind,] < cutoffs[1,match(lowest_consistent_rank, colnames(cutoffs))] & p_dist[row_ind,] < opt_thresh)] # Record the sequence names within the patristic distance cutoff for that rank
                  ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences

                  # If any are Betaproteobacteria, replace class with "Gammaproteobacteria" since they now appear to be an order
                  # within the Gammaproteobacteria (Parks et al., 2018)
                  if("Betaproteobacteria" %in% ref_tax$Class) {
                    ref_tax$Class[which(ref_tax$Class == "Betaproteobacteria")] <- "Gammaproteobacteria"
                  }

                  # Determine all the ranks at which those taxon names are consistent, neglecting blank records
                  df_consistent <- as.data.frame(matrix(nrow = 1, ncol = length(possible_ranks)))
                  colnames(df_consistent) <- possible_ranks

                  df_consistent[1,] <- lapply(possible_ranks,consistency.check)

                  # From the set of sequences at the lowest rank at which we can hope to infer, determine the lowest rank at which all the names are consistent
                  consistent_ranks <- colnames(df_consistent[which(df_consistent == 1)])

                  if("Phylum" %in% consistent_ranks == FALSE) { # If there are no ranks at phylum or lower that are consistent, do not infer taxonomy and move to the next amplicon.
                    df_inferences$Potential_HGT[row_ind] <- 'X' # Mark that this sequence is in a region of the tree with potential HGT and move on to the next sequence
                    #print('Done.')

                  } else if(length(lowest_consistent_rank) > 0) { # But if there is at least one rank that is taxonomically consistent below phylum
                    if(lowest_consistent_rank == consistent_ranks[length(consistent_ranks)]) { # And if the consistent rank is the same threshold rank we used for creating the sequence subset AND the closest sequence is above the required percent identity threshold
                      ref_tax <- taxonomy[which(taxonomy$Tip_label %in% lowest_ids),1:6] # Retrieve the full taxonomies for those sequences
                      df_inferences[row_ind,match('Domain', colnames(df_inferences)):match(lowest_consistent_rank, colnames(df_inferences))] <- as.matrix(ref_tax[1,1:match(lowest_consistent_rank, colnames(ref_tax))]) # Record the taxonomic profile of those in-group as the taxonomic inference

                      ref_species <- taxonomy$Species[which(taxonomy$Tip_label %in% lowest_ids)] # Record references' species names
                      df_inferences$Num_taxa[row_ind] <- length(unique(ref_species)) # Record how many unique species were used to infer the identity

                      location <- taxonomy$Gene_location[match(lowest_ids, taxonomy$Tip_label)] # Record where these sequences were located (i.e., chromosome, plasmid, undetermined)

                      df_inferences$Percent_id[row_ind] <- est_max_id
                      df_inferences$Pat_dist[row_ind] <- min_pdist
                      # If one of the reference sequences is located on a plasmid,
                      if("plasmid" %in% location) {
                        df_inferences$Plasmid[row_ind] <- 'X' # Mark that taxonomic inference in the 'Plasmid' column with an 'X'
                      }
                      #print('Done.')
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(df_inferences)
}
