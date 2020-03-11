#' Taxonomic and supplementary information for \emph{nifH} reference sequences.
#'
#' The dataset of \emph{nifH} sequences used for: 1) constructing reference alignment and tree, 2) evaluating \code{ppit} accuracy, and 3) taxonomic inferencing. Sequences curated from GenBank using ARBitrator (Heller \emph{et al}., 2014).
#'
#' @format \code{Data frame} containing 8876 rows and 26 columns.
#' \describe{
#'   \item{Domain}{Domain of source organism}
#'   \item{Phylum}{Phylum of source organism}
#'   \item{Class}{Class of source organism}
#'   \item{Order}{Order of source organism}
#'   \item{Family}{Family of source organism}
#'   \item{Genus}{Genus of source organism}
#'   \item{Species}{Species of source organism}
#'   \item{Strain}{Strain of source organism}
#'   \item{Reference_strain}{Type/reference strains marked with "X"}
#'   \item{pid_ref}{Sequences used in percent identity calculation marked with "X"}
#'   \item{Used_in_reference_checking}{Sequences used during taxonomic inferencing marked with "X"}
#'   \item{Version_added}{Database version in which sequence was added}
#'   \item{Source_organism}{Source organism of sequence}
#'   \item{Tip_label}{Tip label on \code{nifH_reference_tree_v2}}
#'   \item{Nucleotide_accession}{Nucleotide accession of source scaffold, genome, etc.}
#'   \item{Nucleotide_accession_sequence_length..bp.}{Length of source scaffold, genome, etc. (bp)}
#'   \item{Creation_date}{Date when sequence was deposited into GenBank}
#'   \item{Nucleotide_sequence}{Nucleotide sequence of reference nifH}
#'   \item{CDS_start}{Coding start position in nucleotide accession}
#'   \item{CDS_stop}{Coding stop position in nucleotide accession}
#'   \item{Nucleotide_sequence_length}{Length of nifH reference sequence (bp)}
#'   \item{Gene_location}{Location of \emph{nifH} (i.e., chromosome, plasmid, undetermined)}
#'   \item{Protein_accession}{Protein accession number}
#'   \item{ARBitrator_search_set}{Sequences used for initial ARBitrator search marked with "X"}
#'   \item{Alignment_seed_set}{Sequences used MAFFT-DASH seed alignment}
#'   \item{Suspected_NifH_homolog}{Suspected NifH homologs marked with "X"}
#'   ...
#'   }
#'
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In. prep.
"nifH_reference_taxonomy_v2"

#' Reference \emph{nifH} gene tree (v.1.2.0)
#'
#' Reference \emph{nifH} gene tree to use for inserting \emph{nifH} query sequences. For details about tree construction,
#' see Kapili & Dekas, in prep.
#'
#' @format \code{phylo} object with 8876 tips.
#' @source PPIT: an R package for inferring microbial taxonomy from nifH sequences
"nifH_reference_tree_v2"

#' Reference \emph{nifH} nucleotide alignment (v.1.2.0)
#'
#' Reference \emph{nifH} nucleotide alignment to use for inserting \emph{nifH} query sequences. For details about alignment creation,
#' see Kapili & Dekas, in prep.
#'
#' @format \code{DNAMultipleAlignment} object with 8876 sequences
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In prep.
"nifH_reference_alignment_v2"

#' RAxML info file for SEPP placement (v.1.2.0)
#'
#' RAxML info file generated during construction of \code{nifH_reference_tree_v2}. To be used with SEPP (Mirarab et al., 2012) for insertion of \emph{nifH} query sequences
#' into reference \emph{nifH} alignment and tree
#'
#' @format \code{Data frame} object with 11477 rows and 1 column.
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In prep.
"nifH_reference_RAxML_info_v2"

#' Patristic distance and pairwise percent identity rank cutoffs for \emph{nifH}.
#'
#' Patristic distance and pairwise percent identity rank cutoffs for \emph{nifH}. For details
#' on calculations, see Kapili & Dekas in prep.
#'
#' @format \code{Data frame} object with 2 rows and 5 columns.
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In prep.
"nifH_cutoffs_v2"

#' Phyloseq object containing \emph{nifH} sequences analyzed in Kapili & Dekas, in prep.
#'
#' Phyloseq object containing 1245 \emph{nifH} sequences and 13 samples. For details on sequence prep, see Kapili & Dekas, in prep.
#'
#' @format \code{phyloseq-class} object with 1245 sequences and 13 samples.
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In prep.
"nifH_example_ps"

#' SEPP output alignment for PPIT tutorial.
#'
#' SEPP output alignment to be used in PPIT tutorial. Produced by inserting the \emph{nifH} ASVs from \code{nifH_example_ps} into \code{nifH_reference_alignment_v2} using SEPP (v.4.3.5).
#'
#' @format \code{DNAMultipleAlignment} containing 10121 total sequences (8876 reference \emph{nifH} sequences; 1245 \emph{nifH} ASVs).
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In prep.
"nifH_example_SEPP_alignment"

#' SEPP output tree for PPIT tutorial.
#'
#' SEPP output tree to be used in PPIT tutorial. Produced by inserting the \emph{nifH} ASVs from \code{nifH_example_ps} into \code{nifH_reference_alignment_v2} and \code{nifH_reference_tree_v2} using SEPP (v.4.3.5).
#'
#' @format \code{DNAMultipleAlignment} containing 10121 total sequences (8876 reference \emph{nifH} sequences; 1245 \emph{nifH} ASVs).
#' @source BJ Kapili and AE Dekas. PPIT: an R package for inferring microbial taxonomy from nifH sequences. In prep.
"nifH_example_SEPP_tree"
