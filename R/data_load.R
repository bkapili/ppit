#nifH_example_SEPP_tree <- read.tree(file = '~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/tax_inference_pipeline/field_data/Analysis/output_placement_mehta_20200128.txt')
#use_data(nifH_example_SEPP_tree)

#nifH_example_SEPP_alignment <- readDNAMultipleAlignment(file = '~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/tax_inference_pipeline/field_data/Analysis/output_alignment.fasta')
#use_data(nifH_example_SEPP_alignment)

#nifH_example_ps <- readRDS(file = '~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/tax_inference_pipeline/field_data/Analysis/psMehta_filtered_20200128.rds')
#use_data(nifH_example_ps)

#nifH_reference_taxonomy_v2 <- read.csv(file ="~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/tax_inference_pipeline/pid/ref_thresholds_v101/nifH_ppit_references_v2.csv", fill = TRUE, header = TRUE)
#use_data(nifH_reference_taxonomy_v2)

#nifH_reference_alignment_v2 <- readDNAMultipleAlignment(file = '~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/nifH_database/trees/v2/PAL2NAL/nifH_full_PAL2NAL_20191104')
#use_data(nifH_reference_alignment_v2)

#nifH_reference_tree_v2 <- read.tree(file = "~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/nifH_database/trees/v2/reference_tree/RAxML_ref_tree_v102_bipartitionsBranchLabels.result")
#use_data(nifH_reference_tree_v2)

#nifH_reference_RAxML_info_v2 <- read.delim(file = "~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/nifH_database/trees/v2/reference_tree/RAxML_ref_tree_v102_info.result", sep = "\n", fill = TRUE, blank.lines.skip = FALSE)
#use_data(nifH_reference_RAxML_info_v2)

#write.table(nifH_reference_RAxML_info_v1.0.0, file = "~/Desktop/test_info_file",
#            row.names = FALSE, col.names = FALSE, quote = FALSE)

#nifH_cutoffs_v2 <- readRDS(file = '~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/tax_inference_pipeline/pid/ref_thresholds_v101/nifH_cutoffs_v2.rds')
#use_data(nifH_cutoffs_v2)

#desulfuro_asvs <- read.csv(file = "~/Documents/STANFORD/DEKAS/nifH_Methods_Dev/tax_inference_pipeline/field_data/Analysis/Desulfuromonadales_ASVs_20200205.csv", header = TRUE, row.names = 1)
#colnames(desulfuro_asvs) <- gsub("\\.", "-", colnames(desulfuro_asvs))
#use_data(desulfuro_asvs)


