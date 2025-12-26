library(stringi)
#run genes_that_are_cancer_markers_stability_test before!

stability_threshold_share <- 4
# if a gene appears for 1/stability_results_share
# runs of the stability test, it is stable
stability_results <- readRDS(
    "results/GSE112026/stabiblity_test.rds"
)
#creates a list of stable cancer markers
stability_results_table <- table(unlist(readRDS(
    "results/GSE112026/stabiblity_test.rds"
)))

threshold_to_be_stable <- length(stability_results) %/%
    stability_threshold_share #250 if 4

stable_markers <- names(stability_results_table)[
    which(stability_results_table >= threshold_to_be_stable)
]

sink("results/GSE112026/stable_markers.txt")
cat(stable_markers, sep = "\n")
sink()

#Now, let's check wheter it exists in usual hg19
if (!require(differential.coverage)) {
  devtools::install_github("favorov/differential.coverage")
  library(differential.coverage)
}

genes_19 <- differential.coverage::gencode_hs19_genes

stable_markers_gene_ids <- stri_sub(
  stable_markers, 1,
  stri_locate_first_fixed(stable_markers, "|")[, 1] - 1
)

#let's fix it to match gencode hs 19 names
#(as provided by differential.coverage::gencode_hs19_genes)
#C16orf73|254528 is MEIOB (NCBI 254528, ENSG00000162039)
#CXorf59|286464 is CXorf22 (CFAP47, NCBI 286464, ENSG00000165164)
#LST-3TM12|338821 is SLCO1B7 (NCBI 338821, ENSG00000205754)

gene_id_fix <- stable_markers_gene_ids
names(gene_id_fix) <- gene_id_fix
gene_id_fix["C16orf73"] <- "MEIOB"
gene_id_fix["CXorf59"] <- "CXorf22"
gene_id_fix["LST-3TM12"] <- "SLCO1B7"

#now, we have gencode hs19 gene names
stable_markers_gene_ids <- gene_id_fix[stable_markers_gene_ids] |> unname()

#which list to write to?
if_cancer <- stable_markers |> stri_detect_fixed("cancer")

sink("results/GSE112026/stable_markers_cancer_hs19.txt")
cat(stable_markers_gene_ids[if_cancer], sep = "\n")
sink()

sink("results/GSE112026/stable_markers_normal_hs19.txt")
cat(stable_markers_gene_ids[!if_cancer], sep = "\n")
sink()
