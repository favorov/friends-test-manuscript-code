library(dplyr)
library(friends.test)
library(tictoc)
library(purrr)
library(parallel)

# Get number of available processors
num_cores <- (detectCores() * 2) %/% 3
#we are unmpolite here,
#using 2/3 of cores

cat("ncores: ", num_cores, "\n")


#actally, we need it only for sample names and for sanity check
expression_filename <- "data/GSE112026/GSE112026_HPVOP_RSEMNorm.txt"
if (file.exists(expression_filename)) {
    tic("reading GSE112026...")
    GSE112026_expression_data <-
        as.matrix(
            read.table(
                expression_filename,
                header = TRUE,
                row.names = 1,
                check.names = FALSE
            )
        )
    # remove the last 4 columns using %>% filter
    nc <- ncol(GSE112026_expression_data)
    GSE112026_expression_data <-
        GSE112026_expression_data[, -((nc - 3):nc)]
    toc()
} else {
    stop("File not found: ", expression_filename)
}


samples_are_genes_friends_file <-
    "results/GSE112026/GSE112026_samples_are_genes_friends.rds"

# Check if the RDS friends file exists, if not,
# calculate and save it as RDS

if (file.exists(samples_are_genes_friends_file)) {
    # Load the friends test from the RDS file
    tic("loading friends...")
    samples_are_genes_friends <-
        readRDS(samples_are_genes_friends_file)
    toc()
} else {
    stop(paste0("Cannot load ", samples_are_genes_friends_file,
                "run genes_that_are_markers.r first"))
}

#convert the samples_are_genes_friends data frame
#to list of lists grouped by the marker
samples_are_genes_friends_list <-
    split(
        samples_are_genes_friends,
        samples_are_genes_friends$marker
    )

samples <- colnames(GSE112026_expression_data)

cancer_samples_no <- length(grep("^HNSCC", samples))

friends_share_threshold <- 4

min_cancer_friends <- cancer_samples_no %/% friends_share_threshold
min_normal_friends <- (length(samples) - cancer_samples_no) %/%
    friends_share_threshold

#permutations: we permute the samples and we trace whether
#the cancer genes population decreses in permuted

null_distribution_file <-
    "results/GSE112026/gene_thet_are_markers_no_null_distribution.rds"

if (!file.exists(null_distribution_file)) {
    permutation_number <- 1000000
    perm_reports <- 100


    perm_report_interval <- permutation_number / perm_reports

    set.seed(42) # Set seed for reproducibility

    #make it perrallel
    null_distribution <-
        seq_len(permutation_number) |> mclapply(function(nper) {
            # Perform permutation test
            permuted_samples <-
                sample(samples, length(samples), replace = FALSE)
            #the list is like:
            #$`GPR77|27202`
            #          marker  friend marker.index friend.index friend.rank
            #6936 GPR77|27202 HNSCC36         7364           34           1
            #6937 GPR77|27202 UPPP108         7364           55           2
            #6938 GPR77|27202 HNSCC28         7364           27           3

            marker.genes.in.perm <-
                purrr::imap(samples_are_genes_friends_list, ~{
                    # .x is the element value
                    # .y is the element name
                    permuted.friends <- permuted_samples[.x$friend.index]
                    if (length(permuted.friends) >= min_cancer_friends &&
                            length(grep("^UPPP", permuted.friends)) == 0
                    ) {
                        paste0(.y, ": cancer")
                    } else if (length(permuted.friends) >= min_normal_friends &&
                            length(grep("^HNSCC", permuted.friends)) == 0
                    ) {
                        paste0(.y, ": normal")
                    } else {
                        NULL
                    }
                })

            if (nper %% perm_report_interval == 0) {
                cat(
                    paste0(
                        "Permutation ",
                        nper, "\n"
                    )
                )
            }
            length(unlist(marker.genes.in.perm))
        }, mc.cores = num_cores)

    null_distribution <- unlist(null_distribution)

    saveRDS(
        table(null_distribution),
        null_distribution_file
    )
} else {
    cat(paste0("The file ", null_distribution_file, " already exists"))
}
