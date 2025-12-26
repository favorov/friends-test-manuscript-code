library(dplyr)
library(friends.test)
library(tictoc)
library(purrr)


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

if (!file.exists(samples_are_genes_friends_file)) {
    # Calculate the friends test
    tic("friending..")
    samples_are_genes_friends <-
        friends.test(
            GSE112026_expression_data, max.friends = "all"
        )
    toc()
    tic("saving...")
    saveRDS(
        samples_are_genes_friends,
        samples_are_genes_friends_file
    )
    toc()
} else {
    # Load the friends test from the RDS file
    tic("loading friends...")
    samples_are_genes_friends <-
        readRDS(samples_are_genes_friends_file)
    toc()
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


#the list is like:
#$`GPR77|27202`
#          marker  friend marker.index friend.index friend.rank
#6936 GPR77|27202 HNSCC36         7364           34           1
#6937 GPR77|27202 UPPP108         7364           55           2
#6938 GPR77|27202 HNSCC28         7364           27           3

marker.genes <- purrr::imap(samples_are_genes_friends_list, ~{
    # .x is the element value
    # .y is the element name
    friends <- .x$friend
    if (length(friends) >= min_cancer_friends &&
            length(grep("^UPPP", friends)) == 0
    ) {
        paste0(.y, ": cancer")
    } else if (length(friends) >= min_normal_friends &&
            length(grep("^HNSCC", friends)) == 0
    ) {
        paste0(.y, ": normal")
    } else {
        NULL
    }
})

marker.genes <- unlist(marker.genes)

saveRDS(
    marker.genes,
    "results/GSE112026/marker_genes_one_run.rds"
)