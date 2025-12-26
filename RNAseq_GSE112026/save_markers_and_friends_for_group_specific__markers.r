library(tictoc)
library(stringi)
library(purrr)


samples_are_genes_friends_file <-
    "results/GSE112026/GSE112026_samples_are_genes_friends.rds"

tic("loading friends git in the first run...")
samples_are_genes_friends <-
    readRDS(samples_are_genes_friends_file)
toc()

tic("loading stable markers ...")
stable_markers <-
    read.csv("results/GSE112026/stable_markers.txt", header = FALSE)

stable_markers <-
    map_chr(stri_split_fixed(stable_markers[[1]],": "),\(x) x[1])
toc()

tic("indexing")
stable_markers_with_friends <-
    samples_are_genes_friends[
        samples_are_genes_friends$marker %in% stable_markers, ]
toc()


tic("writing")
write.csv2(
    stable_markers_with_friends,
    "results/GSE112026/stable_markers_with_friends_in_The_Run.csv",
    row.names = FALSE
)
toc()
