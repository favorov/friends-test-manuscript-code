library(ggplot2)
library(dplyr)

marker_genes_one_run <-
  readRDS("results/GSE112026/cancer_marker_genes_one_run.rds")
dist_null_lengths <- 
  readRDS("results/GSE112026/gene_thet_are_markers_no_null_distribution.rds")


stability_test <- readRDS("results/GSE112026/stabiblity_test.rds")

dist_not_null_lengths <-
  table(vapply(stability_test, length, integer(1)))



# Prepare data for three histograms:
#null
dist_null_lengths <- dist_null_lengths / sum(dist_null_lengths)
df_null <- data.frame(dist_null_lengths, colour = "white")
#not null (stability test is used)
dist_not_null_lengths <- dist_not_null_lengths / sum(dist_not_null_lengths)
df_not_null <- data.frame(dist_not_null_lengths, colour = "grey75")
#the third is even not a histogram. It is a mark or the original run
df_one_run <- data.frame(dist_not_null_lengths, colour = "grey25") |>
    filter(Var1 == length(marker_genes_one_run))
#combine them to plot
clmn <- c("Length", "Frequency", "colour")
colnames(df_null) <- colnames(df_not_null) <- colnames(df_one_run) <- clmn
df_combined <- rbind(
    df_null,
    df_not_null,
    df_one_run,df_one_run
)
# Get common y-axis limit
df_combined$Length <- df_combined$Length |> as.character() |> as.numeric()
max_len <- max(df_combined$Length)



# Create histogram
p_a <- ggplot(df_combined, aes(x = as.numeric(Length), y = Frequency)) +
    geom_col(
        fill = df_combined$colour,
        color = "black",
        alpha = 1
    ) +
    xlim(-.5, max_len+.5) +
    ylim(0, 1) +
    labs(x = "Number of marker genes", y = "Frequency") +
    theme_minimal()


# Save p_a to A6 PDF
out_dir <- "results/GSE112026"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
ggsave(
  filename = file.path(out_dir, "p_a_a6.pdf"),
  plot = p_a,
  device = "pdf",
  width = 148, height = 105, units = "mm"
)

#form a list of genes that are at least once in the stability test results
sttgenes <- stability_test |> unlist() |> sort() |> unique()
#is each in each experiment? (map)
#ok, summarise each gene over experiments (reduce)
#and, make a table

gene_occurences <- stability_test |>
    purrr::map(\(x) sttgenes %in% x) |>
    purrr::reduce(\(x,y) x+y)

#great map-reduce code.
#Actually, table(unlist(stability_test)) do the same.

df_gene_occurences <- data.frame(
    gene_occurences |> table(), 
    colour = "grey75"
)

df_gene_occurences$gene_occurences <- 
    df_gene_occurences$gene_occurences |> as.character() |> as.numeric()


# Create histogram
p_b <- ggplot(df_gene_occurences, aes(x = as.numeric(gene_occurences), y = Freq)) +
    geom_col(
        fill = df_gene_occurences$colour,
        color = "black",
        alpha = 1
    ) +
    geom_vline(xintercept = 250, linetype = "dashed") +
    labs(x = "#of tests the gene is identified in", y = "#of genes") +
    scale_x_continuous(
      limits = c(.5, 1000.5),
      breaks = seq(0, 1000, 100),
      minor_breaks = seq(50, 950, 100)
    ) +
    theme_minimal()


# Save p_b to A6 PDF
ggsave(
  filename = file.path(out_dir, "p_b_a6.pdf"),
  plot = p_b,
  device = "pdf",
  width = 148, height = 105, units = "mm"
)