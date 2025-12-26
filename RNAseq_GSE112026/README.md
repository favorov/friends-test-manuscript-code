**RNAseq_GSE112026**

## genes_that_are_markers.r
- Runs the friends test once (or reloads cached result) on the expression matrix.
- Identify group markers for cancer or normal samples, and saves that single-run marker list.

## genes_that_are_markers_stability_test.r
- Repeats the friends test 1,000 times in parallel on the same expression matrix.
- For each run, labels markers as cancer or normal based on friend composition, and saves all runs to stabiblity_test.rds for stability assessment. The runs differ in the random seed, and therefore for the initial ranking ties resolution.

## create_list_of_group_specific_markers.r
- Loads stability test runs, tallies genes occurring in at least 1/4 of runs as “stable,” and writes stable marker IDs.
- Fixes a few gene aliases to gencode hs19 names, splits markers into cancer vs normal, and saves hs19 gene lists.

## genes_that_are_markers_null_hypothesis.r
- Loads expression matrix and previously computed friends results from the first fun..
- Permutes sample labels $10^6$ times in parallel to build a null distribution of how many group markers would appear by chance, then saves the distribution table.

## plot_distribs.r
- Reads single-run markers, null distribution, and stability-test runs.
- Plots histogram PDFs comparing null vs observed marker counts and gene occurrence frequencies across runs, saving p_a_a6.pdf and p_b_a6.pdf.

## save_markers_and_friends_for_group_specific__markers.r
- Loads friends results and the stable marker list, filters to friends of stable markers, and exports them to stable_markers_with_friends_in_The_Run.csv for the notebook run below.

## jaccard_weighted.ipynb
- Notebook calculates the friend set - based weighted Jaccard distance for each pair of markers, and shows the heatmap and runs the clustering using the distance.
