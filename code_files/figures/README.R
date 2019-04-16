source("figures.R")

# To make a plot of age ranges
outfile = "path/to/plot"
pbdbfile = "path/to/ranges/RData" #ranges should be stored as a table with minimum (1st col) and maximum (2nd col) age 
figure_ranges(pbdbfile, outfile)

# To make a plot comparing relative errors and coverage on simulated data
summary_file = "file/with/estimate/error/and/coverage" # as produced by mammals_summary
plotfolder = "path/to/figures/folder"
plot_comparison(summary_file, plotfolder)

# To make a plot comparing estimates of FBD parameters on cetaceans
outfolder = "path/to/folder/with/logfiles"
plotfile = "path/to/output/figure"
cetaceans_plot(outfolder, plotfile)

# To make a plot comparing posteriors on one node
outfolder = "path/to/folder/with/logfiles"
string = "name_of_column_to_use" #should be mrca.date.backward. + node/clade name
plotfile = "path/to/output/figure"
plot_node_comparison(outfolder, string, plotfile)