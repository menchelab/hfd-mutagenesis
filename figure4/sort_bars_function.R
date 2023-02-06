#this function takes a refit object and sample names and
#returns a table that is sorted by total number of mutations per sample
#so that the plot can be shown with samples in decreasing order

sort_bars <- function(refit_object, names_of_samples){
#need to generate dummy plot
p_tmp = plot_contribution(refit_object,
  coord_flip = FALSE,
  mode = "absolute", palette = mm_colors
)
#pull out total contribution data to reorder bars
bar_data = p_tmp$data %>% group_by(Sample) %>% mutate(sum(Contribution))
#make a dummy sort frame to find the correct order of bars
sorting_data = bar_data$`sum(Contribution)` %>% unique()
sorting_frame = cbind(names_of_samples, sorting_data)
sorted = data.frame(sorting_frame[order(sorting_data, decreasing = TRUE), ])

#pull out dataframe to be plotted
plotting_data <- data.frame(refit_object)
colnames(plotting_data) <- names_of_samples
#resort dataframe to be plotted
plotting_data_sorted <- setcolorder(plotting_data, sorted[,1])

return(plotting_data_sorted)
}
