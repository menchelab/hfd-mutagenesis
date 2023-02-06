#Suite of Functions needed to add accureate error bars to 96-channel mutational profile 
#when aggregate data is shown

#function that pulls out plotting data from many runs an stored them in a table
pull_plot_data <- function(mutation_matrix, y_max_param=0.1){
  tmp_plot <- plot_96_profile(mutation_matrix, ymax = y_max_param)
  plot_data <- tmp_plot$data$freq
  return(plot_data)
}

#function to reformat the plotting data into appropriate columns
sort_into_columns <- function(plot_data){
  total = length(plot_data)
  iterations = total/96
  start_index=1
  stop_index=96
  df_plot_data <- data.frame(matrix(ncol=iterations, nrow=96))
  for(i in c(1:iterations)){
    df_plot_data[,i]= plot_data[start_index:stop_index]
    start_index=start_index+96
    stop_index=stop_index+96
  }
  return(df_plot_data)
}

#fuctions to calculate standard deviation and y min max of error bars based on sd
calculate_error_bar_position_min <- function(df_plot_data, mutation_matrix_to_plot, y_max_param =0.1){
  sd = apply(df_plot_data,1,sd)
  bar_position = (plot_96_profile(mutation_matrix_to_plot, ymax= y_max_param))$data$freq
  bar_min = bar_position-sd
  return(bar_min)
}

calculate_error_bar_position_max <- function(df_plot_data, mutation_matrix_to_plot, y_max_param =0.1){
  sd = apply(df_plot_data,1,sd)
  bar_position = (plot_96_profile(mutation_matrix_to_plot, ymax= y_max_param))$data$freq
  bar_max = bar_position+sd
  return(bar_max)
  
}

#function for calculating and plotting error bars
add_error_bars <- function(source_mutation_matrix, mutation_matrix_to_plot, y_max_param=0.1) {
  tmp1 = pull_plot_data(source_mutation_matrix, y_max_param)
  tmp2 = sort_into_columns(tmp1)
  tmp_min=calculate_error_bar_position_min(tmp2, mutation_matrix_to_plot, y_max_param = y_max_param)
  tmp_max=calculate_error_bar_position_max(tmp2, mutation_matrix_to_plot, y_max_param = y_max_param)
  
  ptmp <- plot_96_profile(mutation_matrix_to_plot, ymax=y_max_param)
  ptmp <- ptmp + theme(axis.text.y = element_text(size=12),axis.title.y = element_text(size = 14),
                       axis.title.x = element_text(size = 14))
  ptmp <- ptmp + geom_errorbar(aes(ymax=tmp_max, ymin=tmp_min),width=.1,
                               position=position_dodge(.9))
  ptmp
}


####################################################################
#functions for adding error bars to main indel profiles
#function that pulls out plotting data from many runs an stored them in a table
pull_plot_data_indel_main <- function(mutation_matrix){
  tmp_plot <- plot_main_indel_contexts(mutation_matrix)
  plot_data <- tmp_plot$data$count
  return(plot_data)
}

#function to reformat the plotting data into appropriate columns
sort_into_columns_general <- function(plot_data, channel_no){
  total = length(plot_data)
  iterations = total/channel_no
  start_index=1
  stop_index=channel_no
  df_plot_data <- data.frame(matrix(ncol=iterations, nrow=channel_no))
  for(i in c(1:iterations)){
    df_plot_data[,i]= plot_data[start_index:stop_index]
    start_index=start_index+channel_no
    stop_index=stop_index+channel_no
  }
  return(df_plot_data)
}

#fuctions to calculate standard deviation and y min max of error bars based on sd
calculate_error_bar_position_min_indel_main <- function(df_plot_data, mutation_matrix_to_plot){
  sd = apply(df_plot_data,1,sd)
  bar_position = (plot_main_indel_contexts(mutation_matrix_to_plot))$data$count
  bar_min = bar_position-sd
  return(bar_min)
}

calculate_error_bar_position_max_indel_main <- function(df_plot_data, mutation_matrix_to_plot){
  sd = apply(df_plot_data,1,sd)
  bar_position = (plot_main_indel_contexts(mutation_matrix_to_plot))$data$count
  bar_max = bar_position+sd
  return(bar_max)

}

#function for calculating and plotting error bars
add_error_bars_indel_main <- function(source_mutation_matrix, mutation_matrix_to_plot, channel_no) {
  tmp1 = pull_plot_data_indel_main(source_mutation_matrix)
  tmp2 = sort_into_columns_general(tmp1, channel_no = channel_no)
  tmp_min=calculate_error_bar_position_min_indel_main(tmp2, mutation_matrix_to_plot)
  tmp_max=calculate_error_bar_position_max_indel_main(tmp2, mutation_matrix_to_plot)

  ptmp <- plot_main_indel_contexts(mutation_matrix_to_plot)
  ptmp <- ptmp + theme(axis.text.y = element_text(size=12),axis.title.y = element_text(size = 14),
                       axis.title.x = element_text(size = 14))
  ptmp <- ptmp + geom_errorbar(aes(ymax=tmp_max, ymin=tmp_min),width=.1,
                               position=position_dodge(.9))
  ptmp
}


###########################################################################
#rewrite error bar function for indels
#--------------------------------------------------------
#function that pulls out plotting data from many runs an stored them in a table
pull_plot_data_indel_all <- function(mutation_matrix){
  tmp_plot <- plot_indel_contexts(mutation_matrix)
  plot_data <- tmp_plot$data$count
  return(plot_data)
}

#fuctions to calculate standard deviation and y min max of error bars based on sd
calculate_error_bar_position_min_indel_all <- function(df_plot_data, mutation_matrix_to_plot){
  sd = apply(df_plot_data,1,sd)
  bar_position = (plot_indel_contexts(mutation_matrix_to_plot))$data$count
  bar_min = bar_position-sd
  return(bar_min)
}

calculate_error_bar_position_max_indel_all <- function(df_plot_data, mutation_matrix_to_plot){
  sd = apply(df_plot_data,1,sd)
  bar_position = (plot_indel_contexts(mutation_matrix_to_plot))$data$count
  bar_max = bar_position+sd
  return(bar_max)

}

#function for calculating and plotting error bars
add_error_bars_indel_all <- function(source_mutation_matrix, mutation_matrix_to_plot, channel_no) {
  tmp1 = pull_plot_data_indel_all(source_mutation_matrix)
  tmp2 = sort_into_columns_general(tmp1, channel_no = channel_no)
  tmp_min=calculate_error_bar_position_min_indel_all(tmp2, mutation_matrix_to_plot)
  tmp_max=calculate_error_bar_position_max_indel_all(tmp2, mutation_matrix_to_plot)

  ptmp <- plot_indel_contexts(mutation_matrix_to_plot)
  ptmp <- ptmp + theme(axis.text.y = element_text(size=12),axis.title.y = element_text(size = 14),
                       axis.title.x = element_text(size = 14))
  ptmp <- ptmp + geom_errorbar(aes(ymax=tmp_max, ymin=tmp_min),width=.1,
                               position=position_dodge(.9))
  ptmp
}
