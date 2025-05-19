library('tidyverse')
source('helper.R')

input_folder <- '../Data/CV_L_2000_CDS_loss_plots/'

loss_files <- list.files(path = input_folder, full.names = TRUE, pattern = ".*_loss_.*")
metric_files <- list.files(path = input_folder, full.names = TRUE, pattern = ".*_metric_.*")

loss_df <- bind_rows(map(loss_files, read_delim, delim = '\t'), .id = "file_id")
metric_df <- bind_rows(map(metric_files, read_delim, delim = '\t'), .id = "file_id")

# loss (mse) vs epoch line plot
###--- Supplementary Fig. 2c ---###
df.plot <- loss_df %>% dplyr::rename('Validation' = 'Test') %>% 
	pivot_longer(cols = c('Train', 'Validation'), names_to = 'group', values_to = 'value') %>%
		group_by(Epochs, group) %>% summarize(avg = mean(value), std = sd(value), n =n(), .groups = 'drop') %>%
			mutate(sem = std/sqrt(n)) 

line_plot(df.plot, col_x = 'Epochs', col_y = 'avg', col_y_delt = 'sem', col_group = 'group',
					xbreaks = seq(5, 40, 5), ybreaks = seq(200, 1000, 100), y_title = 'Loss (MSE)', x_title = 'Epochs', fn = 'Loss_vs_epochs_line_plot')

# metric (mae) vs epoch line plot
###--- Supplementary Fig. 2d ---###
df.plot <- metric_df %>% dplyr::rename('Validation' = 'Test') %>% 
	pivot_longer(cols = c('Train', 'Validation'), names_to = 'group', values_to = 'value') %>%
		group_by(Epochs, group) %>% summarize(avg = mean(value), std = sd(value), n =n(), .groups = 'drop') %>%
			mutate(sem = std/sqrt(n))

line_plot(df.plot, col_x = 'Epochs', col_y = 'avg', col_y_delt = 'sem', col_group = 'group',
					xbreaks = seq(5, 40, 5), ybreaks = seq(10, 22, 2), y_title = 'Metric (MAE)', x_title = 'Epochs', fn = 'Metric_vs_epochs_line_plot')