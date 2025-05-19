library('tidyverse')
source('helper.R')

###----------------------------------------------------------------------------------------------------------
##-- 
tb <- read.csv('../Data/INN_train_test/INN_XL_L_2000_CDS_hyperparameter_serch.csv')

# Pearson R versus n_Conv1D_MaxPool_block
###--- Supplementary Fig. 2a ---###
df.plot <- tb %>% select(value, params_n_Conv1D_MaxPool_block) %>% 
	mutate(r = sqrt(-value), 
				 n_repeats = params_n_Conv1D_MaxPool_block - 1)
violin_plot(df_in = df.plot, 
						col_x = 'n_repeats', col_y = 'r',
						label_x = 'Number of repeats\n(Conv1D MaxPool block)',
						ybreaks = seq(0.6,0.9,0.1),
						col_fill = brewer.pal(9, 'GnBu')[c(6,6,6)],
						fn = 'INN_r_vs_n_Conv1D_MaxPool_block')

# Pearson R versus dropout_rate
###--- Supplementary Fig. 2a ---###
df.plot <- tb %>% select(value, params_dropout_rate) %>% mutate(r = sqrt(-value))
violin_plot(df_in = df.plot, 
						col_x = 'params_dropout_rate', col_y = 'r',
						label_x = 'Dropout rate\n',
						ybreaks = seq(0.6,0.9,0.1),
						col_fill = brewer.pal(9, 'GnBu')[c(6,6,6)],
						fn = 'INN_r_vs_n_dropout_rate')

# Pearson R versus n_Conv1D_filter
###--- Supplementary Fig. 2a ---###
df.plot <- tb %>% select(value, params_n_Conv1D_filter) %>% mutate(r = sqrt(-value))
violin_plot(df_in = df.plot, 
						col_x = 'params_n_Conv1D_filter', col_y = 'r',
						label_x = 'Number of filters\n(Conv1D)',
						ybreaks = seq(0.6,0.9,0.1),
						col_fill = brewer.pal(9, 'GnBu')[c(6,6,6,6)],
						fn = 'INN_r_vs_n_Conv1D_filter')

# Pearson R versus n_Conv1D_kernel
###--- Supplementary Fig. 2a ---###
df.plot <- tb %>% select(value, params_n_Conv1D_kernel) %>% mutate(r = sqrt(-value))
violin_plot(df_in = df.plot, 
						col_x = 'params_n_Conv1D_kernel', col_y = 'r',
						label_x = 'Kernel size\n(Conv1D)',
						ybreaks = seq(0.6,0.9,0.1),
						col_fill = brewer.pal(9, 'GnBu')[c(6,6,6,6)],
						fn = 'INN_r_vs_n_Conv1D_kernel')


# Pearson R versus n_rnn_units
###--- Supplementary Fig. 2a ---###
df.plot <- tb %>% select(value, params_n_rnn_units) %>% mutate(r = sqrt(-value))
violin_plot(df_in = df.plot, 
						col_x = 'params_n_rnn_units', col_y = 'r',
						label_x = 'Number of units\n(GRU)',
						ybreaks = seq(0.6,0.9,0.1),
						col_fill = brewer.pal(9, 'GnBu')[c(6,6,6,6)],
						fn = 'INN_r_vs_n_rnn_units')


# Pearson R versus n_dense_neuron
###--- Supplementary Fig. 2a ---###
df.plot <- tb %>% select(value, params_n_dense_neuron) %>% mutate(r = sqrt(-value))
violin_plot(df_in = df.plot, 
						col_x = 'params_n_dense_neuron', col_y = 'r',
						label_x = 'Number of units\n(Dense layer)',
						ybreaks = seq(0.6,0.9,0.1),
						col_fill = brewer.pal(9, 'GnBu')[c(6,6,6,6)],
						fn = 'INN_r_vs_n_dense_neuron')

# Pearson R versus learning_rate
###--- Supplementary Fig. 2b ---###
df.plot <- tb %>% select(value, params_learning_rate) %>% mutate(r = sqrt(-value))
scatter_plot(df_in = df.plot, 
						col_x = 'params_learning_rate', col_y = 'r',
						label_x = 'Learning rate (log10)',
						xbreaks = c(0.00001, 0.0001, 0.001, 0.01),
						ybreaks = seq(0.4,0.9,0.1),
						fn = 'INN_r_vs_learning_rate')

# Pearson R versus l1_reg
###--- Supplementary Fig. 2b ---###
df.plot <- tb %>% select(value, params_l1_reg) %>% mutate(r = sqrt(-value))
scatter_plot(df_in = df.plot, 
						 col_x = 'params_l1_reg', col_y = 'r',
						 label_x = 'L1 regularization',
						 xbreaks = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1),
						 ybreaks = seq(0.4,0.9,0.1),
						 fn = 'INN_r_vs_l1_reg')

# Pearson R versus l2_reg
###--- Supplementary Fig. 2b ---###
df.plot <- tb %>% select(value, params_l2_reg) %>% mutate(r = sqrt(-value))
scatter_plot(df_in = df.plot, 
						 col_x = 'params_l2_reg', col_y = 'r',
						 label_x = 'L2 regularization',
						 xbreaks = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1),
						 ybreaks = seq(0.4,0.9,0.1),
						 fn = 'INN_r_vs_l2_reg')


