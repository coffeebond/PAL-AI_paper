library('tidyverse')
source('helper.R')

###----------------------------------------------------------------------------------------------------------
##-- Load data
info.r <- read.table('../Data/INN_train_test/Info_INN_models_stats.txt', stringsAsFactors = F, header = T)
r_lst <- list()
for (i in 1:nrow(info.r)){
	f_name <- paste0(info.r[i,'Folder'], info.r[i, 'File'])
	temp <- read.table(f_name, stringsAsFactors = F, header = T)
	r_lst[[i]] <- temp %>% 
		mutate(group = info.r[i, 'Name'],
					 type = info.r[i, 'Type'])
}
r_df <- r_lst %>% bind_rows()

###----------------------------------------------------------------------------------------------------------
###--- Fig. 1f ---###
###--- Fig. 1g ---###
##-- violin plot for distributions of r, compare all INN models

df.plot <- r_df %>% filter(type == 'inn') %>%
	mutate(group = factor(group, levels = unique(.$group)))
res <- r_distribution_violin_plot(df_in = df.plot, 
																	groups = levels(df.plot$group),
																	ybreaks = seq(0.3,0.9,0.1), group_colors = rep(brewer.pal(12, 'Set3')[1], length(levels(df.plot$group))),
																	fn = 'XL_oocyte_maturation_INN_prediction_r_violin_plot')

# examine p values of t-tests
print(res, n = Inf)

##-- heatmap for p values 
pval_heatmap(res, col_x = 'group_x', col_y = 'group_y', col_val = 'pval',
						 fn = 'XL_oocyte_maturation_INN_models_r_value_comparion_pval_heatmap')

###----------------------------------------------------------------------------------------------------------
##-- violin plot for distributions of r, compare ResNet models and corresponding INN models
###--- Supplementary Fig. 2e ---###

group_sele <- r_df %>% filter(type == 'resnet') %>% pull(group) %>% unique
df.plot <- r_df %>% filter(group %in% group_sele)

# t test 
pval <- t.test(x = df.plot %>% filter(type == 'inn') %>% pull(r_pearson),
							 y = df.plot %>% filter(type == 'resnet') %>% pull(r_pearson),
							 alternative = 'greater')$p.value

# violin plot
violin_plot(df_in = df.plot, 
						col_x = 'type', col_y = 'r_pearson', 
						label_x = NULL,
						ybreaks = seq(0.5,0.9,0.1),
						fill_color = brewer.pal(8, 'Set2')[c(1,2)],
						fn = 'XL_oocyte_maturation_INN_vs_ResNet_model_prediction_r_violin_plot')

	
