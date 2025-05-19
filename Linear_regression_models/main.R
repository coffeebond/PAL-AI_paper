library('tidyverse')
source('helper.R')
library('matrixStats')

###----------------------------------------------------------------------------------------------------------
##-- violin plot for distributions of r
###--- Supplementary Fig. 1a ---###
info.r <- read.table('../Data/Linear_regression/Info_linear_regression_models_stats.txt', stringsAsFactors = F, header = T)
r_lst <- list()
for (i in 1:nrow(info.r)){
	f_name <- paste0(info.r[i,'Folder'], info.r[i, 'File'])
	temp <- read.table(f_name, stringsAsFactors = F, header = T)
	r_lst[[info.r[i, 'Name']]] <- temp %>% mutate(group = info.r[i, 'Name'])
}
r_df <- r_lst %>% bind_rows()

r_df <- r_df %>% mutate(group = factor(group, levels = unique(r_df$group)))
res <- r_distribution_violin_plot(df_in = r_df, 
																	groups = levels(r_df$group),
																	ybreaks = seq(0.4,0.7,0.1), group_colors = rep(brewer.pal(12, 'Set3')[1], length(levels(r_df$group))),
																	fn = 'XL_oocyte_maturation_linear_regression_model_prediction_r_violin_plot')


# examine p values of t-tests
print(res, n = Inf)

##-- heatmap for p values 
pval_heatmap(res, col_x = 'group_x', col_y = 'group_y', col_val = 'pval',
						 fn = 'XL_oocyte_maturation_linear_regression_models_r_value_comparison_pval_heatmap')

###----------------------------------------------------------------------------------------------------------
###--- Fig. 1c ---###
##-- scatter plot comparing between measured and predicted results
pred_df <- read_delim('../Data/Linear_regression/XL_oocyte_pp_7h_tail_length_change_tag_cutoff_10_train_test_CV_prediction_results.txt', delim = '\t')
res <- prediction_w_sele_scatter_plot(pred_df, 
																			xy_breaks = seq(-60,180,30), 
																			xaxis = 'y', yaxis = 'y_pred_avg',
																			gene_id = 'id', 
																			fn = 'XL_oocyte_tl_change_linear_regression_pred_scatter_plot', ft = 'pdf')

##----------------------------------------------------------------------------------------------------------
##-- bar plot for coefficients of top metrics
coef_df <- read_csv('../Data/Linear_regression/Linear_Regression_coefficient_XL_oocyte_pp_7h_tail_length_change_tag_cutoff_10_by_CV.csv') %>%
	mutate(avg = select(., -feature) %>% rowMeans,
				 std = select(., -feature) %>% as.matrix(.) %>% rowSds) %>% 
		select(feature, avg, std) %>%
			separate(feature, c('feature', 'kmer'), sep = '_')

###--- Supplementary Fig. 1b ---###
res <- coefficient_bar_plot(coef_df %>% filter(feature == 'num'), x = 'kmer', y = 'avg', y_std = 'std', 
														ybreaks = seq(-1.2, 1.2, 0.4),
														fn = 'XL_oocyte_tl_change_linear_regression_top_20_kmer_by_num_coefficients_bar_plot')
###--- Supplementary Fig. 1c ---###
res <- coefficient_bar_plot(coef_df %>% filter(feature == 'pos'), x = 'kmer', y = 'avg', y_std = 'std', 
														ybreaks = seq(-1.2, 1.2, 0.4),
														fn = 'XL_oocyte_tl_change_linear_regression_top_20_kmer_by_pos_coefficients_bar_plot')

