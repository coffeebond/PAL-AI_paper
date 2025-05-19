library(tidyverse)
library(Biostrings)
source('helper.R')

###----------------------------------------------------------------------------------------------------------
###--- Fig. 1e ---###
##-- prediction versus measured scatter plot 
# Training data: Frog endogenous mRNA tail-length change during oocyte maturation
# Predictions: same as the training set, cross validation

pred_df <- read_delim('../Data/INN_model_predictions/XL_oo_pp_7h_tail_length_change_train_test_CV_prediction_results.txt')
res <- prediction_w_sele_scatter_plot(pred_df, xy_breaks = seq(-60,180,30), gene_id = 'id', fn = 'XL_oocyte_INN_pred_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
###--- Fig. 3c ---###
##-- prediction versus measured scatter plot
# Training data: Frog endogenous mRNA tail-length change during oocyte maturation
# Predictions: N60(LC)-PASmos library tail-length change during frog oocyte maturation
df_pred <- read_delim('../Data/INN_model_predictions/XL_INN_model_predicts_N60_PASmos_library.txt')
res <- prediction_w_sele_scatter_plot(df_pred, gene_id = 'id', xy_breaks = seq(-40,200,40), fn = 'XL_oocyte_INN_pred_N60_PASmos_pp_7h_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
###--- Fig. 3d ---###
##-- prediction versus measured scatter plot
# Training data: N60(LC)-PASmos library tail-length change during frog oocyte maturation
# Predictions: same as the training set, cross validation
df_pred <- read_delim('../Data/INN_model_predictions/XL_oo_N60_PASmos_pp_7h_tail_length_change_train_test_CV_prediction_results.txt')
res <- prediction_w_sele_scatter_plot(df_pred, gene_id = 'id', xy_breaks = seq(-40,200,40), fn = 'XL_oocyte_N60_PASmos_INN_pred_N60_PASmos_pp_7h_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
##-- prediction versus measured scatter plot
###--- Supplementary Fig. 4b ---###
###--- Supplementary Fig. 4c ---###
# Training data: Frog endogenous mRNA and N60(LC)-PASmos library tail-length change during frog oocyte maturation
# Predictions: same as the training set, cross validation
df_pred <- read_delim('../Data/INN_model_predictions/XL_endogenous_mRNA_and_N60_PASmos_library_MINN_train_test_CV_prediction_results.txt')
res <- prediction_w_sele_scatter_plot(df_pred %>% filter(!grepl('^b5v5', id)), gene_id = 'id', xy_breaks = seq(-40,200,40), 
																			fn = 'XL_endogenous_mRNA_and_N60_PASmos_library_MINN_train_test_CV_prediction_endogenous_mRNA_scatter_plot', ft = 'pdf')
res <- prediction_w_sele_scatter_plot(df_pred %>% filter(grepl('^b5v5', id)), gene_id = 'id', xy_breaks = seq(-40,200,40), 
																			fn = 'XL_endogenous_mRNA_and_N60_PASmos_library_MINN_train_test_CV_prediction_N60_PASmos_library_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
##-- prediction versus measured scatter plot
###--- Supplementary Fig. 8b ---###
# Training data: Frog endogenous mRNA tail-length change during oocyte maturation
# Predictions: human endogenous mRNA tail-length change during oocyte maturation (MII over GV)
df_pred <- read_delim('../Data/INN_model_predictions/HS_oo_MII_vs_GV_tail_length_change_prediction_results.txt')
res <- prediction_w_sele_scatter_plot(df_pred, gene_id = 'id', xy_breaks = seq(-120,160,40), 
																			fn = 'XL_oocyte_INN_pred_HS_endogenous_mRNA_M2_vs_GV_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
##-- prediction versus measured scatter plot
###--- Supplementary Fig. 8a ---###
# Training data: Frog endogenous mRNA tail-length change during oocyte maturation
# Predictions: mouse endogenous mRNA tail-length change during oocyte maturation (MII over GV)
df_pred <- read_delim('../Data/INN_model_predictions/MM_oo_MII_vs_GV_tail_length_change_prediction_results.txt')
res <- prediction_w_sele_scatter_plot(df_pred, gene_id = 'id', xy_breaks = seq(-120,160,40), 
																			fn = 'XL_oocyte_INN_pred_MM_endogenous_mRNA_M2_vs_GV_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
##-- prediction versus measured scatter plot
# Training data: Frog, mouse, human endogenous mRNA tail-length change during oocyte maturation
# Predictions: same as the training set, cross validation
df_pred <- read_delim('../Data/INN_model_predictions/Datasets_XL_HS_MM_wCDS_train_test_CV_prediction_results.txt')

###--- Supplementary Fig. 8c ---###
res <- prediction_w_sele_scatter_plot(df_pred %>% filter(data_group == 0), gene_id = 'id', xy_breaks = seq(-60,180,30), 
																			fn = 'MINN_pred_XL_endogenous_mRNA_scatter_plot', ft = 'pdf')
###--- Fig. 6b ---###
res <- prediction_w_sele_scatter_plot(df_pred %>% filter(data_group == 1), gene_id = 'id', xy_breaks = seq(-120,160,40), 
																			fn = 'MINN_pred_HS_endogenous_mRNA_scatter_plot', ft = 'pdf')
###--- Fig. 6a ---###
res <- prediction_w_sele_scatter_plot(df_pred %>% filter(data_group == 2), gene_id = 'id', xy_breaks = seq(-120,160,40), 
																			fn = 'MINN_pred_MM_endogenous_mRNA_scatter_plot', ft = 'pdf')

###----------------------------------------------------------------------------------------------------------
##-- scatter plot comparing predicted tail-length changes and measured translation changes of reporter mRNA in Xiong et al. 2022 (PMID: 35697785)
###--- Fig. 6c ---###
tl_pred <- read_delim('../Data/Xiong_2022/Datasets_Xiong_2022_Fig4g_prediction_results.txt')
te_change <- read_delim('../Data/Xiong_2022/Xiong_2022_Fig4g_results_24h.txt')

df_all <- left_join(te_change, tl_pred, by = c('gene_id' = 'id'))
res <- scatter_plot(df_all, col_x = 'te_change_avg', col_y = 'y_pred',
										xlabel = 'Measured translational change (log2)', ylabel = 'Predicted tail-length change',
										xbreaks = seq(0,3,0.5), ybreaks = seq(-20, 40, 10),
										fn = 'Xiong_2022_predicted_tail_length_vs_TE_scatter_plot')



