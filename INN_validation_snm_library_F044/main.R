library(tidyverse)
library(Biostrings)
library(RcppAlgos)
source('helper.R')

###----------------------------------------------------------------------------------------------------------
##-- single-nucleotide mutagenesis library F044
# prediction with the model trained on: Frog endogenous mRNA tail-length change during oocyte maturation
# run this part before making figures

# minn-predicted data (average)
df_pred <- read_delim('../Data/INN_model_predictions/XL_endogenous_mRNA_and_N60_PASmos_library_MINN_predicts_F044_results.txt') %>%
	separate(col = 'id', into = c('pa_id', 'gene_name', 'variant'), sep = "_(?![0-9])")
df_pred_wt <- df_pred %>% filter(variant == 'WT') %>% mutate(y_pred_wt = y_pred) %>% select(gene_name, y_pred_wt)
df_pred_diff <- left_join(df_pred, df_pred_wt, by = 'gene_name') %>% mutate(y_pred_diff = y_pred - y_pred_wt)

# minn-predicted data (all models, trained with XL endogenous mRNAs and N60-PASmos(LC) injected mRNAs)
df_pred_minn_all <- read_delim('../Data/INN_model_predictions/XL_endogenous_mRNA_and_N60_PASmos_library_MINN_predicts_F044_all_models.txt') %>%
	separate(col = 'id', into = c('pa_id', 'gene_name', 'variant'), sep = "_(?![0-9])")

# inn-predicted data (all models, trained with XL endogenous mRNAs)
df_pred_inn_all <- read_delim('../Data/INN_model_predictions/XL_INN_predicts_F044_all_models.txt') %>%
	separate(col = 'id', into = c('pa_id', 'gene_name', 'variant'), sep = "_(?![0-9])")

# measured data
df_measured <- read_delim('../Data/INN_model_predictions/XL_oocyte_pp_7h_library_F044_KXL001_tail_length_change_tag_cutoff_1.txt') %>%
	separate(col = 'gene_id', into = c('pa_id', 'gene_name', 'variant'), sep = "_(?![0-9])") %>%
		dplyr::rename(y_measured = diff)
df_measured_wt <- df_measured %>% filter(variant == 'WT') %>% mutate(y_measured_wt = y_measured) %>% select(gene_name, y_measured_wt)
df_measured_diff <- left_join(df_measured, df_measured_wt, by = 'gene_name') %>% mutate(y_measured_diff = y_measured - y_measured_wt)

# combine minn-predicted result (average) with measured data
df_all <- inner_join(df_pred_diff, df_measured_diff, by = c('gene_name','variant')) 

###----------------------------------------------------------------------------------------------------------
##-- barplots comparing MINN and INN predictive performance (r value with respect to measured data)

# calculate r values between predictions by each MINN model and measurements
df_r_minn <- lapply(0:9, function(i){
	df_sele <- df_pred_minn_all %>% mutate(y_pred = .data[[paste0('model_', i)]]) %>% select(gene_name, variant, y_pred)
	df_wt <- df_sele %>% filter(variant == 'WT') %>% mutate(y_pred_wt = y_pred) %>% select(gene_name, y_pred_wt)
	df_diff <- left_join(df_sele, df_wt, by = 'gene_name') %>% mutate(y_pred_diff = y_pred - y_pred_wt) %>%
		filter(variant != 'WT') 
	df_combined <- inner_join(df_diff, df_measured_diff, by = c('gene_name', 'variant')) %>% 
		group_by(gene_name) %>%
			summarize(rp = round(cor(y_measured_diff, y_pred_diff, method = 'p'), 2), 
								rs = round(cor(y_measured_diff, y_pred_diff, method = 's'), 2), .groups = 'drop') %>%
				pivot_longer(cols = -gene_name, names_to = 'r_type', values_to = 'r_value') %>%
					mutate(model_idx = i)
}) %>% bind_rows() %>% mutate(model = 'minn')

# calculate r values between predictions by each INN model and measurements
df_r_inn <- lapply(0:9, function(i){
	df_sele <- df_pred_inn_all %>% mutate(y_pred = .data[[paste0('model_', i)]]) %>% select(gene_name, variant, y_pred)
	df_wt <- df_sele %>% filter(variant == 'WT') %>% mutate(y_pred_wt = y_pred) %>% select(gene_name, y_pred_wt)
	df_diff <- left_join(df_sele, df_wt, by = 'gene_name') %>% mutate(y_pred_diff = y_pred - y_pred_wt) %>%
		filter(variant != 'WT') 
	df_combined <- inner_join(df_diff, df_measured_diff, by = c('gene_name', 'variant')) %>% 
		group_by(gene_name) %>%
		summarize(rp = round(cor(y_measured_diff, y_pred_diff, method = 'p'), 2), 
							rs = round(cor(y_measured_diff, y_pred_diff, method = 's'), 2), .groups = 'drop') %>%
		pivot_longer(cols = -gene_name, names_to = 'r_type', values_to = 'r_value') %>%
		mutate(model_idx = i)
}) %>% bind_rows() %>% mutate(model = 'inn')

# combine two dataframes
df.plot <- bind_rows(df_r_inn, df_r_minn)
gene.order <- df.plot %>% group_by(gene_name, r_type, model) %>% 
	summarize(r_avg = mean(r_value), .groups = 'drop') %>% filter(r_type == 'rp' & model == 'minn') %>%
		arrange(desc(r_avg)) %>% pull(gene_name)
df.plot <- df.plot %>% mutate(gene_name = factor(gene_name, levels = gene.order))

# t-test comparing minn and inn
df_t_test <- lapply(c('rp', 'rs'), function(r_sele){
	lapply(gene.order, function(gene){
		temp <- df.plot %>% filter(r_type == r_sele & gene_name == gene)
		tibble(
			gene_name = gene,
			n_minn = temp %>% filter(model == 'minn') %>% nrow,
			n_inn = temp %>% filter(model == 'minn') %>% nrow,
			pval = t.test(x = temp %>% filter(model == 'minn') %>% pull(r_value), 
										y = temp %>% filter(model == 'inn') %>% pull(r_value), 
										alternative = 'greater')$p.value
		)
	}) %>% bind_rows() %>% mutate(r_type = r_sele)
}) %>% bind_rows()

# output test results
write_delim(df_t_test, file = 'XL_oocyte_F044_pearson_R_values_MINN_vs_INN_test_results.txt', delim = '\t')

# pearson R barplot with errorbar and jitter points
###--- Fig. 4b ---###
multi_sample_r_value_jitter_boxplot(df.plot %>% filter(r_type == 'rp') %>% mutate(model = factor(model, levels = c('inn', 'minn'))),
																		x = 'gene_name', y = 'r_value', dodge = 'model', dodge_width = 0.7, 
																		mycolor = brewer.pal(8, 'Set2')[c(1,2)], #brewer.pal(9, 'Set1')[c(2,5)], 
																		ybreaks = seq(0, 1, 0.1), dodge_label = c('INN', 'MINN'),
																		xlab = NULL, ylab = 'Pearson R', ft = 'pdf', fn = 'XL_oocyte_F044_pearson_R_values_MINN_vs_INN_boxplot')

# spearman R barplot with errorbar and jitter points
multi_sample_r_value_jitter_boxplot(df.plot %>% filter(r_type == 'rs') %>% mutate(model = factor(model, levels = c('inn', 'minn'))),
																		x = 'gene_name', y = 'r_value', dodge = 'model', dodge_width = 0.7,  
																		mycolor = brewer.pal(8, 'Set2')[c(1,2)], #brewer.pal(9, 'Set1')[c(2,5)], 
																		ybreaks = seq(0, 0.8, 0.1), dodge_label = c('INN', 'MINN'),
																		xlab = NULL, ylab = 'Spearman R', ft = 'pdf', fn = 'XL_oocyte_F044_spearman_R_values_MINN_vs_INN_boxplot')

###----------------------------------------------------------------------------------------------------------
##-- Bar plot showing Pearson and Spearman R for each individual mRNA
df.plot <- df_all %>% filter(variant != 'WT') %>% group_by(gene_name) %>%
	summarize(rp = round(cor(y_measured_diff, y_pred_diff, method = 'p'), 2), 
						rs = round(cor(y_measured_diff, y_pred_diff, method = 's'), 2), .groups = 'drop') %>%
	pivot_longer(cols = -gene_name, names_to = 'r_type', values_to = 'r_value')
gene.order <- df.plot %>% filter(r_type == 'rp') %>% arrange(desc(r_value)) %>% pull(gene_name)
df.plot <- df.plot %>% mutate(gene_name = factor(gene_name, levels = gene.order))
multi_sample_r_value_barplot(df.plot, 
														 x = 'gene_name', y = 'r_value', dodge = 'r_type', 
														 mycolor = brewer.pal(8, 'Set2')[c(1,2)], #brewer.pal(9, 'Set1')[c(2,5)], 
														 ybreaks = seq(0, 0.9, 0.1), dodge_label = c('Pearson R', 'Spearman R'),
														 xlab = NULL, ylab = NULL, ft = 'pdf', fn = 'XL_oocyte_F044_R_values_predicted_vs_measured')

##-- bar plot with jitter showing all-minn-model generated Pearson and Spearman R values for each individual mRNA 
gene.order <- df_r_minn %>% filter(r_type == 'rp') %>% group_by(gene_name) %>% 
	summarize(avg = mean(r_value)) %>% arrange(desc(avg)) %>% pull(gene_name) 
multi_sample_r_value_jitter_barplot(df_r_minn %>% mutate(r_type = factor(r_type, levels = c('rp', 'rs')),
																												 gene_name = factor(gene_name, levels = gene.order)),
																		x = 'gene_name', y = 'r_value', dodge = 'r_type', 
																		mycolor = brewer.pal(8, 'Set2')[c(1,2)], #brewer.pal(9, 'Set1')[c(2,5)], 
																		ybreaks = seq(0, 0.9, 0.1), dodge_label = c('Rp', 'Rs'),
																		xlab = NULL, ylab = 'Correlation coefficient', ft = 'pdf', fn = 'XL_oocyte_F044_R_values_MINN_predicted_vs_measured_all_models')

###----------------------------------------------------------------------------------------------------------
##-- scatter plot comparing mean predicted and measured differences in tail-length change of PAS and CPE mutants
###--- Fig. 4c ---###

fa <- readDNAStringSet('../Data/Sequences/XL_oocyte_ISM_reporter_library_original_mRNAs.fa')
fa_gene_name <- str_split_i(names(fa), '__', 2)
motifs <- c("TTTTA", "AATAAA")

res <- lapply(motifs, function(motif){
	lapply(1:length(fa), function(i){
		motif_start_vec <- start(vmatchPattern(motif, fa))[[i]]
		lapply(1:length(motif_start_vec), function(j){
			if (length(motif_start_vec[j]) > 0){
				df_all %>% filter(variant != 'WT') %>% 
					filter(gene_name == fa_gene_name[i]) %>% # select gene
					extract(variant, into = c("ref", "pos", "alt"), regex = "([A-Z])([0-9]+)([A-Z])", remove = FALSE) %>%
					mutate(pos = as.integer(pos)) %>% # convert to integer
					mutate(rel_pos = pos - motif_start_vec[j]) %>% # relative position to the motif
					filter(rel_pos >= 0 & rel_pos < nchar(motif)) %>% # retain positions overlapping with the motif
					mutate(mt_motif = { # mutate the motif 
						temp <- rep(motif, nrow(.))
						substr(temp, rel_pos + 1, rel_pos + 1) <- alt
						temp
					})
			}
		}) %>% bind_rows()
	}) %>% bind_rows() %>% mutate(wt_motif = motif)
}) %>% bind_rows() 

# calculate the statistics
df.plot <- res %>% select(gene_name, variant, y_pred_diff, y_measured_diff, mt_motif, wt_motif) %>%
	group_by(mt_motif, wt_motif) %>% summarize(pred_avg = mean(y_pred_diff),
																						 pred_std = sd(y_pred_diff),
																						 pred_n = n(),
																						 meas_avg = mean(y_measured_diff),
																						 meas_std = sd(y_measured_diff), 
																						 meas_n = n(), .groups = 'drop') %>%
		mutate(pred_sem = pred_std/(pred_n ** 0.5), meas_sem = meas_std / (meas_n ** 0.5)) %>%
			mutate(mt_motif = factor(mt_motif),
						 wt_motif = factor(wt_motif, levels = motifs))

# make a scatter plot
plot_scatter_w_sem(df.plot,
									 col_x = 'meas_avg', col_y = 'pred_avg', 
									 col_x_sem = 'meas_sem', col_y_sem = 'pred_sem',
									 col_fill = 'wt_motif', col_label = 'mt_motif',
									 label_sele = c('ATTAAA', 'TATAAA', 'CATAAA', 'AGTAAA', 'TTTTT', 'TTTAA'),
									 stat_legend_position = 'bottomright', legend_position = c(0.1,0.9),
									 xbreaks = seq(-70,0,10), ybreaks = NULL, xy_equal = TRUE,
									 labs = c('Measured poly(A) tail-length change (nt)', 'Predicted poly(A) tail-length change (nt)'),
									 fn = 'XL_oocyte_F044_prog_7h_predicted_vs_measured_tail_length_change_average_select_motifs_scatter_plot')


###----------------------------------------------------------------------------------------------------------
# individual mRNAs
###--- Fig. 4d ---###
###--- Fig. 4e ---###
###--- Fig. 5b ---###
###--- Fig. 5c ---###
###--- Supplementary Fig. 5a ---###
###--- Supplementary Fig. 5b ---###
###--- Supplementary Fig. 5c ---###
###--- Supplementary Fig. 6b ---###
###--- Supplementary Fig. 6e ---###
###--- Supplementary Fig. 7a ---###

sele.vec <- df_all %>% arrange(gene_name) %>% pull(gene_name) %>% unique(.)
breaks.lst <- list(seq(-120, 80, 40), seq(-10, 15, 5), seq(-100, 100, 25), seq(-90, 60, 30), seq(-20, 15, 5),
									 seq(-100, 50, 25), seq(-40, 80, 20), seq(-100, 50, 25), seq(-60, 30, 15), seq(-60, 30, 15))

out_folder <- './Scatter_plots/'
ifelse(!dir.exists(file.path(out_folder)), dir.create(file.path(out_folder)), 'Directory exists!')
res <- lapply(1:length(sele.vec), function(x){
	plot_scatter_w_side_decor(df_in = df_all %>% filter(grepl(sele.vec[x], gene_name)) %>% filter(variant != 'WT'), 
														col_x = 'y_measured_diff', col_y = 'y_pred_diff', col_label = 'variant',
														stat_legend_position = 'bottomright',
														xbreaks = breaks.lst[[x]], ybreaks = NULL, xy_equal = TRUE,
														labs = c('Measured poly(A) tail-length change (nt)', 'Predicted poly(A) tail-length change (nt)'), 
														fn = paste0(out_folder, 'XL_oocyte_F044_prog_7h_', sele.vec[x], '_predicted_vs_measured_tail_length_change_scatter_plot'))
})

out_folder <- './Scatter_plots_no_side_decor/'
ifelse(!dir.exists(file.path(out_folder)), dir.create(file.path(out_folder)), 'Directory exists!')
res <- lapply(1:length(sele.vec), function(x){
	plot_scatter(df_in = df_all %>% filter(grepl(sele.vec[x], gene_name)) %>% filter(variant != 'WT'),
							 col_x = 'y_measured_diff', col_y = 'y_pred_diff', col_label = 'variant',
							 stat_legend_position = 'bottomright',
							 xbreaks = breaks.lst[[x]], ybreaks = NULL, xy_equal = TRUE,
							 labs = c('Measured poly(A) tail-length change (nt)', 'Predicted poly(A) tail-length change (nt)'),
							 fn = paste0(out_folder, 'XL_oocyte_F044_prog_7h_', sele.vec[x], '_predicted_vs_measured_tail_length_change_scatter_plot'))
})

###----------------------------------------------------------------------------------------------------------
##-- heatmaps for tail-length change of all mutants (predicted and measured)
###--- Fig. 4f ---###
###--- Fig. 4g ---###
###--- Fig. 5a ---###
###--- Fig. 5d ---###
###--- Supplementary Fig. 5d ---###
###--- Supplementary Fig. 5e ---###
###--- Supplementary Fig. 5f ---###
###--- Supplementary Fig. 6a ---###
###--- Supplementary Fig. 6d ---###
###--- Supplementary Fig. 7b ---###

# get all mRNA names
gene_vec <- df_all %>% arrange(gene_name) %>% pull(gene_name) %>% unique(.)
out_folder = './Heatmaps/'
ifelse(!dir.exists(file.path(out_folder)), dir.create(file.path(out_folder)), 'Directory exists!')

lapply(gene_vec, function(gene_sele){
	# convert the dataframe to the compatible matrix format
	mat.lst <- lapply(c('y_measured_diff', 'y_pred_diff'), function(x){
		# get the wild-type sequence
		df_wt_seq <- df_all %>% filter(gene_name == gene_sele & variant != 'WT') %>% mutate(variant) %>% 
			separate(., col = 'variant', into = c('wt_nt', 'pos', 'mt_nt'), sep = "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)") %>%
				select(wt_nt, pos) %>% unique %>% mutate(pos = as.numeric(.$pos)) %>% arrange(pos) 
		
		# convert the dataframe to matrix
		mat <- df_all %>% filter(gene_name == gene_sele & variant != 'WT') %>% 
			mutate(tl_val = !!as.symbol(x)) %>% select(variant, tl_val) %>%
				separate(., col = 'variant', into = c('wt_nt', 'pos', 'mt_nt'), sep = "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)") %>%
					mutate(pos = as.numeric(.$pos)) %>% arrange(pos) %>% select(-wt_nt) %>%
						pivot_wider(names_from = pos, values_from = tl_val) %>% 
							arrange(mt_nt) %>% select(-mt_nt) %>% mutate_all(~replace_na(., 0)) %>% as.matrix(.)
		rownames(mat) <- c('A', 'C', 'G', 'U')
		colnames(mat) <- df_wt_seq %>% pull(wt_nt) %>% gsub('T', 'U', .)
		return(mat)
	})
	
	# plot the complex heatmaps
	max_min <- max(sapply(mat.lst, function(x){max(max(x), abs(min(x)))}))
	p1 <- select_gene_ism_score_heatmap_w_logo(mat.lst[[2]], utr_label = NULL, pos_range = c(-100, -1), fn = paste0(out_folder, gene_sele, '_predicted'), offset = 0, 
																						 xbreak = 10, g_fsize = 5, block_size = 0.07, nt_label_scaler = 1, line_w = 0.75/.pt, max_min = max_min)
	p2 <- select_gene_ism_score_heatmap_w_logo(mat.lst[[1]], utr_label = NULL, pos_range = c(-100, -1), fn = paste0(out_folder, gene_sele, '_measured'), offset = 0, 
																						 xbreak = 10, g_fsize = 5, block_size = 0.07, nt_label_scaler = 1, line_w = 0.75/.pt, max_min = max_min)
	CairoPDF(file=paste0(out_folder, gene_sele, '_logo_lineplot_heatmaps_combined.pdf'), width=7.54, height=2.2)
	plots <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,1), align = 'v', axis = 'lr')
	print(plots)
	dev.off()
})




