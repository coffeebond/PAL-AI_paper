library(tidyverse)
library(parallel)
library(RcppAlgos)
source('helper.R')

# get prediction results for all variants
input_folder <- '/lab/solexa_bartel/coffee/Sequencing/Gnomad/PolyA/Gnomad_v4_variants_TL_prediction_CV_MINN_XL_HS_MM_L_2000_CDS/'
pred_files <- list.files(path = input_folder, full.names = TRUE, pattern = ".*_INN_predictions.txt")
pred_df <- bind_rows(map(pred_files, read_delim, delim = '\t'), .id = 'file_id')

# filter variants
pred_df_filtered <- pred_df %>% 
	mutate(file_id = NULL, var_id = paste0(chromosome, '_', position, '_', ref, '_', alt)) %>%
	
		# only keep pA sites that are either not "uncertain" or "in_PolyA_DB"
		filter((!grepl('uncertain', ref_pa)) | grepl('in_PolyA_DB', ref_pa)) %>%
		
			# For each variant, only keep one pA_site with the largest absolute tail-length difference
			group_by(var_id) %>% filter(abs(diff_tl) == max(abs(diff_tl))) %>% ungroup(.)


#----------------------------------------
# scatter plot of tail-length difference for all variants
###--- Supplementary Fig. 9c ---### 
dis_cutoff = 199
af_cutoff = 0.001
df.plot <- pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>% # filter(pas != 0) %>%
	mutate(group = case_when(ac == 1 ~ 'g1',
													 ac > 1 & af > af_cutoff ~ 'g2',
													 TRUE ~ 'g0')) %>%
		filter(group != 'g0') %>%
			mutate(pos = -dis2end - 1) %>% arrange(af) # dis2end 0 is pos -1

res <- tl_diff_vs_position_scatter_plot(df.plot, x_col = 'pos', y_col = 'diff_tl', fill_col = 'af', 
																				 xinterval = 40, pos_lims = c(-200, -1), ybreaks = seq(-90,90,15), ft = 'png',
																				 fn = 'Gnomad_v4_variants_tl_diff_scatter_plot')

#----------------------------------------
# compare fraction of singletons and frequent allele in each tail-length bin
###--- Fig. 7d ---###
###--- Supplementary Fig. 9d ---### need to uncomment "#filter(pas == 0) %>%" for "df.freq" and "df.sg"

dis_cutoff = 99
af_cutoff_vec <- 10^seq(-4,-1,1)
af_groups <- c('AF > 0.01%', 'AF > 0.1%', 'AF > 1%', 'AF > 10%')
diff_breaks <- c(-Inf, seq(-25,25,10), Inf)

df.freq <- lapply(1:length(af_cutoff_vec), function(x){
	pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>% 
		filter(ac > 1 & af > af_cutoff_vec[x]) %>% #filter(pas == 0) %>%
		mutate(bin = cut(diff_tl, breaks = diff_breaks, labels = paste0("(", head(diff_breaks, -1), ", ", tail(diff_breaks, -1), "]"), include.lowest = TRUE, right = FALSE)) %>%
		group_by(bin) %>% summarize(n1 = n()) %>% 
		mutate(frac = n1/sum(n1), n2 = sum(n1) - n1, label = af_groups[x])
}) %>% bind_rows(.)

df.sg <- pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>%	
	filter(ac == 1) %>% #filter(pas == 0) %>%
	mutate(bin = cut(diff_tl, breaks = diff_breaks, labels = paste0("(", head(diff_breaks, -1), ", ", tail(diff_breaks, -1), "]"), include.lowest = TRUE, right = FALSE)) %>%
	group_by(bin) %>% summarize(n1 = n()) %>% 
	mutate(frac = n1/sum(n1), n2 = sum(n1) - n1, label = 'Singleton')

# combine two dataframes
df.combined <- bind_rows(df.sg, df.freq) %>% 
	mutate(label = factor(label, levels = c('Singleton', af_groups)))

# calculate ratio of fraction and fisher test
df_ratio_pval <- lapply(af_groups, function(x){
	lapply(df.combined %>% pull(bin) %>% unique, function(y){
		tibble(
			'group_1' = 'Singleton',
			'group_2' = x,
			'bin' = y,
			'n_group_1_in_bin' = df.combined %>% filter(bin == y & label == 'Singleton') %>% pull(n1),
			'n_group_1_out_bin' = df.combined %>% filter(bin == y & label == 'Singleton') %>% pull(n2),
			'n_group_2_in_bin' = df.combined %>% filter(bin == y & label == x) %>% pull(n1),
			'n_group_2_out_bin' = df.combined %>% filter(bin == y & label == x) %>% pull(n2),
			'pval' = fisher.test(
				df.combined %>% filter(bin == y & label %in% c('Singleton', x)) %>%
					arrange(label) %>% select(n1, n2) %>% as.matrix(.) %>% t(.),
				alternative = 'g')$p.value,
			'ratio_of_frac' = df.combined %>% filter(bin == y & label %in% c('Singleton', x)) %>%
				arrange(label) %>% pull(frac) %>% {.[1]/.[2]}
		)
	}) %>% bind_rows
}) %>% bind_rows

# output the test results
write_delim(df_ratio_pval, file = paste0('Gnomad_v4_fraction_of_variants_by_allele_frequency_and_tail_length_change_bin_dis2end_', dis_cutoff + 1, '_test_results.txt'), delim ='\t')
#write_delim(df_ratio_pval, file = paste0('Gnomad_v4_fraction_of_variants_by_allele_frequency_and_tail_length_change_bin_no_pas_change_dis2end_', dis_cutoff + 1, '_test_results.txt'), delim ='\t')

# prepare data for plot
df.plot <- df_ratio_pval %>% mutate(r = 1/ratio_of_frac,
																		group = group_2,
																		pval_bin = case_when(
																			-log10(pval) > 3 ~ 3,
																			-log10(pval) <= 3 & -log10(pval) > 2 ~ 2,
																			-log10(pval) <= 2 & -log10(pval) > 1 ~ 1,
																			TRUE ~ 0
																		)
) %>%
	select(bin, group, pval, pval_bin, r) %>% 
	bind_rows(., tibble('group' = rep('Singleton', length(levels(.$bin))),
											'bin' = levels(.$bin),
											'pval' = NA,
											'pval_bin' = NA, 
											r = 1)) %>% 
	mutate(group = factor(group, levels = c('Singleton', af_groups)),
				 bin = factor(bin, levels = paste0("(", head(diff_breaks, -1), ", ", tail(diff_breaks, -1), "]")),
				 pval_bin = factor(pval_bin))

# make a barplot
fraction_barplot_by_bin(df.plot, 
												col_x = 'bin', col_y = 'r', col_dodge = 'group', col_anno = 'pval_bin',
												ybreaks = seq(0,1.4,0.1), 
												x_title = 'Tail-length difference bin (nt)', y_title = 'Relative frequency', 
												dodge_color = brewer.pal(9, 'OrRd')[c(3, 6:9)], 
												fn = paste0('Gnomad_v4_fraction_of_variants_by_allele_frequency_and_tail_length_change_bin_dis2end_', dis_cutoff + 1)
												#fn = paste0('Gnomad_v4_fraction_of_variants_by_allele_frequency_and_tail_length_change_bin_no_pas_change_dis2end_', dis_cutoff + 1)
)




