library(tidyverse)
library(parallel)
library(RcppAlgos)
source('helper.R')

# get prediction results for all variants
input_folder <- '/lab/solexa_bartel/coffee/Sequencing/AllofUs/Exome_v8_3UTR_variants/All_of_us_predictions/CV_MINN_XL_HS_MM_L_2000_CDS/'
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
###--- Fig. 7b ---###
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
																				fn = 'All_of_Us_v8_variants_tl_diff_scatter_plot')

#----------------------------------------
# compare fraction of singletons and frequent allele in each tail-length bin
###--- Fig. 7c ---###
###--- Supplementary Fig. 9b ---### need to uncomment "#filter(pas == 0) %>%" for "df.freq" and "df.sg"
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
	filter(ac == 1) %>% # filter(pas == 0) %>%
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
			'pval' = fisher.test(
				df.combined %>% filter(bin == y & label %in% c('Singleton', x)) %>%
					arrange(label) %>% select(n1, n2) %>% as.matrix(.) %>% t(.),
				alternative = 'g')$p.value,
			'ratio_of_frac' = df.combined %>% filter(bin == y & label %in% c('Singleton', x)) %>%
				arrange(label) %>% pull(frac) %>% {.[1]/.[2]}
		)
	}) %>% bind_rows
}) %>% bind_rows

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
												fn = paste0('All_of_Us_fraction_of_variants_by_allele_frequency_and_tail_length_change_bin_dis2end_', dis_cutoff + 1)
												#fn = paste0('All_of_Us_fraction_of_variants_by_allele_frequency_and_tail_length_change_bin_no_pas_change_dis2end_', dis_cutoff + 1)
)






#----------------------------------------
# CDF plot of predicted tail-length differences, comparing between singletons and common alleles
# including PAS
dis_cutoff = 100
af_cutoff = 0.001
df.plot <- pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>%
	mutate(group = case_when(ac == 1 ~ 'g1',
													 ac > 1 & af > af_cutoff ~ 'g2',
													 TRUE ~ 'g0')) %>%
	filter(group != 'g0')
res <- tl_diff_cdf_plot(df_in = df.plot, groups = c('g1', 'g2'), group_col = 'group', val_col = 'diff_tl', 
												group_label = c('Singleton', 'AF > 0.1%'), xbreaks = seq(-15,6,3), ybreaks = seq(0,1,0.1),
												fn = 'AllofUs_variants_predicted_tl_diff_cdf_plot')


#----------------------------------------
# Examine the impact on the p-value and the enrichment of severe outcomes
# by different region of 3'-UTRs, different AF cutoff, and different tail length difference cutoff 
# This part may take quite a while to calculate

dis_cutoff_vec <- seq(50,1000,50)
af_cutoff_vec <- 10^seq(-1,-5,-1)
tl_diff_vec <- seq(-25,-5,5)
df.combo <- expand.grid(dis = dis_cutoff_vec, af = af_cutoff_vec)

df.res <- mclapply(1:nrow(df.combo), function(x){
	temp <- pred_df_filtered %>% filter(dis2end <= df.combo[,'dis'][x]) %>%
		mutate(group = case_when(ac == 1 ~'g1',
														 ac > 1 & af >= df.combo[,'af'][x] ~ 'g2',
														 TRUE ~ 'g0'),
					 af_group = df.combo[,'af'][x],
					 dis2end_group = df.combo[,'dis'][x]) %>%
		filter(group != 'g0')
	
	# calculate p-value for wilcox test
	pval <- wilcox.test(x = temp %>% filter(group == 'g1') %>% pull(diff_tl),
											y = temp %>% filter(group == 'g2') %>% pull(diff_tl), alternative = 'less')$p.value
	
	# calculate percentage of severe outcomes
	g1_frac <- sapply(tl_diff_vec, function(y){
		temp %>% filter(group == 'g1') %>% {sum(.$diff_tl < y)/nrow(.)}
	})
	g2_frac <- sapply(tl_diff_vec, function(y){
		temp %>% filter(group == 'g2') %>% {sum(.$diff_tl < y)/nrow(.)}
	})
	fold_lst <- as.list(ifelse(g2_frac == 0, NA, g1_frac/g2_frac))
	names(fold_lst) <- paste0('fold_', tl_diff_vec)
	return(c(list('dis_cutoff' = df.combo[,'dis'][x],
								'af_cutoff' = df.combo[,'af'][x],
								'pval' = pval),
					 fold_lst))
}, mc.cores = 20) %>% bind_rows(.)

#----------------------------------------
# bar plot for relative frequency of highly disruptive variants (tail-length shortening, loss-of-function)
dis_cutoff = 100
af_cutoff_vec <- 10^seq(-4,-1,1)
af_groups <- c('AF > 0.01%', 'AF > 0.1%', 'AF > 1%', 'AF > 10%')
tl_diff_cutoff <- (-30)

df.plot <- lapply(1:length(af_cutoff_vec), function(x){
	pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>%
		filter(ac > 1 & af > af_cutoff_vec[x]) %>% 
		summarize(n1 = sum(diff_tl < tl_diff_cutoff), n = n()) %>%
		mutate(frac = n1/n, n2 = n - n1, label = af_groups[x])
}) %>% bind_rows(pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>%
								 	filter(ac == 1) %>% summarize(n1 = sum(diff_tl < tl_diff_cutoff), n = n()) %>%
								 	mutate(frac = n1/n, n2 = n - n1, label = 'Singleton'), .) %>%
	mutate(frac_nor = frac / max(.$frac), label = factor(label, levels = c('Singleton', af_groups)))

# calculate fisher test
fisher.pvals <- sapply(2:nrow(df.plot), function(x){
	fisher.test(
		df.plot %>% slice(1,x) %>% select(n1, n2) %>% as.matrix(.) %>% t(.),
		alternative = 'g')$p.value
})

fraction_barplot(df.plot, y_col = 'frac_nor', x_col = 'label', ybreaks = seq(0,1,0.2),
								 ylabel = paste0('Relative frequency of variants\n(\u0394L < \u2212', abs(tl_diff_cutoff), ' nt)'),
								 fn = paste0('AllofUs_disruptive_variants_TL_decrease_', tl_diff_cutoff, 'nt_fraction_barplot'))

#----------------------------------------
# compare fraction of singletons and frequent allele by different tail-length cutoff
dis_cutoff = 200
af_cutoff_vec <- 10^seq(-4,-1,1)
af_groups <- c('AF > 0.01%', 'AF > 0.1%', 'AF > 1%', 'AF > 10%')
diff_breaks <- list(c(-Inf, -25), c(-Inf, -20), c(-Inf, -15), c(-Inf, -10), c(10, Inf), c(15, Inf), c(20, Inf), c(25, Inf))

df.freq <- lapply(1:length(diff_breaks), function(x){
	lapply(1:length(af_cutoff_vec), function(y){
		pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>%
			filter(ac > 1 & af > af_cutoff_vec[y]) %>%
				mutate(in_tl_bin = ifelse(diff_tl >= diff_breaks[[x]][1] & diff_tl <= diff_breaks[[x]][2], 'in_tl_bin', 'out_tl_bin')) %>%
					group_by(in_tl_bin) %>% summarize(n = n(), .groups = 'drop') %>%
						pivot_wider(names_from = in_tl_bin, values_from = n) %>%
							mutate(frac = in_tl_bin/(in_tl_bin + out_tl_bin),
										 af_group = af_groups[y])
	}) %>% bind_rows(.) %>% 
		mutate(tl_group = paste0("[", diff_breaks[[x]][1], ", ", diff_breaks[[x]][2], "]"))
}) %>% bind_rows(.)

df.sg <- lapply(1:length(diff_breaks), function(x){
		pred_df_filtered %>% filter(dis2end <= dis_cutoff) %>%
			filter(ac == 1) %>%
				mutate(in_tl_bin = ifelse(diff_tl >= diff_breaks[[x]][1] & diff_tl <= diff_breaks[[x]][2], 'in_tl_bin', 'out_tl_bin')) %>%
					group_by(in_tl_bin) %>% summarize(n = n(), .groups = 'drop') %>%
						pivot_wider(names_from = in_tl_bin, values_from = n) %>%
							mutate(frac = in_tl_bin/(in_tl_bin + out_tl_bin)) %>%
								mutate(tl_group = paste0("[", diff_breaks[[x]][1], ", ", diff_breaks[[x]][2], "]"))
	}) %>% bind_rows(.) %>% mutate(af_group = 'Singleton')

df.combined <- bind_rows(df.sg, df.freq) %>% 
	mutate(label = factor(af_group, levels = c('Singleton', af_groups)))

# calculate ratio of fraction and fisher test
df_ratio_pval <- lapply(af_groups, function(x){
	lapply(df.combined %>% pull(tl_group) %>% unique, function(y){
		tibble(
			'group_1' = 'Singleton',
			'group_2' = x,
			'bin' = y,
			'pval' = fisher.test(
				df.combined %>% filter(tl_group == y & label %in% c('Singleton', x)) %>%
					arrange(label) %>% select(in_tl_bin, out_tl_bin) %>% as.matrix(.) %>% t(.),
				alternative = 'g')$p.value,
			'ratio_of_frac' = df.combined %>% filter(tl_group == y & label %in% c('Singleton', x)) %>%
				arrange(label) %>% pull(frac) %>% {.[1]/.[2]}
		)
	}) %>% bind_rows
}) %>% bind_rows

# prepare data for plot
df.plot <- df_ratio_pval %>% mutate(r = 1/ratio_of_frac,
																		group = group_2,
																		pval_bin = case_when(
																			-log10(pval) >= 3 ~ 3,
																			-log10(pval) < 3 & -log10(pval) >= 2 ~ 2,
																			-log10(pval) < 2 & -log10(pval) >= 1 ~ 1,
																			TRUE ~ 0
																		)
) %>%
	select(bin, group, pval, pval_bin, r) %>% 
		bind_rows(., tibble('group' = rep('Singleton', length(unique(.$bin))),
												'bin' = unique(.$bin),
												'pval' = NA,
												'pval_bin' = NA, 
												r = 1)) %>% 
			mutate(group = factor(group, levels = c('Singleton', af_groups)),
						 bin = factor(bin, levels = unique(bin)),
						 pval_bin = factor(pval_bin))

# make a barplot
fraction_barplot_by_bin(df.plot, 
												col_x = 'bin', col_y = 'r', col_dodge = 'group', col_anno = 'pval_bin',
												ybreaks = seq(0,1.4,0.1), 
												x_title = 'Tail-length difference bin (nt)', y_title = 'Relative frequency', 
												dodge_color = brewer.pal(9, 'OrRd')[c(3, 6:9)], 
												fn = paste0('All_of_Us_fraction_of_variants_by_allele_frequency_and_tail_length_change_cutoff_dis2end_', dis_cutoff)
)
