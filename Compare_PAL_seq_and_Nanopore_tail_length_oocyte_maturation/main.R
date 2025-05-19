library('tidyverse')
source('helper.R')

count_cutoff <- 50 

# Nanopore data from Lee et al. (PMID: 38306272)
ont <- read_delim('../Data/Mouse_oocytes/MM_oocyte_tail_lengths_Nanopore_Lee_2024_PMID_38306272.txt', delim = '\t') %>%
	rename_with(~ c('stage', 'gene_id', 'gene_name', 'count', 'tl_gmean', 'tl_mean', 'tl_med')) 

tl_diff_ont <- ont %>% filter(stage %in% c('GV', 'MII')) %>% select(stage, gene_id, gene_name, count, tl_med) %>%
	pivot_wider(names_from = stage, values_from = c('count', 'tl_med')) %>% filter(if_all(starts_with('count'), ~ .x >= count_cutoff)) %>%
		mutate(diff = tl_med_MII - tl_med_GV)

# HiSeq data from Xiang et al. (PMID: 38460509)
gv <- bind_rows(read_delim('../Data/Mouse_oocytes/GSM7716872_Mouse_oocyte_mRNA_GV_rep1_PAL_seq_v4_processed.txt.gz', delim = '\t'),
								read_delim('../Data/Mouse_oocytes/GSM7716873_Mouse_oocyte_mRNA_GV_rep2_PAL_seq_v4_processed.txt.gz', delim = '\t'))

# poly(A) site annotation from Xiang et al. (PMID: 38460509)
pa_site <- read_delim('../Data/Annotation/MM_pA_site_PMID_38460509_Table_S2.txt', delim = '\t')

# get primary isoform for each gene (Nanopore data don't have 3' end isoform information)
frac_cutoff <- 0.9 # fraction of the primary isoform
df_pa_by_gene <- gv %>% 
	group_by(Gene_id_or_pA_id) %>% summarize(med = median(tail_length), n =n(), .groups = 'drop') %>% 
		dplyr::rename('pA_id' = 'Gene_id_or_pA_id') %>%
			left_join(., pa_site %>% filter(confidence == 'high') %>% select(pA_id, gene_name), by = 'pA_id') %>%
				group_by(gene_name) %>% mutate(frac = n/sum(n)) %>% ungroup(.) %>%
					filter(frac >= frac_cutoff & n >= count_cutoff)

# tail-length change between GV and MII stages calculated from Xiang et al. (PMID: 38460509)
tl_diff_hiseq <- read_delim('../Data/Mouse_oocytes/MM_oocyte_maturation_MII_over_GV_tail_length_change_tag_cutoff_10.txt', delim = '\t') %>%
	filter(count_min >= count_cutoff) %>% dplyr::rename('pA_id' = 'gene_id')

# compare between nanopore data and hiseq data
df_all <- inner_join(tl_diff_hiseq %>% select(pA_id, diff), 
										 df_pa_by_gene %>% select(pA_id, gene_name), 
										 by = 'pA_id') %>%
	inner_join(., tl_diff_ont %>% select(gene_name, diff), by = 'gene_name', suffix = c('_hiseq', '_ont'))

# scatter plot comparing two datasets
###--- Supplementary Fig. 8e ---###
res <- my_scatter_plot(df_all, 
											 gene_id = 'pA_id', xaxis = 'diff_hiseq', yaxis = 'diff_ont', 
											 xy_breaks = seq(-150, 150, 50), label_sele = NULL, color_sele = NULL, 
											 title_x = 'Tail-length difference (nt, by HiSeq)', 
											 title_y = 'Tail-length difference (nt, by Nanopore)', 
											 fn = 'MM_oocyte_maturation_MII_over_GV_tail_length_difference_compare_HiSeq_and_Nanopore_scatter_plot')

# INN-predicted results
pred <- read_delim('../Data/INN_model_predictions/Datasets_XL_HS_MM_wCDS_train_test_CV_prediction_results.txt', delim = '\t') %>% 
	dplyr::rename('pA_id' = 'id') %>% filter(data_group == 2)
df.plot <- inner_join(tl_diff_ont %>% select(gene_name, diff),
											df_pa_by_gene %>% select(gene_name, pA_id), by = 'gene_name') %>%
	inner_join(., pred, by = 'pA_id')

# scatter plot comparing ONT dataset to predicted results
###--- Supplementary Fig. 8d ---###
res <- my_scatter_plot(df.plot, 
											 gene_id = 'pA_id', xaxis = 'diff', yaxis = 'y_pred', 
											 xy_breaks = seq(-150, 150, 50), label_sele = NULL, color_sele = NULL, 
											 title_x = 'Measured tail-length change (nt, by Nanopore)', 
											 title_y = 'Predicted tail-length change (nt)', 
											 fn = 'MM_oocyte_maturation_MII_over_GV_tail_length_difference_compare_Nanopore_and_prediction_scatter_plot')
