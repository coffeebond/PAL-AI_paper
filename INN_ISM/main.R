library(tidyverse)
library(Biostrings)
library(rhdf5)
source('helper.R')

#------------------------------------------------------------------------------------------------------------------------------------
#--- examine predicted tail-length changes in ISM results (not iterative exclusion)
kmer.gain <- read_delim('../Data/INN_ISM/XL_oocyte_ISM_average_utr_len_limit_300_kmer_gain_pA_weight_avg_dis.txt', col_names= F)
colnames(kmer.gain) <- c('kmer', 'n', 'avg', 'dis')
kmer.gain <- kmer.gain %>% mutate(std = (dis/(n-1)) ** 0.5) %>% mutate(sem = std / (n ** 0.5))

kmer.loss <- read_delim('../Data/INN_ISM/XL_oocyte_ISM_average_utr_len_limit_300_kmer_loss_pA_weight_avg_dis.txt', col_names = F)
colnames(kmer.loss) <- c('kmer', 'n', 'avg', 'dis')
kmer.loss <- kmer.loss %>% mutate(std = (dis/(n-1)) ** 0.5) %>% mutate(sem = std / (n ** 0.5))

#---  scatter plot of difference in mean pA length versus mean count number 
###--- Supplementary Fig. 3a ---###
res <- tl_change_vs_count_splot(df_in = kmer.gain, kmer_len = 6, count_cutoff = 1, 
																label_sele = list(c('TTTTAA','ATTTTA', 'TTTTTA', 'TTTTAT', 'GTTTTA', 'CTTTTA'), c('AATAAA')), 
																ybreaks = seq(-4,8,2),
																fn = 'XL_oocyte_ISM_6mer_gain_pA_difference_vs_count_splot')

###--- Fig. 2a ---###
res <- tl_change_vs_count_splot(df_in = kmer.loss, kmer_len = 6, count_cutoff = 1, 
																label_sele = list(c('TTTTAA', 'ATTTTA', 'TTTTTA', 'TTTTAT', 'GTTTTA', 'CTTTTA'), 
																									c('AATAAA', 'TAATAA', 'ATAAAG', 'CAATAA', 'AAATAA', 'ATAAAA')), 
																ybreaks = seq(-12,2,2),
																fn = 'XL_oocyte_ISM_6mer_loss_pA_difference_vs_count_splot')

#--- logo plot from 8mer, aligned with custom score, for kmer gain
# filter kmer by performing one sample t-test against global mean
kmer_len <- 8
ism.sele <- kmer.gain %>% filter(nchar(kmer) == kmer_len)
avg_mean <- ism.sele %>% pull(avg) %>% mean(.)
avg_std <- ism.sele %>% pull(avg) %>% sd(.)
ism.sele <- ism.sele %>% mutate(z_score = (avg - avg_mean) / avg_std) %>%
	mutate(p_positive = one_sample_t_test(avg, std, n, diff = avg_mean, lower.tail = F),
				 p_negative = one_sample_t_test(avg, std, n, diff = avg_mean, lower.tail = T)) %>%
	mutate(log10p = ifelse(avg > avg_mean, -log10(pmax(p_positive, .Machine$double.xmin)), -log10(pmax(p_negative, .Machine$double.xmin))))

res <- plot_top_kmer_logo(as.data.frame(ism.sele %>% filter(log10p >= 2)), 
													weight = 'avg', mode = 'positive', n_seed_max = Inf, z_cutoff = 3,
													align_score_cutoff = -0.1, score_match = 1, score_offset = -1, score_mismatch = -1.5, type = 'bits',
													fn = paste0('XL_oocyte_ISM_', kmer_len, 'mer_top_kmers_logo'))

# combined logo plot using the result from 'plot_top_kmer_logo'
###--- Supplementary Fig. 3b ---###
plot_aligned_logos_combined(aln_lst = res[[1]], aln_df = res[[2]], weight_lst = res[[3]], ylab = 'Predicted tail-length change (nt)', fn = paste0('XL_oocyte_ISM_', kmer_len, 'mer_top_kmers_gain_logo'))

#--- logo plot from 8mer, aligned with custom score, for kmer loss
# filter kmer by performing one sample t-test against global mean
kmer_len <- 8
ism.sele <- kmer.loss %>% filter(nchar(kmer) == kmer_len)
avg_mean <- ism.sele %>% pull(avg) %>% mean(.)
avg_std <- ism.sele %>% pull(avg) %>% sd(.)
ism.sele <- ism.sele %>% mutate(z_score = (avg - avg_mean) / avg_std) %>%
	mutate(p_positive = one_sample_t_test(avg, std, n, diff = avg_mean, lower.tail = F),
				 p_negative = one_sample_t_test(avg, std, n, diff = avg_mean, lower.tail = T)) %>%
	mutate(log10p = ifelse(avg > avg_mean, -log10(pmax(p_positive, .Machine$double.xmin)), -log10(pmax(p_negative, .Machine$double.xmin))))

res <- plot_top_kmer_logo(as.data.frame(ism.sele %>% filter(log10p >= 2)), 
													weight = 'avg', mode = 'negative', n_seed_max = Inf, z_cutoff = (-3),
													align_score_cutoff = -0.1, score_match = 1, score_offset = -1, score_mismatch = -1.5, type = 'bits',
													fn = paste0('XL_oocyte_ISM_', kmer_len, 'mer_top_kmers_logo'))

# combined logo plot using the result from 'plot_top_kmer_logo'
###--- Fig. 2b ---###
plot_aligned_logos_combined(aln_lst = res[[1]], aln_df = res[[2]], weight_lst = res[[3]], ylab = 'Predicted tail-length change (nt)', fn = paste0('XL_oocyte_ISM_', kmer_len, 'mer_top_kmers_loss_logo'))


#------------------------------------------------------------------------------------------------------------------------------------
#--- examine predicted tail-length changes in an iterative exclusion manner in ISM results
kmer.gain <- read_delim('../Data/INN_ISM/XL_oocyte_ISM_utr_len_limit_300_gain_kmer_iterative_exclusion_pA_weight_avg_dis.txt')
kmer.gain <- kmer.gain %>% mutate(std = (dis/(n-1)) ** 0.5) %>% mutate(sem = std / (n ** 0.5))

kmer.loss <- read_delim('../Data/INN_ISM/XL_oocyte_ISM_utr_len_limit_300_loss_kmer_iterative_exclusion_pA_weight_avg_dis.txt')
kmer.loss <- kmer.loss %>% mutate(std = (dis/(n-1)) ** 0.5) %>% mutate(sem = std / (n ** 0.5))

#--- barplot of difference in mean pA length versus mean count number 
###--- Supplementary Fig. 3c ---###
kmer_tl_diff_by_ie_bplot(kmer.gain %>% arrange(round) %>% dplyr::slice(1:15), 
												 y_col = 'avg', mybreaks = seq(0, 8, 2), 
												 fn = 'XL_oocyte_ISM_6mer_gain_iterative_exclusion_pA_difference_barplot')

###--- Fig. 2c ---###
kmer_tl_diff_by_ie_bplot(kmer.loss %>% arrange(round) %>% dplyr::slice(1:15), 
												 y_col = 'avg', mybreaks = seq(-10, 0, 2),
												 fn = 'XL_oocyte_ISM_6mer_loss_iterative_exclusion_pA_difference_barplot')

#------------------------------------------------------------------------------------------------------------------------------------
#--- examine predicted tail-length changes for select kmer by each nucleotide
###--- Supplementary Fig. 3g ---###

kmers <- h5read('../Data/INN_ISM/XL_oocyte_ISM_utr_len_limit_300_loss_kmer_iterative_exclusion_pA_weight_avg_dis.h5', 'kmer')
kmer.loss <- h5read('../Data/INN_ISM/XL_oocyte_ISM_utr_len_limit_300_loss_kmer_iterative_exclusion_pA_weight_avg_dis.h5', 'stats')

sele.loss.avg <- t(kmer.loss[2,,,1]) # 3 (count, mean, dis) x 6 (position) x 4 (ACGT) x n (kmer count)
sele.loss.n <- t(kmer.loss[1,,,1]) 
sele.loss.dis <- t(kmer.loss[3,,,1]) 
sele.loss.std <- (sele.loss.dis / (sele.loss.n - 1)) ** 0.5
sele.loss.sem <- sele.loss.std / ((sele.loss.n - 1) ** 0.5)

df.plot <- inner_join(stat_matrix_to_kmer_df(sele.loss.avg) %>% dplyr::rename('avg' = 'val'),
											stat_matrix_to_kmer_df(sele.loss.sem) %>% dplyr::rename('sem' = 'val'),
											by = 'kmer') %>%
	arrange(desc(avg)) %>% mutate(kmer = factor(kmer, levels = .$kmer))

res <- kmer_single_nt_mutation_tl_diff_bplot(df.plot,
																						 col_x = 'kmer', 
																						 col_y = 'avg',
																						 col_y_sem = 'sem',
																						 mybreaks = seq(-12,0,2), pad = 0.1, 
																						 fn = 'XL_oocyte_ISM_utr_len_limit_300_loss_PAS_all_mutants_barplot')

#------------------------------------------------------------------------------------------------------------------------------------
# Examine tail length change by position from ISM results for select genes 
# make a combined plot including a heatmap, a line plot and a logo plot
###--- Fig. 2g ---###
###--- Fig. 2h ---###
###--- Supplementary Fig. 3i ---###
###--- Supplementary Fig. 3j ---###

ism_path <- '../Data/INN_ISM/XL_oocyte_ISM_average.h5'
# h5ls(ism_path)
pa_id_vec <- h5read(ism_path, 'id')
diff_score_mat <- h5read(ism_path, 'diff_scores')
seqs <- h5read(ism_path, 'seq')
pa.site <- read.table('../Data/Annotation/XL_pA_site_PMID_38460509_Table_S2.txt', sep = '\t', header = T)

sele.vec = c('mos.L', 'tpx2.L', 'ccnb1.2.L', 'ccnb2.L', 
						 'aurkaip1.L', 
						 'dbf4.L',
						 'mad2l1.L', 'lima1.L', 'magoh.S',
						 'atp1a1.S')

res.lst = list()
if (startsWith(sele.vec[1], 'Chr')){
	pa_id_sele <- sele.vec
	gene_name_sele <- pa.site[pa.site$pA_id %in% pa_id_vec, 'gene_name']
} else {
	gene_name_sele <- sele.vec
	pa_id_sele <- pa.site[pa.site$gene_name %in% gene_name_sele, 'pA_id']
	gene_name_sele <- pa.site[pa.site$pA_id %in% pa_id_sele, 'gene_name'] # there could be more than on pA sites for each gene
}

res.lst <- lapply(1:length(pa_id_sele), function(i){
	seq_sele <- seqs[pa_id_vec == pa_id_sele[i]]
	score_sele <- diff_score_mat[,,pa_id_vec == pa_id_sele[i]]
	score_sele[is.na(score_sele)] <- 0 # replace NaN with 0
	colnames(score_sele) <- unlist(strsplit(gsub('T', 'U', seq_sele) , ''))
	rownames(score_sele) <- c('A', 'C', 'G', 'U')
	fn_idx <- paste0(gene_name_sele[i], '_', pa_id_sele[i])
	res <- select_gene_ism_score_heatmap_w_logo(score_sele, utr_label = NULL, pos_range = c(-103, -1), fn = fn_idx, offset = 3, 
																							xbreak = 10, g_fsize = 5, block_size = 0.07, nt_label_scaler = 1, line_w = 0.75/.pt)
	return(res)
})
