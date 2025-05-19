source('./helper.R')
library('tidyverse')

###------------------------------------------------------------------------------------------------------------------------------------------------------
# scatter plot for tail change associated with position-less 3-mers when one UUUUA is present
###--- Supplementary Fig. 3d ---###
df_c <- read_delim('../Data/Library_N60_PASmos_data/N60_PASmos_library_frog_oocyte_uninjected_Tail_seq_single_CPE_TL_stat_kmer_position_less.txt.gz',
									 col_names = c('kmer', 'count', 'tl_avg', 'tl_sd', 'tl_log_avg', 'tl_log_sd')) %>% 
	mutate(tl_sd = (tl_sd / (count - 1)) ** 0.5, tl_log_sd = (tl_log_sd / (count - 1)) ** 0.5)
df_e <- read_delim('../Data/Library_N60_PASmos_data/N60_PASmos_library_frog_oocyte_progesterone_7hr_Tail_seq_single_CPE_TL_stat_kmer_position_less.txt.gz',
									 col_names = c('kmer', 'count', 'tl_avg', 'tl_sd', 'tl_log_avg', 'tl_log_sd')) %>%
	mutate(tl_sd = (tl_sd / (count - 1)) ** 0.5, tl_log_sd = (tl_log_sd / (count - 1)) ** 0.5)

res <- tl_change_vs_count_splot(ctrl = df_c, expt = df_e, 
																kmer_len = 3, count_cutoff = 1, ybreaks = seq(12,18,1), xlim = c(100000, 3000000), show_sig = T,
																kmer_sele_lst = list(c('TTT'), c('TGT'), c('GTT')), color_labels = NULL,
																color_sele = brewer.pal(9,'Set1')[1:3], plot_label = T,
																fn = 'XL_oocyte_library_B3_3mer_one_TTTTA_7h_vs_0h_pA_difference_vs_count_splot')

###------------------------------------------------------------------------------------------------------------------------------------------------------
# scatter plot for tail change associated with position 3-mers when one UUUUA is present
###--- Fig. 2d ---###
df_c <- read_delim('../Data/Library_N60_PASmos_data/N60_PASmos_library_frog_oocyte_uninjected_Tail_seq_single_CPE_kmer_rel_pos_pA_weight_avg_dis.txt.gz',
									 col_names = TRUE) %>% dplyr::rename("position" = "rel_pos", "count" = "n", "tl_avg" = "avg", "tl_sd" = "dis") %>%
	mutate(tl_sd = (tl_sd / (count - 1)) ** 0.5)
df_e <- read_delim('../Data/Library_N60_PASmos_data/N60_PASmos_library_frog_oocyte_progesterone_7hr_Tail_seq_TL_stat_single_CPE_kmer_rel_pos_pA_weight_avg_dis.txt.gz',
									 col_names = TRUE) %>% dplyr::rename("position" = "rel_pos", "count" = "n", "tl_avg" = "avg", "tl_sd" = "dis") %>%
	mutate(tl_sd = (tl_sd / (count - 1)) ** 0.5)
df.plot <- kmer_tl_diff_test(ctrl = as.data.frame(df_c), expt = as.data.frame(df_e), kmer_len = 3, count_cutoff = 1, t_mode = 'n', padj_inf = TRUE) %>% 
	mutate(dis2cpe = ifelse(pos < 0, pos + 2, pos - 4)) %>% as_tibble
res <- kmer_tl_diff_by_rel_pos_complex_plot(df.plot, col_x = 'dis2cpe', col_y = 'diff_avg', 
																						kmer_sele = c('TTT', 'TGT', 'GTT'), 
																						xlim = c(-10, 10), ybreaks = seq(10,20,2), 
																						fn = 'XL_oocyte_library_B3_3mer_one_TTTTA_7h_vs_0h_pA_difference_vs_relative_position_to_CPE_scatter_plot')



###------------------------------------------------------------------------------------------------------------------------------------------------------
# Examine tail-length change by position for nucleotide composition

# load data
ctrl <- read_delim('../Data/Library_N60_PASmos_data/N60_PASmos_library_frog_oocyte_uninjected_Tail_seq_TL_stat_kmer_by_position.txt.gz',
									 col_names = c('kmer','position', 'count', 'tl_avg', 'tl_dis', 'tl_log_avg', 'tl_log_dis'))
expt <- read_delim('../Data/Library_N60_PASmos_data/N60_PASmos_library_frog_oocyte_progesterone_7hr_Tail_seq_TL_stat_kmer_by_position.txt.gz',
									 col_names = c('kmer','position', 'count', 'tl_avg', 'tl_dis', 'tl_log_avg', 'tl_log_dis'))

# logo plot for relative tail-length change by position 
res <- inner_join(ctrl %>% filter(nchar(kmer) == 1),
									expt %>% filter(nchar(kmer) == 1),
									by = c('kmer', 'position'), suffix = c('_1', '_2')) %>%
	mutate(diff = tl_avg_2 - tl_avg_1) %>% group_by(position) %>% mutate(diff_nor = diff - mean(diff)) %>%
	select(kmer, position, diff_nor) %>% pivot_wider(names_from = 'position', values_from = 'diff_nor') %>% 
	mutate(kmer = gsub('T', 'U', kmer)) %>% 
	column_to_rownames(var = "kmer")
position_logo_plot(as.matrix(res), 
									 xbreaks = seq(10, 60, 10), xlabels = rev(seq(0, 50, 10)),
									 ybreaks = seq(-25,25,5)/10,
									 fn = 'CPEmos_N60_library_frog_oocyte_progesterone_7hr_over_0hr_TL_difference_by_position_logo_plot')