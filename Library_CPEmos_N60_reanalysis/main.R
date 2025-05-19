source('./helper.R')
library('tidyverse')


###------------------------------------------------------------------------------------------------------------------------------------------------------
# Logo plot for tail-length change by position for all positions in the random region
# load data
ctrl <- read_delim('../Data/Library_CPEmos_N60_data/CPEmos_N60_library_frog_oocyte_uninjected_Tail_seq_TL_stat_kmer_by_position.txt.gz',
									 col_names = c('kmer','position', 'count', 'tl_avg', 'tl_dis', 'tl_log_avg', 'tl_log_dis'))
expt <- read_delim('../Data/Library_CPEmos_N60_data/CPEmos_N60_library_frog_oocyte_progesterone_5hr_Tail_seq_TL_stat_kmer_by_position.txt.gz',
									 col_names = c('kmer','position', 'count', 'tl_avg', 'tl_dis', 'tl_log_avg', 'tl_log_dis'))

# logo plot for relative tail-lenght change by position 
res <- inner_join(ctrl %>% filter(nchar(kmer) == 1),
									expt %>% filter(nchar(kmer) == 1),
									by = c('kmer', 'position'), suffix = c('_1', '_2')) %>%
	mutate(diff = tl_avg_2 - tl_avg_1) %>% group_by(position) %>% mutate(diff_nor = diff - mean(diff)) %>%
	select(kmer, position, diff_nor) %>% pivot_wider(names_from = 'position', values_from = 'diff_nor') %>% 
	mutate(kmer = gsub('T', 'U', kmer)) %>% 
	column_to_rownames(var = "kmer")
position_logo_plot(as.matrix(res), 
									 xbreaks = seq(1, 51, 10), xlabels = seq(0, 50, 10),
									 ybreaks = seq(-25,25,5)/10,
									 fn = 'CPEmos_N60_library_frog_oocyte_progesterone_5hr_over_0hr_TL_difference_by_position_logo_plot')

###------------------------------------------------------------------------------------------------------------------------------------------------------
# scatter plot of difference in mean pA length versus mean count number (exclude AWUAAA)
###--- Supplementary Fig. 3h ---###
# load data
ctrl <- read_delim('../Data/Library_CPEmos_N60_data/CPEmos_N60_library_frog_oocyte_uninjected_Tail_seq_no_AWUAAA_kmer_pos_less_pA_weight_avg_dis.txt.gz', col_names = F) 
colnames(ctrl) = c('kmer', 'count', 'tl_avg', 'tl_sd', 'tl_log_avg', 'tl_log_sd')
expt <- read_delim('../Data/Library_CPEmos_N60_data/CPEmos_N60_library_frog_oocyte_progesterone_5hr_Tail_seq_no_AWUAAA_kmer_pos_less_pA_weight_avg_dis.txt.gz', col_names = F) 
colnames(expt) = c('kmer', 'count', 'tl_avg', 'tl_sd', 'tl_log_avg', 'tl_log_sd')
res <- tl_change_vs_count_splot(ctrl = ctrl, expt = expt, 
																kmer_len = 6, count_cutoff = 1, ybreaks = seq(-6,0,1), xlim = c(10000, 5000000), show_sig = T,
																kmer_sele_lst = list(c('AGTAAA', 'ATAAAG', 'TATAAA')), 
																color_labels = NULL, plot_label = T,
																color_sele = brewer.pal(8,'Set2')[1],
																fn = 'XL_oocyte_library_B2_no_AWUAAA_6mer_5h_vs_0h_pA_difference_vs_count_splot')
