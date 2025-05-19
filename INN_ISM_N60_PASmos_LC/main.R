library(tidyverse)
library(Biostrings)
library(rhdf5)
source('helper.R')

#------------------------------------------------------------------------------------------------------------------------------------
#--- examine predicted tail-length changes in an iterative exclusion manner in ISM results
kmer.gain <- read_delim('../Data/INN_ISM_snm_library_F044/XL_oocyte_N60_PASmos_LC_ISM_gain_kmer_iterative_exclusion_pA_weight_avg_dis.txt')
kmer.gain <- kmer.gain %>% mutate(std = (dis/(n-1)) ** 0.5) %>% mutate(sem = std / (n ** 0.5))

kmer.loss <- read_delim('../Data/INN_ISM_snm_library_F044/XL_oocyte_N60_PASmos_LC_ISM_loss_kmer_iterative_exclusion_pA_weight_avg_dis.txt')
kmer.loss <- kmer.loss %>% mutate(std = (dis/(n-1)) ** 0.5) %>% mutate(sem = std / (n ** 0.5))

#--- barplot of difference in mean pA length versus mean count number 
###--- Supplementary Fig. 4d ---###
kmer_tl_diff_by_ie_bplot(kmer.gain %>% arrange(round) %>% dplyr::slice(1:15), 
												 y_col = 'avg', mybreaks = seq(0, 30, 5), 
												 fn = 'XL_oocyte_N60_PASmos_LC_ISM_6mer_gain_iterative_exclusion_pA_difference_barplot')

###--- Supplementary Fig. 4e ---###
kmer_tl_diff_by_ie_bplot(kmer.loss %>% arrange(round) %>% dplyr::slice(1:15), 
												 y_col = 'avg', mybreaks = seq(-30, 0, 5),
												 fn = 'XL_oocyte_N60_PASmos_LC_ISM_6mer_loss_iterative_exclusion_pA_difference_barplot')

#------------------------------------------------------------------------------------------------------------------------------------
#--- examine motif insertion outcome by position
###--- Supplementary Fig. 4f ---###
df.mi <- read_delim('../Data/INN_ISM_N60_PASmos_LC/Motif_insertion_results_average.txt', delim = '\t')

df.plot <- df.mi %>% filter(motif %in% c('TTTTA')) %>%
	mutate(motif = gsub('T', 'U', motif)) %>%
	mutate(motif = factor(motif, levels = c('UUUUA'))) %>%
	filter(dis2end >= 0) # remove negative value because AAA was appended for the INN input

res <- motif_insertion_line_plot(df_in = df.plot %>% mutate(dis2end = dis2end - 21) %>% filter(dis2end >= 0), 
																 xbreak = 5, ybreaks = seq(0,50,10), pos_lims = c(0, 1), 
																 fn = 'XL_oocyte_INN_motif_insertion_position_effect')
