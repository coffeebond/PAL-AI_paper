library(tidyverse)
library(Biostrings)
source('helper.R')
#------------------------------------------------------------------------------------------------------------------------------------
# PAS and CPE
###--- Fig. 2e ---###
df.mi <- read_delim('../Data/INN_motif_insertion/Motif_insertion_results_average.txt', delim = '\t')

df.plot <- df.mi %>% filter(motif %in% c('TTTTA', 'AATAAA', 'ATTAAA')) %>%
	mutate(motif = gsub('T', 'U', motif)) %>%
		mutate(motif = factor(motif, levels = c('UUUUA', 'AAUAAA', 'AUUAAA'))) %>%
			filter(dis2end >= 0) # remove negative value because AAA was appended for the INN input

res <- motif_insertion_line_plot(df_in = df.plot, 
																 xbreak = 100, pos_lims = c(0.75, 1), 
																 fn = 'XL_oocyte_INN_motif_insertion_position_effect')

res <- motif_insertion_line_plot(df_in = df.plot, 
																 xbreak = 10, pos_lims = c(0.95, 1), 
																 fn = 'XL_oocyte_INN_motif_insertion_position_effect_inlet')

#------------------------------------------------------------------------------------------------------------------------------------
# All PAS motifs other than AAUAAA
###--- Supplementary Fig. 3f ---###
df.pas <- read_delim('../Data/INN_motif_insertion/Motif_insertion_top_PAS_masked_results_average.txt', delim = '\t')

df.plot <- df.pas %>% 
	mutate(motif = gsub('T', 'U', motif)) %>%
		mutate(motif = factor(motif, levels = unique(.$motif))) %>%
			filter(dis2end >= 0) # remove negative value because AAA was appended for the INN input

res <- motif_insertion_line_plot(df_in = df.plot, 
																 xbreak = 100, pos_lims = c(0.75, 1), 
																 col_sele = brewer.pal(12, 'Set3')[c(1:length(levels(df.plot$motif)))],
																 fn = 'XL_oocyte_INN_other_PAS_motif_insertion_position_effect')

res <- motif_insertion_line_plot(df_in = df.plot, 
																 xbreak = 10, pos_lims = c(0.95, 1), 
																 col_sele = brewer.pal(12, 'Set3')[c(1:length(levels(df.plot$motif)))],
																 fn = 'XL_oocyte_INN_other_PAS_motif_insertion_position_effect_inlet')

#------------------------------------------------------------------------------------------------------------------------------------
# relative positional effects between CPE and PAS
###--- Fig. 2f ---###
df.rel <- read_delim('../Data/INN_motif_insertion/Motif_insertion_average_CPE2PAS.txt.gz')

df.plot <- df.rel %>% group_by(dis2pas) %>% summarize(avg = mean(tl_diff), std = sd(tl_diff), n = n(), .groups = 'drop') %>%
	mutate(sem = std / (n ** 0.5), pos = -dis2pas) # positive dis2pas means PAS is downstream of CPE
pas_cpe_line_plot(df.plot, xbreaks = seq(-30,30,5), fn = 'XL_oocyte_INN_motif_insertion_CPE_PAS_relative_position_effect')
