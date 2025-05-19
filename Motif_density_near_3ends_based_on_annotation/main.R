library(tidyverse)
library(Biostrings)
source('helper.R')

#------------------------------------------------------------------------------------------------------------------------------------
# density plot for CPE and PAS in Xenopus
###--- Supplementary Fig. 3e ---###

l_max <- 1000 # max length of 3' UTR to analyze (from 3' ends)

# load 3'-UTR sequence and filter out "uncertain" entries
utr3 <- readDNAStringSet('../Data/Sequences/XENLA_10.1_GCF_XBmodels_annotation_fixed_UTR3_polished_by_PA_site_stitched.fa')
utr3_filtered <- utr3[!grepl('uncertain', names(utr3))]

# get denominator for calculating density (number of sequences with nucleotides at each position)
df_pos_count <- tibble(
	'pos' = c(-l_max:-1),
	'n_total' = sapply(1:l_max, function(x){
		sum(width(utr3_filtered) > (l_max - x))
	})
)

# find CPE motif occurrences
cpe_pos_lst <- start(vmatchPattern('TTTTA', utr3_filtered))
df.cpe <- tibble(
	'pos' = unlist(Map('-', cpe_pos_lst, width(utr3_filtered) + 1)) 
) %>% filter(pos >= (-l_max)) %>% group_by(pos) %>% summarize(n = n()) %>% 
	right_join(., df_pos_count, by = 'pos') %>% replace_na(list(n = 0)) %>%
		mutate(freq = n / n_total, motif = 'TTTTA') 

# find PAS motif occurrences
df.pas <- tibble(
	'pos' = unlist(Map('-', start(vmatchPattern('AWTAAA', utr3_filtered, fixed = FALSE)), width(utr3_filtered) + 1)) 
) %>% filter(pos >= (-l_max)) %>% group_by(pos) %>% summarize(n = n()) %>% 
	right_join(., df_pos_count, by = 'pos') %>% replace_na(list(n = 0)) %>%
	mutate(freq = n / n_total, motif = 'AWTAAA') 

# combine dataframes and remove rows with position beyond the 3'-UTR regions
df.plot <- bind_rows(df.cpe, df.pas) %>% 
	mutate(motif = gsub('T', 'U', motif),
				 dis2end = pos + nchar(motif)) %>%
		mutate(adj_x = dis2end + l_max,
					 motif = factor(motif, levels = c('UUUUA', 'AWUAAA'))) %>%
			filter(adj_x <= l_max)

res <- motif_density_line_plot(df.plot, 
															 col_x = 'adj_x', col_y = 'freq', 
															 col_group = 'motif' , 
															 xbreaks = seq(500,1000,100), xlabels = seq(500, 0, -100), x_lim = c(500, 1000),
															 ybreaks = seq(0, 0.1, 0.02),  
															 x_title = 'Distance to 3\' ends (nt)',
															 y_title = 'Density', 
															 group_color = NULL, 
															 fn = 'XL_oocyte_3UTR_select_motif_density')

res <- motif_density_line_plot(df.plot, 
															 col_x = 'adj_x', col_y = 'freq', 
															 col_group = 'motif' , 
															 xbreaks = seq(950,1000,10), xlabels = seq(50, 0, -10), x_lim = c(950, 1000),
															 ybreaks = seq(0, 0.1, 0.02),  
															 x_title = 'Distance to 3\' ends (nt)',
															 y_title = 'Density', 
															 group_color = NULL, 
															 fn = 'XL_oocyte_3UTR_select_motif_density_inlet')

#------------------------------------------------------------------------------------------------------------------------------------
# density plot for CPE and PAS in mice
l_max <- 1000 # max length of 3' UTR to analyze (from 3' ends)

# load 3'-UTR sequence and filter out "uncertain" entries
utr3 <- readDNAStringSet('/lab/solexa_bartel/coffee/Sequencing/Reference/Transcriptome/MM/Oocyte/UTR3/gencode.vM10.primary_assembly.annotation_annotation_fixed_UTR3_polished_by_PA_site_stitched.fa')
utr3_filtered <- utr3[!grepl('uncertain', names(utr3))]

# get denominator for calculating density (number of sequences with nucleotides at each position)
df_pos_count <- tibble(
	'pos' = c(-l_max:-1),
	'n_total' = sapply(1:l_max, function(x){
		sum(width(utr3_filtered) > (l_max - x))
	})
)

# find CPE motif occurrences
cpe_pos_lst <- start(vmatchPattern('TTTTA', utr3_filtered))
df.cpe <- tibble(
	'pos' = unlist(Map('-', cpe_pos_lst, width(utr3_filtered) + 1)) 
) %>% filter(pos >= (-l_max)) %>% group_by(pos) %>% summarize(n = n()) %>% 
	right_join(., df_pos_count, by = 'pos') %>% replace_na(list(n = 0)) %>%
	mutate(freq = n / n_total, motif = 'TTTTA') 

# find PAS motif occurrences
df.pas <- tibble(
	'pos' = unlist(Map('-', start(vmatchPattern('AWTAAA', utr3_filtered, fixed = FALSE)), width(utr3_filtered) + 1)) 
) %>% filter(pos >= (-l_max)) %>% group_by(pos) %>% summarize(n = n()) %>% 
	right_join(., df_pos_count, by = 'pos') %>% replace_na(list(n = 0)) %>%
	mutate(freq = n / n_total, motif = 'AWTAAA') 

# combine dataframes and remove rows with position beyond the 3'-UTR regions
df.plot <- bind_rows(df.cpe, df.pas) %>% 
	mutate(motif = gsub('T', 'U', motif),
				 dis2end = pos + nchar(motif)) %>%
	mutate(adj_x = dis2end + l_max,
				 motif = factor(motif, levels = c('UUUUA', 'AWUAAA'))) %>%
	filter(adj_x <= l_max)

res <- motif_density_line_plot(df.plot, 
															 col_x = 'adj_x', col_y = 'freq', 
															 col_group = 'motif' , 
															 xbreaks = seq(500,1000,100), xlabels = seq(500, 0, -100), x_lim = c(500, 1000),
															 ybreaks = seq(0, 0.1, 0.02),  
															 x_title = 'Distance to 3\' ends (nt)',
															 y_title = 'Density', 
															 group_color = NULL, 
															 fn = 'MM_oocyte_3UTR_select_motif_density')

res <- motif_density_line_plot(df.plot, 
															 col_x = 'adj_x', col_y = 'freq', 
															 col_group = 'motif' , 
															 xbreaks = seq(950,1000,10), xlabels = seq(50, 0, -10), x_lim = c(950, 1000),
															 ybreaks = seq(0, 0.1, 0.02),  
															 x_title = 'Distance to 3\' ends (nt)',
															 y_title = 'Density', 
															 group_color = NULL, 
															 fn = 'MM_oocyte_3UTR_select_motif_density_inlet')


#------------------------------------------------------------------------------------------------------------------------------------
# density plot for CPE and PAS in humans
l_max <- 1000 # max length of 3' UTR to analyze (from 3' ends)

# load 3'-UTR sequence and filter out "uncertain" entries
utr3 <- readDNAStringSet('/lab/solexa_bartel/coffee/Sequencing/Reference/Transcriptome/HS/Oocyte/UTR3/gencode.v25.annotation_annotation_fixed_UTR3_polished_by_PA_site_stitched.fa')
utr3_filtered <- utr3[!grepl('uncertain', names(utr3))]

# get denominator for calculating density (number of sequences with nucleotides at each position)
df_pos_count <- tibble(
	'pos' = c(-l_max:-1),
	'n_total' = sapply(1:l_max, function(x){
		sum(width(utr3_filtered) > (l_max - x))
	})
)

# find CPE motif occurrences
cpe_pos_lst <- start(vmatchPattern('TTTTA', utr3_filtered))
df.cpe <- tibble(
	'pos' = unlist(Map('-', cpe_pos_lst, width(utr3_filtered) + 1)) 
) %>% filter(pos >= (-l_max)) %>% group_by(pos) %>% summarize(n = n()) %>% 
	right_join(., df_pos_count, by = 'pos') %>% replace_na(list(n = 0)) %>%
	mutate(freq = n / n_total, motif = 'TTTTA') 

# find PAS motif occurrences
df.pas <- tibble(
	'pos' = unlist(Map('-', start(vmatchPattern('AWTAAA', utr3_filtered, fixed = FALSE)), width(utr3_filtered) + 1)) 
) %>% filter(pos >= (-l_max)) %>% group_by(pos) %>% summarize(n = n()) %>% 
	right_join(., df_pos_count, by = 'pos') %>% replace_na(list(n = 0)) %>%
	mutate(freq = n / n_total, motif = 'AWTAAA') 

# combine dataframes and remove rows with position beyond the 3'-UTR regions
df.plot <- bind_rows(df.cpe, df.pas) %>% 
	mutate(motif = gsub('T', 'U', motif),
				 dis2end = pos + nchar(motif)) %>%
	mutate(adj_x = dis2end + l_max,
				 motif = factor(motif, levels = c('UUUUA', 'AWUAAA'))) %>%
	filter(adj_x <= l_max)

res <- motif_density_line_plot(df.plot, 
															 col_x = 'adj_x', col_y = 'freq', 
															 col_group = 'motif' , 
															 xbreaks = seq(500,1000,100), xlabels = seq(500, 0, -100), x_lim = c(500, 1000),
															 ybreaks = seq(0, 0.1, 0.02),  
															 x_title = 'Distance to 3\' ends (nt)',
															 y_title = 'Density', 
															 group_color = NULL, 
															 fn = 'HS_oocyte_3UTR_select_motif_density')

res <- motif_density_line_plot(df.plot, 
															 col_x = 'adj_x', col_y = 'freq', 
															 col_group = 'motif' , 
															 xbreaks = seq(950,1000,10), xlabels = seq(50, 0, -10), x_lim = c(950, 1000),
															 ybreaks = seq(0, 0.1, 0.02),  
															 x_title = 'Distance to 3\' ends (nt)',
															 y_title = 'Density', 
															 group_color = NULL, 
															 fn = 'HS_oocyte_3UTR_select_motif_density_inlet')
