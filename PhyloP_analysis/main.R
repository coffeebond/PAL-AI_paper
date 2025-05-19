library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(rhdf5)
library(parallel)
library(pbmcapply)
library(Biostrings)
source('helper.R')

# get the bed file
utr_bed <- GRanges(import('../Data/Annotation/HS_3UTR_annotation_PMID_38460509_Table_S3.bed', format = "BED"))
utr_bed$pa_id <- str_split_i(utr_bed$name, '__', 1)

# get the ISM results
ism_res_file <- '../Data/INN_ISM/HS_oocyte_ISM_average.h5'
diff_scores <- h5read(ism_res_file, name = 'diff_scores') # n_channel x l_seq_length x m_count_seqs
pa_ids <- h5read(ism_res_file, name = 'id')
seqs <- DNAStringSet(h5read(ism_res_file, name = 'seq'))
names(seqs) <- pa_ids
offset <- 3 # the last 3 nucleotides were appended 'AAA' 
l_max <- 2000

# average outcome for each nucleotide of each sequence
avg_diff_scores <- diff_scores %>% apply(., 3, function(mat){colMeans(mat, na.rm = TRUE)}) 
diff_df <- lapply(seq_along(pa_ids), function(x){
	list('pa_id' = pa_ids[x],
			 'dis2end' = seq((l_max - offset - 1), 0, -1),
			 'tl_diff' = avg_diff_scores[1:(l_max - offset), x]
	)
}) %>% bind_rows %>% filter(!is.na(tl_diff))


# filter the bed file (this may take a few minutes)
utr_bed_filtered <- utr_bed[utr_bed$pa_id %in% pa_ids]
modified_exons <- unlist(GRangesList(pbmclapply(split(utr_bed_filtered, utr_bed_filtered$pa_id), function(x){adjust_exons(x, l_max = (l_max-offset))}, mc.cores = 40)))

# get the per-nucleotide positions (this may take a few minutes)
per_nt_pos <- unlist(GRangesList(pbmclapply(split(modified_exons, modified_exons$pa_id), compute_relative_positions, mc.cores = 40)))

# get the PhyloP scores in a bigwig file
bg <- '/lab/solexa_bartel/coffee/Sequencing/PhyloP/hg38.phyloP100way.bw'
scores <- import(bg, format = 'BigWig', which = modified_exons)

# add PhyloP scores to the per-nucleotide position table
df.pnp <- left_join(as_tibble(per_nt_pos) %>% select(seqnames, start, end, strand, dis2end, pa_id, name), 
										as_tibble(scores) %>% select(seqnames, start, end, score), by = c('seqnames', 'start', 'end'))


# annotate CPE and PAS positions
cpe_df <- tile_motif_positions(seqs, 'TTTTA', fixed = 'subject') %>% 
	dplyr::rename('cpe' = 'flag') %>% mutate(dis2end = l_max - offset - pos)
pas_df <- tile_motif_positions(seqs, 'AWTAAA', fixed = 'subject') %>% 
	dplyr::rename('pas' = 'flag') %>% mutate(dis2end = l_max - offset - pos)


# combine tail-length change and PhyloP scores, and add CPE and PAS annotation/flag
df.combined <- left_join(diff_df, df.pnp %>% select(dis2end, pa_id, score, name),
										 by = c('pa_id', 'dis2end')) %>% filter(!is.na(score)) %>%
	left_join(., cpe_df %>% select(-pos), by = c('pa_id', 'dis2end')) %>%
		left_join(., pas_df %>% select(-pos), by = c('pa_id', 'dis2end')) 

# make a boxplot
###--- Fig. 7a ---###
dis2end_cutoff <- 100
diff_breaks <- c(-Inf, seq(-25,25,10), Inf)
group_labels <- paste0("(", head(diff_breaks, -1), ", ", tail(diff_breaks, -1), "]")
df.plot <- df.combined %>% filter(dis2end < dis2end_cutoff) %>%
	mutate(group = cut(tl_diff, 
										 breaks = diff_breaks,
										 labels = group_labels,
										 include.lowest = TRUE, right = FALSE)) %>%
		mutate(group = factor(group, levels = group_labels))

res <- phylop_box_plot(df.plot,
											 group_col = 'group', y_col = 'score',
											 group_label = levels(df.plot$group), 
											 color_sele = brewer.pal(11, 'RdYlBu')[rev(seq(3,9,1))], 
											 ybreaks = seq(-1,5,0.5), 
											 fn = 'PhyloP_score_by_tl_diff_boxplot')
print(res, n = Inf)

##-- heatmap for p values
###--- Supplementary Fig. 9a ---### 
pval_heatmap(res %>% mutate(group_x = factor(group_x, levels = group_labels),
														group_y = factor(group_y, levels = group_labels)), 
						 col_x = 'group_x', col_y = 'group_y', col_val = 'pval',
						 fn = 'PhyloP_score_by_tl_diff_wilcon_test_pvalue_comparion_heatmap')


# exclude those disrupt a PAS and make a boxplot
dis2end_cutoff <- 100
diff_breaks <- c(-Inf, seq(-25,25,10), Inf)
df.plot <- df.combined %>% filter(dis2end < dis2end_cutoff & is.na(pas)) %>%
	mutate(group = cut(tl_diff, 
										 breaks = diff_breaks,
										 labels = paste0("(", head(diff_breaks, -1), ", ", tail(diff_breaks, -1), "]"), include.lowest = TRUE, right = FALSE))
res <- phylop_box_plot(df.plot,
											 group_col = 'group', y_col = 'score',
											 group_label = levels(df.plot$group), color_sele = NULL, ybreaks = NULL, 
											 fn = 'PhyloP_score_by_tl_diff_no_PAS_disruption_boxplot')