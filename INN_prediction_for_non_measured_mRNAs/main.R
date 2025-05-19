library('biomaRt')
library('gprofiler2')
library('tidyverse')
source('helper.R')

###------------------------------------------------------------------------------------------------------------------------------------------------------
# examine conserved genes with tail lengthening during oocyte maturation in all three species
###------------------------------------------------------------------------------------------------------------------------------------------------------


###-----------------------------------------------------------------------------------------------
##-- This part generates tail-length change (either measured or predicted) by gene (merging mRNA isoforms for each gene)

tag_cutoff <- 50 # minimal number of poly(A) tags (reads) for calculating median values
frac_cutoff <- 0.9 # fraction of the primary isoform must be more than this to calculate tail-length change from measured data

##------------------------------------------------
#- Human data
hs_gv <- read_delim('../Data/Human_oocytes/GV_HRR596534_all_tags.txt.gz')
hs_m2 <- read_delim('../Data/Human_oocytes/MII_HRR596536_all_tags.txt.gz')
hs_pa <- read_delim('../Data/Annotation/HS_pA_site_PMID_38460509_Table_S2.txt')

# group by gene and calculate fraction of isoforms
hs.lst <- lapply(list(hs_gv, hs_m2), function(x){
	x %>% group_by(gene_id) %>% summarize(med = median(tl_by_palseq), n = n(), .groups = 'drop') %>%
		inner_join(., hs_pa %>% dplyr::select(pA_id, gene_name, confidence) %>% filter(confidence == 'high') %>% dplyr::select(-confidence), by = c('gene_id'='pA_id')) %>%
			group_by(gene_name) %>% mutate(frac = n/sum(n)) %>% ungroup(.) 
}) 
names(hs.lst) <- c('gv', 'm2')

# obtain the primary isoform
hs.pri <- 1:length(hs.lst) %>% lapply(function(x){
	hs.lst[[x]] %>% group_by(gene_name) %>% mutate(frac = max(frac)) %>%
		ungroup(.) %>% mutate(stage = names(hs.lst)[x])
}) %>% bind_rows(.)

# calculate tail-length change
df.hs <- lapply(hs.lst, function(x){
	x %>%	filter(frac >= frac_cutoff & n >= tag_cutoff)
}) %>% Reduce(function(x,y){inner_join(x, y, by = c('gene_id', 'gene_name'), suffix = c('_1', '_2'))},.) %>% 
	mutate(tl_diff = med_2 - med_1, tl_init = med_1, tl_diff_log2fc = log2(med_2/med_1)) %>% 
		dplyr::select(gene_name, gene_id, tl_diff, tl_init, tl_diff_log2fc) %>% dplyr::rename('pA_id' = 'gene_id') %>%
			dplyr::rename_with(~paste0('hs_', .x))

##------------------------------------------------
#- Mouse data
mm_gv <- bind_rows(read_delim('../Data/Mouse_oocytes/GSM7716872_Mouse_oocyte_mRNA_GV_rep1_PAL_seq_v4_processed.txt.gz'),
									 read_delim('../Data/Mouse_oocytes/GSM7716873_Mouse_oocyte_mRNA_GV_rep2_PAL_seq_v4_processed.txt.gz'))
mm_m2 <- bind_rows(read_delim('../Data/Mouse_oocytes/GSM7716874_Mouse_oocyte_mRNA_MII_rep1_PAL_seq_v4_processed.txt.gz'),
									 read_delim('../Data/Mouse_oocytes/GSM7716875_Mouse_oocyte_mRNA_MII_rep2_PAL_seq_v4_processed.txt.gz'))
mm_pa <- read_delim('../Data/Annotation/MM_pA_site_PMID_38460509_Table_S2.txt')

# group by gene and calculate fraction of isoforms
mm.lst <- lapply(list(mm_gv, mm_m2), function(x){
	x %>% dplyr::rename('gene_id' = 'Gene_id_or_pA_id') %>% group_by(gene_id) %>%
		summarize(med = median(tail_length), n =n(), .groups = 'drop') %>%
			inner_join(., mm_pa %>% dplyr::select(pA_id, gene_name, confidence) %>% filter(confidence == 'high') %>% dplyr::select(-confidence), by = c('gene_id'='pA_id')) %>%
				group_by(gene_name) %>% mutate(frac = n/sum(n)) %>% ungroup(.) 
}) 
names(mm.lst) <- c('gv', 'm2')

# obtain the primary isoform
mm.pri <- 1:length(mm.lst) %>% lapply(function(x){
	mm.lst[[x]] %>% group_by(gene_name) %>% mutate(frac = max(frac)) %>%
		ungroup(.) %>% mutate(stage = names(mm.lst)[x])
}) %>% bind_rows(.)

# calculate tail-length change
df.mm <- lapply(mm.lst, function(x){
	x %>%	filter(frac >= frac_cutoff & n >= tag_cutoff)
}) %>% Reduce(function(x,y){inner_join(x, y, by = c('gene_id', 'gene_name'), suffix = c('_1', '_2'))},.) %>% 
	mutate(tl_diff = med_2 - med_1, tl_init = med_1, tl_diff_log2fc = log2(med_2/med_1)) %>% 
		dplyr::select(gene_name, gene_id, tl_diff, tl_init, tl_diff_log2fc) %>% dplyr::rename('pA_id' = 'gene_id') %>%
			dplyr::rename_with(~paste0('mm_', .x))

##------------------------------------------------
#- Frog data
xl_pp0h <- read_delim('../Data/Frog_oocytes/GSM7716851_Frog_oocyte_mRNA_no_progesterone_PAL_seq_v3_processed.txt.gz')
xl_pp7h <- read_delim('../Data/Frog_oocytes/GSM7716855_Frog_oocyte_mRNA_progesterone_7hr_PAL_seq_v3_processed.txt.gz')
xl_pa <- read_delim('../Data/Annotation/XL_pA_site_PMID_38460509_Table_S2.txt')

# group by gene and calculate fraction of isoforms
xl.lst <- lapply(list(xl_pp0h, xl_pp7h), function(x){
	x %>% dplyr::rename('gene_id' = 'Gene_id_or_pA_id') %>% 
		group_by(gene_id) %>% summarize(med = median(tail_length), n = n(), .groups = 'drop') %>%
			inner_join(., xl_pa %>% dplyr::select(pA_id, gene_name, confidence) %>% filter(confidence == 'high') %>% dplyr::select(-confidence), by = c('gene_id'='pA_id')) %>%
				group_by(gene_name) %>% mutate(frac = n/sum(n)) %>% ungroup(.) 
})
names(xl.lst) <- c('pp0h', 'pp7h')

# obtain the primary isoform
xl.pri <- 1:length(xl.lst) %>% lapply(function(x){
	xl.lst[[x]] %>% group_by(gene_name) %>% mutate(frac = max(frac)) %>%
		ungroup(.) %>% mutate(stage = names(xl.lst)[x])
}) %>% bind_rows(.)

# calculate tail-length change
df.xl	<- lapply(xl.lst, function(x){
	x %>%	filter(frac >= frac_cutoff & n >= tag_cutoff)
}) %>% Reduce(function(x,y){inner_join(x, y, by = c('gene_id', 'gene_name'), suffix = c('_1', '_2'))},.) %>% 
	mutate(tl_diff = med_2 - med_1, tl_init = med_1, tl_diff_log2fc = log2(med_2/med_1)) %>% 
		dplyr::select(gene_name, gene_id, tl_diff, tl_init, tl_diff_log2fc) %>% dplyr::rename('pA_id' = 'gene_id') %>%
			dplyr::rename_with(~paste0('xl_', .x))

##------------------------------------------------
#- combine the measured tail-length-change values 
# read in the conversion table
df.gn.all <- read_delim('../Data/Annotation/HS_MM_XL_homologous_gene_name_table.txt', delim = '\t')

# examine data of homologous genes in all species
df.tl.all <- df.xl %>% dplyr::select(xl_gene_name, xl_tl_diff) %>% mutate(xl_gene_name_long = xl_gene_name) %>%
	separate(xl_gene_name, sep = "\\.", into = c('xl_gene_name', NA), extra = 'drop', fill = 'right') %>% 
		full_join(., df.gn.all %>% na.omit, by = 'xl_gene_name', relationship = "many-to-many") %>% 
			left_join(., df.hs %>% dplyr::select(hs_gene_name, hs_tl_diff), by = 'hs_gene_name') %>%
				left_join(., df.mm %>% dplyr::select(mm_gene_name, mm_tl_diff), by = 'mm_gene_name') %>% 
					filter(!if_any(ends_with('_name'), is.na)) %>% # remove rows that have NA in gene names
						mutate(tl_diff_min = reduce(across(ends_with('_tl_diff')), pmin, na.rm = T),
						 non_na = rowSums(!is.na(across(ends_with('_tl_diff'))))) # get the minimal TL change of all species

# some human genes may have homologous XL genes and/or MM genes, use the average
redundant_genes <- df.tl.all %>% group_by(hs_gene_name) %>% summarize(n = n()) %>% filter(n > 1) %>% pull(hs_gene_name)
df.tl.all.uniq <- lapply(redundant_genes, function(x){
	temp <- df.tl.all %>% filter(hs_gene_name == x)
	temp_out <- temp %>% dplyr::slice(1)
	temp_out$xl_tl_diff <- ifelse(all(is.na(temp$xl_tl_diff)), NA, mean(temp$xl_tl_diff, na.rm = T))
	temp_out$mm_tl_diff <- ifelse(all(is.na(temp$mm_tl_diff)), NA, mean(temp$mm_tl_diff, na.rm = T))
	tl_vec <- c(temp_out$hs_tl_diff, temp_out$mm_tl_diff, temp_out$xl_tl_diff)
	temp_out$tl_diff_min <- ifelse(all(is.na(tl_vec)), NA, min(tl_vec, na.rm = T))
	temp_out$non_na <- sum(!is.na(tl_vec)) # number of non-missing measured values for each gene
	return(temp_out)
}) %>% bind_rows(.) %>%
	bind_rows(., df.tl.all %>% filter(!(hs_gene_name %in% redundant_genes)))

##------------------------------------------------
#- get the predicted tail-length-change values

# the predicted value for the primary isoform will be used for each gene
# the fraction of the primary isoform must be more than this to be included
pri_frac_cutoff <- 0.5  

hs_tl <- read_delim('../Data/INN_model_predictions/XL_HS_MM_MINN_model_predicts_all_HS_mRNAs.txt')
mm_tl <- read_delim('../Data/INN_model_predictions/XL_HS_MM_MINN_model_predicts_all_MM_mRNAs.txt')
xl_tl <- read_delim('../Data/INN_model_predictions/XL_INN_model_predicts_all_XL_mRNAs.txt')

df_pred_lst <- lapply(list(list(hs_pa, hs_tl, hs.pri),
													 list(mm_pa, mm_tl, mm.pri),
													 list(xl_pa, xl_tl, xl.pri)), 
											function(lst){
												lst[[1]] %>% dplyr::select(pA_id, gene_name, confidence) %>%
													filter(confidence == 'high') %>% dplyr::select(-confidence) %>%
														inner_join(lst[[2]], ., by = c('id' = 'pA_id')) %>% 
															right_join(., lst[[3]] %>% filter(frac >= pri_frac_cutoff & (stage == 'gv' | stage == 'pp0h')), 
																		 by = c('id' = 'gene_id', 'gene_name' = 'gene_name'))
											})

# examine data of homologous genes in all species
df.tl.pred <- df_pred_lst[[3]] %>% dplyr::select(gene_name, y_pred) %>% dplyr::rename(xl_gene_name = gene_name, xl_tl_pred = y_pred) %>%
	mutate(xl_gene_name_long = xl_gene_name) %>%
		separate(xl_gene_name, sep = "\\.", into = c('xl_gene_name', NA), extra = 'drop', fill = 'right') %>% 
			full_join(., df.gn.all %>% na.omit, by = 'xl_gene_name', relationship = "many-to-many") %>% 
				left_join(., df_pred_lst[[1]] %>% dplyr::select(gene_name, y_pred) %>% dplyr::rename(hs_gene_name = gene_name, hs_tl_pred = y_pred), by = 'hs_gene_name', relationship = "many-to-many") %>%
					left_join(., df_pred_lst[[2]] %>% dplyr::select(gene_name, y_pred) %>% dplyr::rename(mm_gene_name = gene_name, mm_tl_pred = y_pred), by = 'mm_gene_name', relationship = "many-to-many") %>% 
						filter(!if_any(ends_with('_name'), is.na)) %>% # remove rows that have NA in gene names
							mutate(tl_pred_min = reduce(across(ends_with('_tl_pred')), pmin, na.rm = T)) # the minimal TL change of all species

# some human genes may have homologous XL genes and/or MM genes, use the average
redundant_genes <- df.tl.pred %>% group_by(hs_gene_name) %>% summarize(n = n()) %>% filter(n > 1) %>% pull(hs_gene_name)
df.tl.pred.uniq <- lapply(redundant_genes, function(x){
	temp <- df.tl.pred %>% filter(hs_gene_name == x)
	temp_out <- temp %>% dplyr::slice(1)
	temp_out$xl_tl_pred <- ifelse(all(is.na(temp$xl_tl_pred)), NA, mean(temp$xl_tl_pred, na.rm = T))
	temp_out$mm_tl_pred <- ifelse(all(is.na(temp$mm_tl_pred)), NA, mean(temp$mm_tl_pred, na.rm = T))
	tl_vec <- c(temp_out$hs_tl_pred, temp_out$mm_tl_pred, temp_out$xl_tl_pred)
	temp_out$tl_pred_min <- ifelse(all(is.na(tl_vec)), NA, min(tl_vec, na.rm = T))
	return(temp_out)
}) %>% bind_rows(.) %>%
	bind_rows(., df.tl.pred %>% filter(!(hs_gene_name %in% redundant_genes)))

##------------------------------------------------
#- combine measured data with predicted data and export the data
df.combined <- left_join(df.tl.all.uniq %>% dplyr::select(-xl_gene_name_long), 
												 df.tl.pred.uniq %>% dplyr::select(hs_gene_name, hs_tl_pred, mm_tl_pred, xl_tl_pred), by = 'hs_gene_name') %>%
	mutate(hs_tl = ifelse(is.na(hs_tl_diff) & (!is.na(hs_tl_pred)), hs_tl_pred, hs_tl_diff),
				 hs_flag = ifelse(is.na(hs_tl_diff) & (!is.na(hs_tl_pred)), TRUE, FALSE),
				 mm_tl = ifelse(is.na(mm_tl_diff) & (!is.na(mm_tl_pred)), mm_tl_pred, mm_tl_diff),
				 mm_flag = ifelse(is.na(mm_tl_diff) & (!is.na(mm_tl_pred)), TRUE, FALSE),
				 xl_tl = ifelse(is.na(xl_tl_diff) & (!is.na(xl_tl_pred)), xl_tl_pred, xl_tl_diff),
				 xl_flag = ifelse(is.na(xl_tl_diff) & (!is.na(xl_tl_pred)), TRUE, FALSE)) %>%
		dplyr::select(hs_gene_name, hs_tl, mm_tl, xl_tl, hs_flag, mm_flag, xl_flag) %>%
			mutate(tl_min = reduce(across(ends_with('_tl')), pmin, na.rm = T), # minimal of the tail-length change in all 3 species (regardless of measured or predicted)
						 non_na = rowSums(!is.na(across(ends_with('_tl')))), # number of species with either measured or predicted tail-length change
						 n_pred = rowSums(across(ends_with('_flag')))) # number of species with predicted tail-length change

# export the table to a file
out_folder <- './'
write_delim(df.combined, file = paste0(out_folder, 'HS_MM_XL_mRNA_TL_change_measurements_with_predictions_all_genes.txt'), delim = '\t')

# If the previous steps have completed and the table has been written, it can be reloaded back with the following line
#df.combined <- read_delim('HS_MM_XL_mRNA_TL_change_measurements_with_predictions_all_genes.txt', delim = '\t')
#df.gn.all <- read_delim('../Data/Annotation/HS_MM_XL_homologous_gene_name_table.txt', delim = '\t')

# for select genes with increased tail lengths, set the cutoff for filtering
tl_diff_cutoff = 15
non_na_cutoff = 3
df.filtered <- df.combined %>% filter(tl_min >= tl_diff_cutoff, non_na >= non_na_cutoff) %>%
	arrange(n_pred, desc(non_na), desc(tl_min)) %>% dplyr::rename(gene_name = hs_gene_name)

# save the table to a file
out_folder <- './'
write_delim(df.filtered, file = paste0(out_folder, 'HS_MM_XL_mRNA_TL_change_measurements_with_predictions_top_genes.txt'), delim = '\t')


###-----------------------------------------------------------------------------------------------
##-- This part generates changes in translational efficiency by gene 

# load human translational changes
hs.te <- read_delim('../Data/Human_oocytes/HS_oocytes_GV_MII_TE_changes_Zou_2022_Science.txt') %>%
	separate(gene_id, into = c('gene_id', NA), sep = '[.]')

hs.ctb <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
								filters = 'ensembl_gene_id',
								values = hs.te$gene_id,
								mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")) %>% as_tibble

hs.te <- left_join(hs.te, hs.ctb, by = c('gene_id' = 'ensembl_gene_id')) %>%
	dplyr::rename('gene_name' = 'hgnc_symbol') 

hs.te.sele <- hs.te %>% filter(!is.na(gene_name)) %>% dplyr::rename_with(~paste0('hs_', .x), everything())

# load mouse translational changes
mm.te <- read_delim('../Data/Mouse_oocytes/MM_oocytes_GV_MII_TE_changes_Xiong_2022_NCB.txt') %>%
	separate(gene_id, into = c('gene_id', NA), sep = '[.]')

# sometimes this part fails because the server is not availiable. Just re-run it.
mm.ctb <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
								filters = 'ensembl_gene_id',
								values = mm.te$gene_id,
								mart = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")) %>% as_tibble

mm.te <- left_join(mm.te, mm.ctb, by = c('gene_id' = 'ensembl_gene_id')) %>%
	dplyr::rename('gene_name' = 'external_gene_name')

mm.te.sele <- mm.te %>% filter(!is.na(gene_name)) %>% dplyr::rename_with(~paste0('mm_', .x), everything())

# examine data of homologous genes in both species
df.te.all <- hs.te.sele %>% full_join(., df.gn.all %>% na.omit %>% dplyr::select(-starts_with('xl_')), by = 'hs_gene_name') %>%
	left_join(., mm.te.sele, by = 'mm_gene_name') %>% 
		filter(!if_any(ends_with('_name'), is.na)) %>% # remove rows that have NA in gene names
			mutate(te_diff_min = reduce(across(ends_with('_te_change')), pmin, na.rm = T),
				 non_na = rowSums(!is.na(across(ends_with('_te_change'))))) # get the minimal TL change of all species

# some human genes may have homologous XL genes and/or MM genes, use the average
redundant_genes <- df.te.all %>% group_by(hs_gene_name) %>% summarize(n = n()) %>% filter(n > 1) %>% pull(hs_gene_name)
df.te.all.uniq <- lapply(redundant_genes, function(x){
	temp <- df.te.all %>% filter(hs_gene_name == x)
	temp_out <- temp %>% dplyr::slice(1)
	temp_out$mm_te_change <- ifelse(all(is.na(temp$mm_te_change)), NA, mean(temp$mm_te_change, na.rm = T))
	te_vec <- c(temp_out$hs_te_change, temp_out$mm_te_change)
	temp_out$te_change_min <- ifelse(all(is.na(te_vec)), NA, min(te_vec, na.rm = T))
	temp_out$non_na <- sum(!is.na(te_vec))
	return(temp_out)
}) %>% bind_rows(.) %>%
	bind_rows(., df.te.all %>% filter(!(hs_gene_name %in% redundant_genes)))


###-----------------------------------------------------------------------------------------------
##-- This part makes some heatmaps of tail-length changes and translational changes

#-----------------------------------------------------------
#- genes with all measured tail-length-change values (n=19)
df.sele <- df.combined %>% filter(tl_min >= 15, non_na == 3, n_pred == 0) %>%
	arrange(n_pred, desc(non_na), desc(tl_min)) %>% dplyr::rename(gene_name = hs_gene_name)

# heatmap for tail-length change
res <- plot_heatmaps(mat_in = df.sele %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
										 info_df = df.sele %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl),
										 bound = 75, mybreaks = seq(15,75,15), legend_title = 'Tail-length change (nt)',  
										 row_split = NULL, show_row_names = T,
										 show_column_names = T, column_split = NULL,
										 column_labels = c('H.s.', 'M.m.', 'X.l.'),
										 g_fsize = 4.5, anno_width = 0.1, 
										 f_type = 'pdf', f_width = 2, f_height = 2,
										 fn = 'HS_MM_XL_mRNA_TL_change_top_genes_all_measured')

# corresponding heatmap for translation change
df.plot <- df.te.all.uniq %>% filter(hs_gene_name %in% df.sele$gene_name) %>%
	mutate(gene_name = hs_gene_name) %>% dplyr::select(gene_name, hs_te_change, mm_te_change) %>%
		mutate(gene_name = factor(gene_name, levels = df.sele$gene_name)) %>% arrange(gene_name)
res <- plot_heatmaps(mat_in = df.plot %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
										 info_df = df.plot,
										 bound = 3, mybreaks = seq(-3,3,1),legend_title = 'Translational efficiency change (log2)',
										 row_split = NULL, show_row_names = T,
										 show_column_names = T, column_split = NULL,
										 column_labels = c('H.s.', 'M.m.'),
										 g_fsize = 4.5, anno_width = 0.1, main_fill = 'PiYG',
										 f_type = 'pdf', f_width = 2, f_height = 2,
										 fn = 'HS_MM_mRNA_TE_change_sele_genes_all_measured')

#-----------------------------------------------------------
#- genes with one tail-length-change values (57)
###--- Fig. 6d ---###

df.sele <- df.combined %>% filter(tl_min >= 15, non_na == 3, n_pred == 1) %>%
	arrange(n_pred, desc(non_na), desc(tl_min)) %>% dplyr::rename(gene_name = hs_gene_name) 

# heatmap for tail-length change
res <- plot_heatmaps(mat_in = df.sele %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
										 info_df = df.sele %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl),
										 pred_mat = df.sele %>% dplyr::select(gene_name, hs_flag, mm_flag, xl_flag) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
										 bound = 75, mybreaks = seq(15,75,15), legend_title = 'Tail-length change (nt)', 
										 row_split = NULL, show_row_names = T,
										 show_column_names = T, column_split = NULL, 
										 column_labels = c('H.s.', 'M.m.', 'X.l.'),
										 g_fsize = 4.5, anno_width = 0.1, 
										 f_type = 'pdf', f_width = 2, f_height = 4,
										 fn = 'HS_MM_XL_mRNA_TL_change_top_genes_one_predicted')

# corresponding heatmap for translation change
df.plot <- df.te.all.uniq %>% filter(hs_gene_name %in% df.sele$gene_name) %>%
	mutate(gene_name = hs_gene_name) %>% dplyr::select(gene_name, hs_te_change, mm_te_change) %>%
		mutate(gene_name = factor(gene_name, levels = df.sele$gene_name)) %>% arrange(gene_name)
res <- plot_heatmaps(mat_in = df.plot %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
										 info_df = df.plot,
										 bound = 3, mybreaks = seq(-3,3,1),legend_title = 'Translational efficiency change (log2)',
										 row_split = NULL, show_row_names = T,
										 show_column_names = T, column_split = NULL,
										 column_labels = c('H.s.', 'M.m.'),
										 g_fsize = 4.5, anno_width = 0.1, main_fill = 'PiYG',
										 f_type = 'pdf', f_width = 2, f_height = 4,
										 fn = 'HS_MM_mRNA_TE_change_sele_genes_one_predicted')

#-----------------------------------------------------------
#- genes with two tail-length-change values (108)
###--- Supplementary Fig. 8f ---###

df.sele <- df.combined %>% filter(tl_min >= 15, non_na == 3, n_pred == 2) %>%
	arrange(n_pred, desc(non_na), desc(tl_min)) %>% dplyr::rename(gene_name = hs_gene_name) 

# heatmap for tail-length change
temp_lst <- list(df.sele %>% slice(1:floor(nrow(df.sele)/2)), df.sele %>% slice((floor(nrow(df.sele)/2) + 1):nrow(df.sele)))
temp_vec <- c('part_1', 'part_2')
res <- lapply(1:2, function(x){
	plot_heatmaps(mat_in = temp_lst[[x]] %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
								info_df = temp_lst[[x]] %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl),
								pred_mat = temp_lst[[x]] %>% dplyr::select(gene_name, hs_flag, mm_flag, xl_flag) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
								bound = 75, mybreaks = seq(15,75,15), legend_title = 'Tail-length change (nt)',
								row_split = NULL, show_row_names = T,
								show_column_names = T, column_split = NULL, 
								column_labels = c('H.s.', 'M.m.', 'X.l.'),
								g_fsize = 4.5, anno_width = 0.1, 
								f_type = 'pdf', f_width = 2, f_height = 4,
								fn = paste0('HS_MM_XL_mRNA_TL_change_top_genes_two_predicted', temp_vec[[x]]))
})
	
# corresponding heatmap for translation change
df.plot.lst <- lapply(temp_lst, function(x){
	df.te.all.uniq %>% filter(hs_gene_name %in% x$gene_name) %>%
		mutate(gene_name = hs_gene_name) %>% dplyr::select(gene_name, hs_te_change, mm_te_change) %>%
			mutate(gene_name = factor(gene_name, levels = x$gene_name)) %>% arrange(gene_name)
})
	
res <- lapply(1:2, function(x){
	plot_heatmaps(mat_in = df.plot.lst[[x]] %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
								info_df = df.plot.lst[[x]],
								bound = 3, mybreaks = seq(-3,3,1),legend_title = 'Translational efficiency change (log2)',
								row_split = NULL, show_row_names = T,
								show_column_names = T, column_split = NULL,
								column_labels = c('H.s.', 'M.m.'),
								g_fsize = 4.5, anno_width = 0.1, main_fill = 'PiYG',
								f_type = 'pdf', f_width = 2, f_height = 4,
								fn = paste0('HS_MM_mRNA_TE_change_sele_genes_two_predicted', temp_vec[[x]]))
})

#-----------------------------------------------------------
#- genes with three tail-length-change values (99)
###--- Supplementary Fig. 8g ---###

df.sele <- df.combined %>% filter(tl_min >= 15, non_na == 3, n_pred == 3) %>%
	arrange(n_pred, desc(non_na), desc(tl_min)) %>% dplyr::rename(gene_name = hs_gene_name) 

# heatmap for tail-length change
temp_lst <- list(df.sele %>% slice(1:floor(nrow(df.sele)/2)), df.sele %>% slice((floor(nrow(df.sele)/2) + 1):nrow(df.sele)))
temp_vec <- c('part_1', 'part_2')
res <- lapply(1:2, function(x){
	plot_heatmaps(mat_in = temp_lst[[x]] %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
								info_df = temp_lst[[x]] %>% dplyr::select(gene_name, hs_tl, mm_tl, xl_tl),
								pred_mat = temp_lst[[x]] %>% dplyr::select(gene_name, hs_flag, mm_flag, xl_flag) %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
								bound = 75, mybreaks = seq(15,75,15), legend_title = 'Tail-length change (nt)',
								row_split = NULL, show_row_names = T,
								show_column_names = T, column_split = NULL, 
								column_labels = c('H.s.', 'M.m.', 'X.l.'),
								g_fsize = 4.5, anno_width = 0.1, 
								f_type = 'pdf', f_width = 2, f_height = 4,
								fn = paste0('HS_MM_XL_mRNA_TL_change_top_genes_three_predicted', temp_vec[[x]]))
})

# corresponding heatmap for translation change
df.plot.lst <- lapply(temp_lst, function(x){
	df.te.all.uniq %>% filter(hs_gene_name %in% x$gene_name) %>%
		mutate(gene_name = hs_gene_name) %>% dplyr::select(gene_name, hs_te_change, mm_te_change) %>%
		mutate(gene_name = factor(gene_name, levels = x$gene_name)) %>% arrange(gene_name)
})

res <- lapply(1:2, function(x){
	plot_heatmaps(mat_in = df.plot.lst[[x]] %>% column_to_rownames(var = 'gene_name') %>% as.matrix(.),
								info_df = df.plot.lst[[x]],
								bound = 3, mybreaks = seq(-3,3,1),legend_title = 'Translational efficiency change (log2)',
								row_split = NULL, show_row_names = T,
								show_column_names = T, column_split = NULL,
								column_labels = c('H.s.', 'M.m.'),
								g_fsize = 4.5, anno_width = 0.1, main_fill = 'PiYG',
								f_type = 'pdf', f_width = 2, f_height = 4,
								fn = paste0('HS_MM_mRNA_TE_change_sele_genes_three_predicted', temp_vec[[x]]))
})


###-----------------------------------------------------------------------------------------------
##-- This part performs gene set enrichment analysis
###--- Fig. 6e ---###

df.sele <- df.combined %>% filter(tl_min >= 15, non_na == 3) %>% arrange(desc(hs_tl))

gres <- gost(query = df.sele %>% pull(hs_gene_name),
						 organism = "hsapiens", ordered_query = TRUE,
						 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
						 measure_underrepresentation = FALSE, evcodes = FALSE,
						 user_threshold = 0.05, correction_method = "g_SCS",
						 domain_scope = "annotated", custom_bg = NULL,
						 numeric_ns = "", sources = NULL, as_short_link = FALSE)

# get top non-overlapping GO categories
n_cutoff = 8
sele_go = c()
sele_vec = c()
for (i in 1:nrow(gres$result)){
	if (gres$result$source[i] == 'GO:BP'){
		id_all = c(gres$result$term_id[i], gres$result$parents[[i]])
		if (!any(id_all %in% sele_go)){
			sele_vec <- c(sele_vec, i)
			if (length(sele_vec) >= n_cutoff){
				break
			}
		}
		sele_go <- c(sele_go, gres$result$term_id[i])
	}
}

gost_plot <- gost_bar_plot(df_in = gres$result %>% dplyr::slice(sele_vec), padj_cutoff = 0.05, 
													 mycolor = alpha(brewer.pal(8,'Accent')[1], 1), 
													 mybreaks = seq(0,15,1),
													 fn = paste0('HS_MM_XL_mRNA_TL_change_top_genes_gprofiler_sele'))

