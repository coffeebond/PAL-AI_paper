# this script calculate translational efficiency changes for human and mouse mRNAs during oocyte maturation with re-precessed published datasets
library(tidyverse)
library(DESeq2)

size.factor = function(dframe) {
	# a function for normalization using DESeq2
	if (is.data.frame(dframe) == F) {
		stop("Please use dataframe as input for this function!")
	} else {
		#ref = apply(dframe,1,mean)
		#ratio = apply(dframe[ref!=0,],2,function(x){x/ref[ref!=0]})
		#sf = apply(ratio,2,median)
		#sf = sf/max(sf)
		suppressWarnings(suppressMessages(library(DESeq2)))
		colData=as.data.frame(colnames(dframe))
		colnames(colData)=c('condition')
		rownames(colData)=colData[,1]
		dds = DESeqDataSetFromMatrix(countData = dframe, colData=colData, design = ~condition)
		dds = estimateSizeFactors(dds)
		sf = sizeFactors(dds)
	}
	return (sf)
}

###-----------------------------------------------------------------------------------------------------------------------------------------------
##- This part is for human

# load human data
info_hs <- read.table('../Data/Human_oocytes/info_translational_efficiency.txt', header = T)
for (i in 1:nrow(info_hs)){
	temp <- read.table(paste0(info_hs$Folder[i], info_hs$File[i])) %>% as_tibble
	colnames(temp) <- c('gene_id', info_hs$Name[i])
	temp <- temp %>% filter(!grepl('^__', gene_id)) 
	if (info_hs$Type[i] == 'RNA'){
		if (exists('df.hs.rna')){
			if (colnames(temp[2]) %in% colnames(df.hs.rna)){
				# merge technical replicates
				df.hs.rna <- inner_join(df.hs.rna, temp, by = 'gene_id', suffix = c('_TR_1', '_TR_2'))
				df.hs.rna$TR_merged <- df.hs.rna[, paste0(colnames(temp[2]), '_TR_1')] + df.hs.rna[, paste0(colnames(temp[2]), '_TR_2')]
				df.hs.rna[, colnames(temp[2])]	<- df.hs.rna$TR_merged
				df.hs.rna <- df.hs.rna %>% select(-contains('TR_'))
			} else {
				df.hs.rna <- inner_join(df.hs.rna, temp, by = 'gene_id')
			}
		} else {
			df.hs.rna <- temp 
		}
	}
	if (info_hs$Type[i] == 'RPF'){
		if (exists('df.hs.rpf')){
			if (colnames(temp[2]) %in% colnames(df.hs.rpf)){
				# merge technical replicates
				df.hs.rpf <- inner_join(df.hs.rpf, temp, by = 'gene_id', suffix = c('_TR_1', '_TR_2'))
				df.hs.rpf$TR_merged <- df.hs.rpf[, paste0(colnames(temp[2]), '_TR_1')] + df.hs.rpf[, paste0(colnames(temp[2]), '_TR_2')]
				df.hs.rpf[, colnames(temp[2])]	<- df.hs.rpf$TR_merged
				df.hs.rpf <- df.hs.rpf %>% select(-contains('TR_'))
			} else {
				df.hs.rpf <- inner_join(df.hs.rpf, temp, by = 'gene_id')
			}
		} else {
			df.hs.rpf <- temp 
		}
	}
}

###--- calculate TE and TE changes
# GV stage 
rna_cutoff <- 30
df.hs.rna.sub <- df.hs.rna %>% select(contains('GV'), gene_id) %>% filter(if_all(contains('RNA'), function(x){x >= rna_cutoff}))
df.hs.all <- inner_join(df.hs.rna.sub, df.hs.rpf %>% select(contains('GV'), gene_id), by = 'gene_id') %>% as.data.frame(.)
rownames(df.hs.all) <- df.hs.all$gene_id
df.hs.all$gene_id <- NULL
sf <- size.factor(df.hs.all)
sf <- sf / max(sf)
df.hs.all.n <- t(t(df.hs.all) / sf) %>% as.data.frame(.)
df.hs.te.gv <- df.hs.all.n %>% mutate(TE_GV = 0.5 * log2((GV_RPF_rep1+1)*(GV_RPF_rep2+1)/((GV_RNA_rep1+1)*(GV_RNA_rep2+1))))
write_delim(df.hs.te.gv %>% as_tibble(rownames = 'gene_id'), file = 'HS_oocytes_GV_TE_Zou_2022_Science.txt', delim = '\t')

# MII stage 
rna_cutoff <- 30
df.hs.rna.sub <- df.hs.rna %>% select(contains('MII'), gene_id) %>% filter(if_all(contains('RNA'), function(x){x >= rna_cutoff}))
df.hs.all <- inner_join(df.hs.rna.sub, df.hs.rpf %>% select(contains('MII'), gene_id), by = 'gene_id') %>% as.data.frame(.)
rownames(df.hs.all) <- df.hs.all$gene_id
df.hs.all$gene_id <- NULL
sf <- size.factor(df.hs.all)
sf <- sf / max(sf)
df.hs.all.n <- t(t(df.hs.all) / sf) %>% as.data.frame(.)
df.hs.te.m2 <- df.hs.all.n %>% mutate(TE_MII = 0.5 * log2((MII_RPF_rep1+1)*(MII_RPF_rep2+1)/((MII_RNA_rep1+1)*(MII_RNA_rep2+1))))
write_delim(df.hs.te.m2 %>% as_tibble(rownames = 'gene_id'), file = 'HS_oocytes_MII_TE_Zou_2022_Science.txt', delim = '\t')

# calculate TE changes
df.hs.te <- inner_join(df.hs.te.gv %>% as_tibble(rownames = 'gene_id') %>% select(gene_id, TE_GV), 
										df.hs.te.m2 %>% as_tibble(rownames = 'gene_id') %>% select(gene_id, TE_MII), by = 'gene_id') %>%
	mutate(te_change = TE_MII - TE_GV)
write_delim(df.hs.te, file = 'HS_oocytes_GV_MII_TE_changes_Zou_2022_Science.txt', delim = '\t')

###-----------------------------------------------------------------------------------------------------------------------------------------------
##- This part is for mouse

# load mouse data
info_mm <- read.table('../Data/Mouse_oocytes/info_translational_efficiency.txt', header = T)
for (i in 1:nrow(info_mm)){
	temp <- read.table(paste0(info_mm$Folder[i], info_mm$File[i])) %>% as_tibble
	colnames(temp) <- c('gene_id', info_mm$Name[i])
	temp <- temp %>% filter(!grepl('^__', gene_id)) 
	if (info_mm$Type[i] == 'RNA'){
		if (exists('df.mm.rna')){
			df.mm.rna <- inner_join(df.mm.rna, temp, by = 'gene_id')
		} else {
			df.mm.rna <- temp 
		}
	}
	if (info_mm$Type[i] == 'RPF'){
		if (exists('df.mm.rpf')){
			df.mm.rpf <- inner_join(df.mm.rpf, temp, by = 'gene_id')
		} else {
			df.mm.rpf <- temp 
		}
	}
}

###--- calculate TE and TE changes
# GV stage 
rna_cutoff <- 30
df.mm.rna.sub <- df.mm.rna %>% select(contains('GV'), gene_id) %>% filter(if_all(contains('RNA'), function(x){x >= rna_cutoff}))
df.mm.all <- inner_join(df.mm.rna.sub, df.mm.rpf %>% select(contains('GV'), gene_id), by = 'gene_id') %>% as.data.frame(.)
rownames(df.mm.all) <- df.mm.all$gene_id
df.mm.all$gene_id <- NULL
sf <- size.factor(df.mm.all)
sf <- sf / max(sf)
df.mm.all.n <- t(t(df.mm.all) / sf) %>% as.data.frame(.)
df.mm.te.gv <- df.mm.all.n %>% mutate(TE_GV = 0.5 * log2((GV_RPF_rep1+1)*(GV_RPF_rep2+1)/((GV_RNA_rep1+1)*(GV_RNA_rep2+1))))
write_delim(df.mm.te.gv %>% as_tibble(rownames = 'gene_id'), file = 'MM_oocytes_GV_TE_Xiong_2022_NCB.txt', delim = '\t')

# MII stage 
rna_cutoff <- 30
df.mm.rna.sub <- df.mm.rna %>% select(contains('MII'), gene_id) %>% filter(if_all(contains('RNA'), function(x){x >= rna_cutoff}))
df.mm.all <- inner_join(df.mm.rna.sub, df.mm.rpf %>% select(contains('MII'), gene_id), by = 'gene_id') %>% as.data.frame(.)
rownames(df.mm.all) <- df.mm.all$gene_id
df.mm.all$gene_id <- NULL
sf <- size.factor(df.mm.all)
sf <- sf / max(sf)
df.mm.all.n <- t(t(df.mm.all) / sf) %>% as.data.frame(.)
df.mm.te.m2 <- df.mm.all.n %>% mutate(TE_MII = 0.5 * log2((MII_RPF_rep1+1)*(MII_RPF_rep2+1)/((MII_RNA_rep1+1)*(MII_RNA_rep2+1))))
write_delim(df.mm.te.m2 %>% as_tibble(rownames = 'gene_id'), file = 'MM_oocytes_MII_TE_Xiong_2022_NCB.txt', delim = '\t')

# calculate TE changes
df.mm.te <- inner_join(df.mm.te.gv %>% as_tibble(rownames = 'gene_id') %>% select(gene_id, TE_GV), 
										df.mm.te.m2 %>% as_tibble(rownames = 'gene_id') %>% select(gene_id, TE_MII), by = 'gene_id') %>%
	mutate(te_change = TE_MII - TE_GV)
write_delim(df.mm.te, file = 'MM_oocytes_GV_MII_TE_changes_Xiong_2022_NCB.txt', delim = '\t')

