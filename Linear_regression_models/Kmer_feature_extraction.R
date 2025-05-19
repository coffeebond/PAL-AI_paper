library(tidyverse)
library(Biostrings)
library(parallel)
library(RcppAlgos)
n_core = 60

kmer_len_min = 3
kmer_len_max = 7

fa <- readDNAStringSet('../Data/Linear_Regression/XENLA_10.1_GCF_XBmodels_annotation_fixed_UTR3_polished_by_PA_site_stitched.fa')
len_max <- 300
len_min <- 10
seq_append <- 'AAA'

# append 3A after each 3'-UTR
fa_3a <- xscat(fa, seq_append) 
names(fa_3a) <- names(fa)

# filter by minimal length
fa_filtered <- fa_3a[width(fa_3a) >= len_min]

utr3 <- DNAStringSet(fa_filtered, start = pmax(width(fa_filtered) - len_max + 1, 1), end = width(fa_filtered))
utr3_id <- str_split_fixed(names(utr3), '__',2)[,1]

# kmer count
df.num <- lapply(kmer_len_min:kmer_len_max, function(x){
	as_tibble(oligonucleotideFrequency(utr3, x, as.prob = F))
}) %>% bind_cols %>% 
	rename_with(function(x){paste0('num_', x)}) %>%
		mutate(id = utr3_id) %>% relocate(id, .before = everything())
	
# output the kmer count table 
write_delim(df.num,	file = paste0('XL_v10_UTR3_last_', len_max, '_nt_kmer_',kmer_len_min, '_to_', kmer_len_max,  '_count_table.gz'), delim = '\t')

# kmer distance (this may take a while)
kmer_lst <- lapply(kmer_len_min:kmer_len_max, function(x){
	apply(permuteGeneral(c('A', 'C', 'G', 'T'), m = x, repetition = T), 1, paste, collapse = '')
}) %>% Reduce(function(x,y){c(x,y)}, .)
	
df.pos <- mclapply(kmer_lst, function(x){
	# find all kmer positions and their distances to 3' ends
	# use "reverse" to find the distance
	v <- sapply(vmatchPattern(reverse(x), reverse(utr3)), function(y){
			if (length(y)==0){
				return(0)
			} else {
				d <- pmax(start(y) - 1 - length(seq_append), 1)
				return(sum(1/d))
			}
		})
	v <- list(v)
	names(v) <- x
	return(v)
}, mc.cores = min(n_core, length(kmer_lst))) %>% bind_cols %>% 
	rename_with(function(x){paste0('pos_', x)}) %>%
		mutate(id = utr3_id) %>% relocate(id, .before = everything())

# output the kmer count table 
write_delim(df.pos,	file = paste0('XL_v10_UTR3_last_', len_max, '_nt_kmer_',kmer_len_min, '_to_', kmer_len_max,  '_position_table.gz'), delim = '\t')


