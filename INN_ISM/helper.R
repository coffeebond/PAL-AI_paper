library('signs')
library('Cairo')
library('RColorBrewer')
library('ggrepel')
library('ggseqlogo')
library('cowplot')

custom_breaks <- function(bmin, bmax, digits = 0, length.out = 8, zero = TRUE) {
	bmin = floor(bmin * (10 ^ digits)) / (10 ^ digits)
	bmax = ceiling(bmax * (10 ^ digits)) / (10 ^ digits)
	if (bmin > 0 | bmax < 0) zero = FALSE
	d = round((bmax - bmin) / (length.out - 1), digits)
	if (d == 0) d = round((bmax - bmin) / (length.out - 1), digits+1)
	if (zero) {
		return(unique(c(rev(seq(0, bmin, -d)), seq(0, bmax, d))))
	} else {
		return(seq(bmin, bmax, d))
	}
}

custom_quantile <- function(x) {
	r <- c(quantile(x, 0.1), quantile(x, 0.25), median(x), quantile(x, 0.75), quantile(x, 0.9))
	names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
	r
}

subset_DNAStringSet <- function(s, l){
	# a function to subset a DNAStringSet to the last l nucleotides
	# s is a DNAStringSet
	# l is the length of 3'-end region to keep
	DNAStringSet(s, start = pmax(1, width(s) - l + 1), end = width(s))
}

encode_seq_vec <- function(seq_vec){
	# one-hot encode sequence
	sapply(seq_vec, function(x){
		if (x == 'A'){
			return(c(1,0,0,0))
		} else if (x == 'C') {
			return(c(0,1,0,0))
		} else if (x == 'G') {
			return(c(0,0,1,0))
		} else if (x == 'T' | x == 'U') {
			return(c(0,0,0,1))
		} else return(c(0,0,0,0))
	}) %>% cbind(.)
}

stat_matrix_to_kmer_df <- function(m, nt_vec = c('A', 'C', 'G', 'U')){
	# this converts a mutation statistic matrix to a dataframe by kmer
	wt_kmer = sapply(1:ncol(m), function(x){
		nt_vec[is.na(m[,x])]
	}) %>% paste(collapse ='')
	
	lapply(1:ncol(m), function(x){
		lapply(c(1:nrow(m))[!is.na(m[,x])], function(y){
			mt_kmer <- wt_kmer
			substr(mt_kmer, x, x) <- nt_vec[y] 
			list('kmer' = mt_kmer, 'val' = m[y,x])
		}) %>% bind_rows()
	}) %>% bind_rows()
}

one_sample_t_test <- function(m, sd, count, diff = 0, lower.tail = TRUE) {
	# perform one-sample t-test against a set value 
	# null hypothesis: m >= diff
	t <- (m - diff) / (sd / (count ** 0.5)) 
	pval <- pt(t, count-1, lower.tail = lower.tail)
	#return(list('t' = t, 'pval' = pval))
	return(pval)
}

tl_change_vs_count_splot <- function(df_in, kmer_len = 5, count_cutoff = 1, label_sele = NULL, color_sele = NULL, ybreaks = NULL, fn = 'plot'){
	# This function takes in a position-less kmer dataset makes a scatter plot comparing the mean tail length change (avg) to count number (n)
	df_in <- df_in %>% filter(nchar(kmer) == kmer_len & n >= count_cutoff)
	avg_mean <- mean(df_in$avg)
	avg_std <- sd(df_in$avg)
	df_in <- df_in %>% mutate(z_score = (avg - avg_mean) / avg_std) %>%
		arrange(desc(z_score)) %>% mutate(rank = 1:nrow(.))
	
	if (length(color_sele) != length(label_sele)){
		mycolor = c(brewer.pal(9, 'Set1')[1:length(label_sele)], t_col('gray50', 0.3))
	} else {
		mycolor = c(color_sele, t_col('gray50', 0.3))
	}
	df_in$color <- mycolor[length(mycolor)]
	
	if (is.list(label_sele)){ # labels must be a list of vectors
		for (i in 1:length(label_sele)){
			df_in <- df_in %>% mutate(color = ifelse(kmer %in% label_sele[[i]], mycolor[i], color))
		}
		df_in <- df_in %>% mutate(label = ifelse(kmer %in% unlist(label_sele), gsub('T', 'U', kmer), ''))
	} else {
		stop('"label_sele" must be a list!')
	}
	
	df_in$color <- factor(df_in$color, levels = mycolor)
	df_in <- df_in %>% arrange(desc(color))
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	xmin = min(df_in$n)
	xmax = max(df_in$n)
	xbreaks = 10 ** seq(floor(log10(xmin)), ceiling(log10(xmax)),1)
	xlabels = paste0('10^', log10(xbreaks))
	if(is.null(ybreaks)){
		ymin <- min(df_in$avg)
		ymax <- max(df_in$avg)
		ybreaks = custom_breaks(ymin, ymax, digits = 1, length.out = 7)
	} else {
		ymin <- min(df_in$avg, ybreaks)
		ymax <- max(df_in$avg, ybreaks)
	}
	
	z_min = min(min(df_in$z_score), min((ybreaks - avg_mean)/avg_std))
	z_max = max(max(df_in$z_score), max((ybreaks - avg_mean)/avg_std))
	z_label = unique(c(rev(seq(0, z_min, by=-3)), 0, seq(0, z_max, by=3)))
	n_kmer = nrow(df_in)
	
	gplot <- ggplot(df_in, aes(x = n, y = avg)) +
		geom_point(size = 1, aes(fill = color), shape = 21, color = t_col('black', 0.5), stroke = 0.2) + 
		labs(x = 'Number of sequences', y = "Mean predicted tail-length difference (nt)") + 
		coord_cartesian(xlim = c(min(xbreaks), max(xbreaks)), ylim = c(ymin, ymax)) +
		scale_x_continuous(trans = 'log10', breaks = xbreaks, labels = do.call(expression, rlang::parse_exprs(xlabels))) +
		scale_y_continuous(breaks = ybreaks, labels = signs(ybreaks, accuracy = 0.1),
											 sec.axis = sec_axis(~ ((. - avg_mean) / avg_std), name = 'Z-score', 
											 										breaks = z_label, labels = signs(z_label))) +
		scale_fill_manual(values = mycolor) + 
		geom_hline(yintercept = avg_mean, linetype='dashed', size = line_w, color = t_col('black', 0.5)) +
		annotate('text', x = xbreaks[1], y = ybreaks[length(ybreaks)], hjust = 0, vjust = 1, label = paste0('n = ', signs(n_kmer, format = scales::comma)), size = (font_size - font_size_offset) * 5 /14) + 
		theme(
			legend.position='none',
			legend.justification = c(0,1),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.title = element_text(size=font_size,color='black', hjust = 0.5, face = 'bold'),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			
			axis.title.y = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.y = element_text(size=font_size,color='black', hjust = 1, margin = margin(0,1,0,0)),
			axis.line.y = element_line(color='black',size=line_w),
			axis.ticks.y = element_line(color='black',size=line_w),
			axis.text.y.right = element_text(hjust = 1, margin = margin(0,0,0,1)),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5, margin = margin(0,0,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin = margin(0,0,0,0)),
			axis.line.x = element_line(color='black',size=line_w),
			axis.ticks.x = element_line(color='black',size=line_w),
			axis.ticks.length = unit(2, "pt"),
			
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(font_size - font_size_offset, 'pt'),
			legend.text = element_text(size = font_size - font_size_offset, margin = margin(0,0,1,-3)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(-3,0,0,3),
			aspect.ratio = 1
		)
	if(!is.null(label_sele)){
		gplot <- gplot + geom_text_repel(data = df_in, aes(x = n, y = avg, label = label), max.overlaps = Inf,
																		 color = 'black', size = (font_size - font_size_offset) * 5 / 14, seed = 6,
																		 min.segment.length = 0, segment.size = 0.5/.pt, box.padding = 0.8)
	}
	#png(file=paste0(fn, '.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	CairoPDF(paste0(fn, '.pdf'), width = 2.3, height = 2, family = 'Arial')
	print(gplot)
	dev.off()
	invisible(df_in)
}

seq_align <- function(ref_string, ali_string, score_match = 1, score_offset = -1, score_mismatch = -1.5) {
	# a function to manually align two strings by scanning all offset positions with no gap allowed
	# returns the highest score and corresponding offset position
	# if there is a tie, output all offset positions
	ref = paste0(strrep('N', (nchar(ali_string)-1)), ref_string, strrep('N', (nchar(ali_string)-1)))
	s.df = data.frame('score' = -Inf, 'offset' = (-(nchar(ali_string)-1)):(nchar(ref_string)-1))
	for (i in 1:(nchar(ref) - nchar(ali_string) + 1)) {
		ref.sub <- substr(ref, i, (i + nchar(ali_string) - 1))
		ref.vec <- unlist(strsplit(ref.sub, ''))
		ali.vec <- unlist(strsplit(ali_string, ''))
		s.df$score[i] <- sapply(1:length(ref.vec), function(x) {
			if (ref.vec[x] == ali.vec[x]) {
				return(score_match)
			} else if (ref.vec[x] == 'N') {
				return(score_offset)
			} else return(score_mismatch)
		}) %>% sum(.)
	}
	max_score <- max(s.df$score)
	offset_pos <- s.df[s.df$score == max_score, 'offset']
	return(list('max_score' = max_score, 'offset_pos' = offset_pos))
}

plot_top_kmer_logo <- function(df_in, weight = 'diff_avg', mode = 'positive', z_cutoff = 3,
															 align_score_cutoff = -0.1, score_match = 1, score_offset = -1, score_mismatch = -1.5, n_seed_max = 1, type = 'bits', ft = 'pdf', fn = 'plot') {
	# this function aligns kmers passing z_cutoff with my custom method and makes a logo plot with specified weight
	# weight: column name of selection of weight for ranking and computing weight matrix, can either be 't' or 'diff_avg'
	# mode: to examine positive or negative tail length changes, must be either 'positive' or 'negative'
	# types: can be either 'bits' or 'probability' (from 'ggseqlogo')
	
	# make sure 'z_score' column is the input dataframe
	df_in <- as.data.frame(df_in)
	if ('z_score' %in% colnames(df_in) == FALSE) {
		warning('No z_score in the columns! \"z_cutoff\" is not applied.')
	} else {
		# subset the table
		if (z_cutoff > 0 ) {
			df.sub <- df_in[df_in$z_score > z_cutoff, ]
		} else if (z_cutoff < 0) {
			df.sub <- df_in[df_in$z_score < z_cutoff, ]
		} else print('\"z_cutoff\" must be either positive or negtive values!')
	}
	
	if (mode == 'positive') {
		df.sub <- df.sub[df.sub[, weight] > 0, ]
	} else if (mode == 'negative') {
		df.sub <- df.sub[df.sub[, weight] < 0, ]
	} else print('\"diff_group\" must be either \"positive\" or \"negtive\"!')
	
	if(nrow(df.sub) == 0) {
		stop('No kmers meets the required condition!')
	}
	
	# sort the dataframe in case it's not pre-sorted (based on weight)
	df.sub <- df.sub[order(abs(df.sub[ ,weight]), decreasing = TRUE), ]
	
	print(paste0('Start alignment of kmers passing the cutoff, with a total number of ', nrow(df.sub)))
	
	# use my custom way to aligning kmers
	aln.lst <- list() # the list for aligned kmers, each seed has a dataframe
	seeds <- c() # all the seeds that kmers aligned to
	for (i in 1:nrow(df.sub)) {
		if (i == 1) {
			# initialize the first seed 
			seeds[df.sub$kmer[i]] <- df.sub$kmer[i]
			aln.lst[[df.sub$kmer[i]]] <- data.frame('kmer' = df.sub$kmer[i],
																							'kmer_align' = paste0(strrep('-', nchar(df.sub$kmer[i])-1), df.sub$kmer[i], strrep('-', nchar(df.sub$kmer[i])-1)),
																							'weight' = df.sub[, weight][i],
																							'kmer_count' = 1,
																							stringsAsFactors = F)
		} else {
			aln.scores <- sapply(1:length(seeds), function(x){
				seq_align(ref_string = seeds[x], ali_string = df.sub$kmer[i], score_match = score_match, score_offset = score_offset, score_mismatch = score_mismatch)[[1]]
			})
			
			if (max(aln.scores) > align_score_cutoff | length(seeds) >= n_seed_max ) {	
				# equally distribute the 'weight' to best matched offset positions for aligned kmer seeds
				for (j in 1:length(aln.scores)) {
					if (aln.scores[j] == max(aln.scores)) {
						offset_pos <- seq_align(ref_string = seeds[j], ali_string = df.sub$kmer[i], score_match = score_match, score_offset = score_offset, score_mismatch = score_mismatch)[[2]]
						partial_weight <- df.sub[, weight][i] / sum(aln.scores == max(aln.scores)) / length(offset_pos)
						for (pos in offset_pos) {
							aln.lst[[j]][nrow(aln.lst[[j]])+1, ] <- list('kmer' = df.sub$kmer[i],
																													 'kmer_align' = paste0(strrep('-', nchar(df.sub$kmer[i])-1+pos), df.sub$kmer[i], strrep('-', nchar(df.sub$kmer[i])-1-pos)),
																													 'weight'= partial_weight,
																													 'kmer_count' = 1 / sum(aln.scores == max(aln.scores)) / length(offset_pos))
						}
					}
				}
			} else {
				# if the kmer cannot be aligned to any kmer seeds, ceate a new seed 
				seeds[df.sub$kmer[i]] <- df.sub$kmer[i]
				aln.lst[[df.sub$kmer[i]]] <- data.frame('kmer' = df.sub$kmer[i],
																								'kmer_align' = paste0(strrep('-', nchar(df.sub$kmer[i])-1), df.sub$kmer[i], strrep('-', nchar(df.sub$kmer[i])-1)),
																								'weight' = df.sub[, weight][i], 
																								'kmer_count' = 1,
																								stringsAsFactors = F)
			}
		}
	}
	# combine the list into a dataframe
	aln.df <- lapply(1:length(aln.lst), function(x) {
		aln.lst[[x]] %>% mutate(., seed = names(aln.lst)[x])
	}) %>% bind_rows(.)
	aln.df$seed <- factor(aln.df$seed, levels = names(aln.lst))
	
	# calculate nucleotide weight matrices for each position of each aligned kmer group
	weight.lst <- lapply(1:length(aln.lst), function(x) {
		temp.lst <- lapply(1:nrow(aln.lst[[x]]), function(y) {
			temp.df <- matrix(0, nrow = 5, ncol = nchar(aln.lst[[x]]$kmer_align[y]))
			rownames(temp.df) <- c('A', 'C', 'G', 'T', 'N')
			for (i in 1:ncol(temp.df)) {
				temp.df[gsub('-', 'N', unlist(str_split(aln.lst[[x]]$kmer_align[y], '')))[i], i] <- aln.lst[[x]]$weight[y]
			}
			return(temp.df)
		})
		
		# sum up all the weights for each base at each position
		weight.df <- Reduce('+', temp.lst)
		
		# remove positions that all weights are in base 'N' (no kmer aligned to these positions)
		# it assumes that this only happens on two ends
		weight.df <- weight.df[, weight.df['N', ] != apply(weight.df, 2, sum)]
		
		# equally distribute weights from 'N' to four bases and normalize for each position
		weight.df <- apply(weight.df, 2, function(x) {
			temp <- x[1:4] + x['N']/4
			return(temp/sum(temp))
		}) 
		rownames(weight.df)[rownames(weight.df) == 'T'] = 'U'
		
		return(weight.df)
	})
	names(weight.lst) <- names(aln.lst)
	
	# make logo plots with custom color scheme
	cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'U'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
												cols=brewer.pal(9,'Set1')[c(1:3,5)])
	line_w = 0.75/.pt
	font_size = 5.5
	for (i in 1:length(weight.lst)){
		file_name = paste0(fn, '_weight_', weight, '_', mode, '_top_', nrow(df.sub), '_kmer_logo_', i)
		gplot <- ggplot()+
			geom_logo(weight.lst[[i]], method = type, seq_type = 'rna', col_scheme=cs1)+
			theme(
				legend.position="none",
				axis.title.y = element_text(size=font_size,color='black',vjust=1),
				axis.text.y = element_text(size=font_size,color='black'),
				axis.line.y = element_line(color='black',size=line_w),
				axis.ticks.y = element_line(color='black',size=line_w),
				
				axis.title.x = element_text(size=font_size,color='black',vjust=1),
				axis.text.x = element_text(size=font_size,color='black',vjust=1),
				axis.line.x = element_blank(),
				axis.ticks.x = element_line(color='black',size=line_w),
				
				plot.background = element_rect(fill = 'transparent',color=NA),
				panel.background = element_blank(),
				panel.grid = element_blank()
			)
		if (ft != 'none' & ft != 'None'){
			if (ft == 'pdf') {
				CairoPDF(file=file_name, width=0.8*ncol(weight.lst[[i]]), height=3, family = 'Arial')
			} else if (ft == 'png') {
				png(file=file_name, width=0.8*ncol(weight.lst[[i]]), height=3,unit='in',res=600,bg='transparent')
			}
			print(gplot)
			dev.off()
		}
	}
	invisible(list('aln_lst' = aln.lst, 'aln_df' = aln.df,'weight_lst' = weight.lst))
}

plot_aligned_logos_combined <- function(aln_lst, aln_df, weight_lst, type = 'bits', ylab = 'Z-score', f_width = NULL, f_height = NULL, fn = 'logo'){
	# This function takes in output from 'plot.top.kmer.logo' and it makes a combined plot including:
	# 1) logo plots; 2) barplot; 3) pie charts
	
	# calculate weighted average of tail-length difference for each group
	weight.mean.df <- data.frame('avg' = sapply(1:length(aln_lst), function(x) {
		sum(aln_lst[[x]]$weight * abs(aln_lst[[x]]$weight) / sum(abs(aln_lst[[x]]$weight)))
	}),
	'seed' = names(aln_lst), stringsAsFactors = F)
	weight.mean.df$seed <- factor(weight.mean.df$seed, levels = rev(weight.mean.df$seed))
	
	# make logo plots with custom color scheme
	cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'U'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
												cols=brewer.pal(9,'Set1')[c(1:3,5)])
	gplot.logo.lst <- list()
	for (i in 1:length(weight_lst)) {
		gplot.logo.lst[[i]] <- ggplot()+
			geom_logo(weight_lst[[i]], method = type, seq_type = 'rna', col_scheme=cs1)+
			scale_y_continuous(limits = c(0,2))+
			theme(
				legend.position="none",
				axis.title.y = element_blank(),#element_text(size=12,color='black',vjust=1),
				axis.text.y = element_blank(),#element_text(size=12,color='black'),
				axis.line.y = element_blank(),#element_line(color='black',size=0.5),
				axis.ticks.y = element_blank(),#element_line(color='black',size=0.5),
				
				axis.title.x = element_blank(),#element_text(size=12,color='black',vjust=1),
				axis.text.x = element_blank(),# element_text(size=12,color='black',vjust=1),
				axis.line.x = element_blank(),
				axis.ticks.x = element_blank(),#element_line(color='black',size=0.5),
				
				plot.title = element_blank(),#element_text(color='black',size=14,hjust=0.5),
				plot.background = element_rect(fill = 'transparent',color=NA),
				panel.background = element_blank(),
				plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'inch'),
				panel.grid = element_blank()
			)
	}
	
	# make a barplot
	if (aln_df[1,'weight'] > 0) {
		mycolor = brewer.pal(11, 'PuOr')[9]
		xpos = 'bottom'
		ypos = 'right'
		ymin = 0
		ymax = 1.05 * max(aln_df$weight) #ceiling(max(weight.mean.df$avg))
		ybreaks = custom_breaks(ymin, ymax, length.out=6)
	} else {
		mycolor = brewer.pal(8, 'Set2')[1]
		xpos = 'top'
		ypos = 'right'
		ymin = 1.05 * min(aln_df$weight) #floor(min(weight.mean.df$avg))
		ymax = 0
		ybreaks = custom_breaks(ymin, ymax, length.out=6)
	}
	
	# correct weights for multi-aligned kmers 
	aln_df <- aln_df %>% select(!kmer_align) %>% group_by(kmer, seed) %>% 
		summarize(weight = sum(weight), kmer_count = sum(kmer_count), .groups = 'drop') %>%
		mutate(weight = ifelse(kmer_count != 1, weight / kmer_count, weight))
	
	line_w = 0.75/.pt
	font_size = 5.5
	gplot.bar <- ggplot(weight.mean.df, aes(x = seed, y = avg)) +
		geom_bar(stat = 'identity', position = position_dodge(), width = 0.6, fill = mycolor) + 
		geom_jitter(data = aln_df, aes(x = seed, y = weight), shape = 21, position = position_jitter(height = 0, width = .2), size = 0.8, stroke = 0.3) +
		labs(y = ylab) +
		scale_y_continuous(position = ypos, limits = c(ymin, ymax), breaks = ybreaks, expand = expansion(mult = 0)) +
		scale_x_discrete(position = xpos, expand = expansion(add = 0.6)) +
		coord_flip(expand = T, clip = 'off') +
		theme(
			legend.position='bottom',
			legend.text = element_text(size=8,color='black'),
			legend.title = element_blank(), 
			legend.background = element_blank(),
			legend.key.size = unit(0.25, "cm"),
			legend.margin=margin(0,0,0,0),
			legend.box.margin=margin(-8,0,0,0),
			legend.key = element_blank(),
			
			axis.title.y = element_blank(),
			axis.text.y = element_blank(), #element_text(size=8, color='black', family = 'Courier'),
			axis.ticks.y = element_blank(),
			axis.line.y = element_line(color='black',size=line_w),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.x = element_text(size=font_size,color='black'),
			axis.line.x = element_line(color='black',size=line_w),
			axis.ticks.x = element_line(color='black',size=line_w),
			
			plot.title = element_blank(), #element_text(color='black',size=10,hjust=0.5),
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), 'inch'),
			panel.background = element_blank(),
			panel.grid = element_blank()
		)
	
	# make a pie chart for each logo plot
	count.df = data.frame('counts' = sapply(aln_lst, function(x) {
		sum(x$kmer_count)
	}),
	'seed' = names(aln_lst), stringsAsFactors = F)
	gplot.pie.lst <- list()
	for (i in 1:nrow(count.df)) {
		pie.df <- data.frame('group' = c('sele', 'other'),
												 'counts' = c(count.df$counts[i], sum(count.df$counts)-count.df$counts[i]),
												 stringsAsFactors = F)
		pie.df$group <- factor(pie.df$group, levels = c('sele', 'other'))
		
		gplot.pie.lst[[i]] <- ggplot(pie.df, aes(x = '', y = counts, fill = group)) +
			geom_bar(stat = "identity") +
			coord_polar('y') +
			labs(y = paste0(round(count.df$counts[i],1),'/',sum(count.df$counts))) + 
			scale_fill_manual(values = c(mycolor, 'gray80')) + 
			theme(
				legend.position = 'none',
				axis.title.x = element_text(size=font_size-1.5,color='black', margin =  margin(0,0,0,0)),
				axis.title.y = element_blank(),
				axis.text.x = element_blank(),
				axis.text.y = element_blank(),
				panel.border = element_blank(),
				panel.grid = element_blank(),
				axis.ticks = element_blank(),
				plot.title = element_blank(),
				panel.background = element_blank(),
				plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'inch')
			)
	}
	
	# make a figure with all plots using the 'cowplot' package
	plot.all <- ggdraw() +
		draw_plot(gplot.bar, x = 0.43, y = 0, width = 0.57, height = 1) 
	for (i in 1:length(aln_lst)) {
		plot.all <- plot.all + 
			draw_plot(gplot.logo.lst[[i]], x = 0.08, y = (length(aln_lst)-i)/(length(aln_lst)+1), width = 0.35, height = 1/(length(aln_lst)+1)) +
			draw_plot(gplot.pie.lst[[i]], x = 0.01, y = (length(aln_lst)-i)/(length(aln_lst)+1), width = 0.07, height = 1/(length(aln_lst)+1)) 
	}
	
	if (is.null(f_width)){
		f_width <- 3
	}
	if (is.null(f_height)){
		f_height <- 0.75+0.2*length(aln_lst)
	}
	ggsave(file = paste0(fn, '_combined_plot.pdf'), plot = plot.all, width = f_width, height = f_height, units = "in")
}

kmer_tl_diff_by_ie_bplot <- function(df_in, y_col, mybreaks = NULL, pad = 0.1, fn = 'barplot'){
	# this function takes the kmer stat table output from iterative exclusion and makes a barplot for tail length difference
	df.plot <- df_in %>% mutate(round = factor(round, levels = rev(df_in %>% pull(round))))
	df.plot$y <- df_in %>% pull(y_col)
	kmer.len <- nchar(df.plot %>% pull(kmer) %>% {.[1]})
	
	# make a barplot
	if (df.plot$y[1] > 0) {
		my_palette = 'OrRd'
		palette_direction = 1
		xpos = 'bottom'
		ypos = 'right'
		ymin = 0
		if (is.null(mybreaks)){
			ymax = 1.05 * max(df.plot$y) 
			mybreaks = custom_breaks(ymin, ymax, length.out=6)
		}else{
			ymax = max(1.05 * max(df.plot$y), mybreaks)
		}
		label_pos <- ymax
		label_hjust <- 0
		plot_margin <- unit(c(0.05, 0.2+0.05*kmer.len, 0.05, 0.05), 'inch')
	} else {
		my_palette = 'GnBu'
		palette_direction = (-1)
		xpos = 'top'
		ypos = 'right'
		ymax = 0
		if (is.null(mybreaks)){
			ymin = 1.05 * min(df.plot$y)
			mybreaks = custom_breaks(ymin, ymax, length.out=6)
		}else{
			ymin = min(1.05 * min(df.plot$y), mybreaks)
		}
		label_pos <- ymin
		label_hjust <- 1
		plot_margin <- unit(c(0.05, 0.05, 0.05, 0.2+0.05*kmer.len), 'inch')
	}
	
	line_w = 0.75/.pt
	font_size = 6
	gplot <- ggplot(df.plot, aes(x = round, y = y, fill = y)) +
		geom_bar(stat = 'identity', width = 0.7) + 
		geom_errorbar(aes(ymin = y - sem, ymax = y + sem), width = 0.5, size = 0.1) +
		labs(x = 'Round (iterative exclusion)', y = 'Mean predicted tail-length difference (nt)') +
		scale_y_continuous(position = ypos, limits = c(ymin, ymax), breaks = mybreaks, labels = signs(mybreaks, accuracy = 0.1), expand = expansion(mult = 0)) +
		scale_fill_distiller(palette = my_palette, direction = palette_direction) +
		scale_x_discrete(position = xpos, expand = expansion(add = 0.6)) +
		geom_text(aes(x = round, y = label_pos, label = gsub('T', 'U', kmer)), size = font_size *5/14, hjust = label_hjust, family = 'Courier') +
		coord_flip(expand = T, clip = 'off') +
		theme(
			legend.position='none',
			
			axis.title.y = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.y = element_text(size=font_size,color='black'),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.x = element_text(size=font_size,color='black'),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(2, "pt"),
			#aspect.ratio = 1,		
			plot.background = element_rect(fill = 'transparent', color=NA),
			plot.margin = plot_margin,
			panel.background = element_blank(),
			panel.grid = element_blank()
		)
	#png(file=paste0(fn, '_', f.type, '_top_', top.cutoff, '_barplot.png'), width=3, height=1+0.15*nrow(tab.in),unit='in',res=300,bg='transparent')
	CairoPDF(file = paste0(fn, '.pdf'), width = 2, height = 2)
	print(gplot)
	dev.off()
}

kmer_single_nt_mutation_tl_diff_bplot <- function(df_in, col_x, col_y, col_y_sem = NULL, mybreaks = NULL, pad = 0.1, fn = 'barplot'){
	# this function takes a dataframe and makes barplot
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x), y = df_in %>% pull(col_y))
	
	if(!is.null(col_y_sem)){
		df.plot <- df.plot %>% mutate(sem = df_in %>% pull(col_y_sem))
	}
	
	# make a barplot
	if (df.plot$y[1] > 0) {
		my_palette = 'OrRd'
		palette_direction = 1
		xpos = 'bottom'
		ypos = 'right'
		ymin = 0
		if (is.null(mybreaks)){
			ymax = 1.05 * max(df.plot$y) 
			mybreaks = custom_breaks(ymin, ymax, length.out=6)
		}else{
			ymax = max(1.05 * max(df.plot$y), mybreaks)
		}
		label_pos <- ymax
		label_hjust <- 0
		plot_margin <- unit(c(0.05, 0.05, 0.05, 0.1), 'inch')
	} else {
		my_palette = 'GnBu'
		palette_direction = (-1)
		xpos = 'top'
		ypos = 'right'
		ymax = 0
		if (is.null(mybreaks)){
			ymin = 1.05 * min(df.plot$y)
			mybreaks = custom_breaks(ymin, ymax, length.out=6)
		}else{
			ymin = min(1.05 * min(df.plot$y), mybreaks)
		}
		label_pos <- ymin
		label_hjust <- 1
		plot_margin <- unit(c(0.05, 0.05, 0.05, 0.1), 'inch')
	}
	
	line_w = 0.75/.pt
	font_size = 6
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = y)) +
		geom_bar(stat = 'identity', width = 0.7) 
	
	if (!is.null(col_y_sem)){
		gplot <- gplot + 
			geom_errorbar(aes(ymin = y - sem, ymax = y + sem), width = 0.5, linewidth = 0.5/.pt) 
	}
	
	gplot <- gplot + 
		labs(y = 'Difference in tail-length change (nt)') +
		scale_y_continuous(position = ypos, limits = c(ymin, ymax), breaks = mybreaks, labels = signs(mybreaks, accuracy = 0.1), expand = expansion(mult = 0)) +
		scale_fill_distiller(palette = my_palette, direction = palette_direction) +
		scale_x_discrete(position = xpos, expand = expansion(add = 0.6)) +
		coord_flip(expand = T, clip = 'off') +
		theme(
			legend.position='none',
			
			axis.title.y = element_blank(),
			axis.text.y = element_text(size=font_size,color='black'),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.x = element_text(size=font_size,color='black'),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(2, "pt"),
			#aspect.ratio = 1,		
			plot.background = element_rect(fill = 'transparent', color=NA),
			plot.margin = plot_margin,
			panel.background = element_blank(),
			panel.grid = element_blank()
		)
	#png(file=paste0(fn, '_', f.type, '_top_', top.cutoff, '_barplot.png'), width=3, height=1+0.15*nrow(tab.in),unit='in',res=300,bg='transparent')
	CairoPDF(file = paste0(fn, '.pdf'), width = 2, height = 2)
	print(gplot)
	dev.off()
}

select_gene_ism_score_heatmap_w_logo <- function(mat_in, utr_label = NULL, pos_range = NULL, xbreak = 10, g_fsize = 5.5, block_size = 0.2, nt_label_scaler = 1.1, line_w = 0.75/.pt, fn = 'result', offset = 0){
	# This function takes in a tail-length change matrix, with x as positions and y as nucleotide identity (ACGU)
	# It subsets the position range given by 'pos_range' and makes 3 plots
	# 1) heatmap of tail-length change of all mutations at all positions
	# 2) line plot of max and min tail-length change for each position
	# 3) logo plot of inverse mean tail-length change of 3 mutations for each position (indicating importance)
	# note 'block_size' has a unit of 'inch'
	
	# subset data with regions of interest
	if (!is.null(pos_range)){
		mat_in = mat_in[, (pos_range[1]+ncol(mat_in)+1):(pos_range[2]+ncol(mat_in)+1)]
		if (!is.null(utr_label)){
			utr_label = utr_label[(pos_range[1]+length(utr_label)+1):(pos_range[2]+length(utr_label)+1)]
		}
	}
	max_min = max(max(mat_in), abs(min(mat_in)))
	
	# main heatmap
	df_mat <- as.data.frame(mat_in)
	
	# set wt values as NA
	#for (i in 1:4) {
	#	df_mat[i, colnames(df_mat) == c('A', 'C', 'G', 'U')[i]] <- NA 
	#}
	
	colnames(df_mat) <- 1:ncol(mat_in)
	df_mat <- df_mat %>% mutate(y = 4:1) %>% pivot_longer(cols = -y, names_to = 'x', values_to = 'val') %>% mutate(x = as.numeric(x))
	df_x_labels <- data.frame('seq' = colnames(mat_in)) %>% mutate('color' = case_when(seq == 'A' ~ brewer.pal(9, 'Set1')[1],
																																										 seq == 'C' ~ brewer.pal(9, 'Set1')[2],
																																										 seq == 'G' ~ brewer.pal(9, 'Set1')[3],
																																										 TRUE ~ brewer.pal(9, 'Set1')[5])) 
	df_na_labels <- data.frame('seq' = colnames(mat_in), 'x' = 1:ncol(mat_in), label = 'X') %>% 
		mutate(y = case_when(seq == 'A' ~ 4,
												 seq == 'C' ~ 3,
												 seq == 'G' ~ 2, 
												 TRUE ~ 1))
	
	htmap <- ggplot(df_mat, aes(x = x, y = y)) +
		geom_tile(aes(fill = val), color = 'white', linewidth = line_w/2) +
		coord_fixed(ratio = 1, clip = 'off', expand = F, xlim = c(0.5,ncol(mat_in)+0.5))+
		scale_x_continuous(breaks = seq(ncol(mat_in), 1, -xbreak),
											 labels = signs(seq(-1, -ncol(mat_in), -xbreak), accuracy = 1)) +
		scale_fill_gradient2(limits = c(-max_min, max_min), na.value = "gray90",
												 low = brewer.pal(9, 'Set1')[2], high = brewer.pal(9, 'Set1')[1], mid = 'white',
												 #low = brewer.pal(11,"RdBu")[11], high = brewer.pal(11,"RdBu")[1], mid = brewer.pal(11,"RdBu")[6],
												 breaks = custom_breaks(-max_min, max_min, length.out = 6, digits = 0),
												 labels = signs(custom_breaks(-max_min, max_min, length.out = 6, digits = 0), accuracy = 1),
												 guide = guide_colorbar(title.position = 'left', title = 'DPTLC (nt)')) +
		annotate('text', x = -0.3, y = 4:1, color = brewer.pal(9, 'Set1')[c(1:3,5)], 'label' = c('A', 'C', 'G', 'U'), size = block_size*25.4*nt_label_scaler) +
		annotate('text', x = seq(1,ncol(mat_in),1), y = -0.3, color = df_x_labels[,'color'], 'label' = df_x_labels[,'seq'], size = block_size*25.4*nt_label_scaler) +
		annotate('text', x = seq(ncol(mat_in), 1, -xbreak) - offset, y = -1.3, 'label' = signs(seq(-1, -ncol(mat_in), -xbreak), accuracy = 1), size = block_size*25.4*nt_label_scaler) +
		geom_text(data = df_na_labels, aes(x = x, y = y, label = label), color = 'gray50', size = block_size*25.4*nt_label_scaler*0.6) +
		theme(
			legend.position='right',
			legend.title.align = 0.5,
			legend.title = element_text(size = g_fsize, angle = 90, margin = margin(r = 2)),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(g_fsize, 'pt'),
			legend.text.align = 1,
			legend.text = element_text(size = g_fsize, margin = margin(l = 1)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(l = -8),
			
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.title = element_blank(),
			plot.margin = unit(c(0.0,0.01,0.0,0.01), "inch"),
			
			axis.title = element_blank(),
			axis.text = element_blank(),
			axis.line = element_blank(),
			axis.ticks = element_blank()
		)
	
	# line plot
	df_mm = bind_rows(data.frame('score' = apply(mat_in, 2, max), 'pos' = 1:ncol(mat_in) - 0.5, group = 'max'),
										data.frame('score' = apply(mat_in, 2, min), 'pos' = 1:ncol(mat_in) - 0.5, group = 'min'))
	line_plot <- ggplot(df_mm, aes(x = pos, y = score, group = group)) +
		geom_line(aes(color = group), size = line_w) +
		scale_color_manual(values = brewer.pal(9, 'Set1')[1:2], labels = c('Max', 'Min')) +
		coord_cartesian(clip = 'on', expand = F, xlim = c(0, ncol(mat_in)), ylim = c(-max_min*1.05, max_min*1.05)) +
		geom_hline(yintercept = 0, size = line_w, lty = 2) + 
		scale_y_continuous(breaks = custom_breaks(-max_min, max_min, length.out = 6, digits = 0),
											 labels = signs(custom_breaks(-max_min, max_min, length.out = 6, digits = 0), accuracy = 1)) +
		theme(
			legend.position='right',
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(g_fsize, 'pt'),
			legend.text = element_text(size = g_fsize, margin = margin(l = 1)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(l = -8),
			
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.title = element_blank(),
			plot.margin = unit(c(0.0,0.01,0.0,0.01), "inch"),
			
			axis.title.y = element_blank(),
			axis.text.y = element_text(size=g_fsize,color='black', hjust = 1, margin = margin(0,1,0,0)),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(block_size/5, "inch"),
			
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			axis.line.x = element_blank(),
			axis.ticks.x = element_blank()
		)
	
	# logo plot
	# weight is the average tail length change when a nucleotide is mutated to the other three
	logo_weight <- encode_seq_vec(colnames(mat_in)) %>% t() %>% {.*apply(mat_in,2,mean)*4/3} %>% t()
	rownames(logo_weight) <- c('A', 'C', 'G', 'U')
	logo_weight <- -logo_weight # inverse the values so that position values mean important
	
	cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'U'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
												cols=brewer.pal(9,'Set1')[c(1:3,5)]) # use custom color scheme for logo plots
	logoplot <- ggplot()+
		geom_logo(logo_weight, method = 'custom', seq_type = 'rna', col_scheme=cs1)+
		coord_cartesian(clip = 'on', expand = F) +
		theme(
			legend.position="none",
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			axis.line.x = element_blank(),
			axis.ticks.x = element_blank(),
			
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin = unit(c(0.0,0.01,0.0,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank()
		)
	
	CairoPDF(file=paste0('INN_ism_', fn,'_logo_lineplot_heatmap.pdf'), width=block_size*(ncol(mat_in) + 5), height=block_size*15)
	plots <- plot_grid(logoplot, line_plot, htmap, ncol = 1, rel_heights = c(1,1.5,2.5), align = 'v', axis = 'lr')
	print(plots)
	dev.off()
	invisible(df_mat)
}
