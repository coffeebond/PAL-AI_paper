library('signs')
library('Cairo')
library('ggseqlogo')
library('ggrepel')
library('RColorBrewer')
library('Biostrings')

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

my_t_test <- function(m1, sd1, num1, m2, sd2, num2, diff = 0, lower.tail = TRUE) {
	# perform one-sided Welch's t-test against a set value 
	# null hypothesis: m1 >= m2 + diff
	se <- sqrt(sd1 ^ 2 / num1 + sd2 ^ 2 / num2)
	df <- ((sd1 ^ 2 / num1 + sd2 ^ 2 / num2) ^ 2 )  /((sd1 ^ 2 / num1) ^ 2 / (num1 - 1) + (sd2 ^ 2 / num2) ^ 2 / (num2 - 1) )
	t <- (m1 - m2 - diff) / se
	pval <- pt(t, df, lower.tail = lower.tail)
	#return(list('t' = t, 'pval' = pval))
	return(pval)
}

position_logo_plot <- function(df_in, ybreaks = NULL, xbreaks = NULL, xlabels = NULL, fn = 'plot'){
	df.plot <- df_in
	if (is.null(ybreaks)) {
		ymin <- min(df.plot)
		ymax <- max(df.plot)
		ybreaks <- custom_breaks(ymin, ymax, digits = 1, length.out = 8)
	} else {
		ymin <- min(ybreaks, df.plot)
		ymax <- max(ybreaks, df.plot)
	}
	
	if (is.null(xbreaks)) {
		xmin <- 1
		xmax <- ncol(df.plot)
		xbreaks <- custom_breaks(ymin, ymax, digits = 0, length.out = 8)
	} else {
		xmin <- min(xbreaks, 1)
		ymax <- max(xbreaks, ncol(df.plot))
	}
	
	if (is.null(xlabels) | (length(xlabels) != length(xbreaks))){
		xlabels <- xbreaks
	}
	
	cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'U'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
												cols=brewer.pal(9,'Set1')[c(1:3,5)])
	line_w = 0.75/.pt
	font_size = 5.5
	gplot <- ggplot()+
		geom_logo(df.plot, method = 'custom', seq_type = 'rna', col_scheme=cs1)+
		scale_x_continuous(breaks = xbreaks, labels = xlabels) +
		scale_y_continuous(breaks = ybreaks, labels = signs(ybreaks, accuracy = 0.1)) +
		labs(y = 'Relative tail-lenght difference (nt)', x = 'Distance to CPE (nt)') +
		geom_hline(yintercept = 0, linewidth = line_w) +
		theme(
			legend.position="none",
			axis.title.y = element_text(size=font_size,color='black',vjust=1),
			axis.text.y = element_text(size=font_size,color='black'),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=1),
			axis.text.x = element_text(size=font_size,color='black',vjust=1),
			axis.line.x = element_blank(),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank()
		)
	CairoPDF(file=paste0(fn,'.pdf'), width=0.1*ncol(df.plot), height=3, family = 'Arial')
	print(gplot)
	dev.off()
}

tl_change_vs_count_splot <- function(ctrl, expt, kmer_len = 5, count_cutoff = 1, kmer_sele_lst = NULL, kmer_labels = NULL, plot_label = F, color_sele = NULL, color_labels = NULL, ybreaks = NULL, xlim = NULL, p_cutoff = 0.01, xpad = 0.05, show_sig = FALSE, fn = 'plot', ft = 'pdf'){
	# This function takes in position-less datasets and calculate difference in mean pA length for each kmer (selected length)
	# and it makes a scatter plot comparing the difference in pA length to mean count number
	
	# get the difference of the mean for all variants
	mean_diff <- expt %>% filter(count == max(count)) %>% dplyr::slice(1) %>% pull(tl_avg) - ctrl %>% filter(count == max(count)) %>% dplyr::slice(1) %>% pull(tl_avg)
	
	# calculate difference of the mean for all kmers and perform Welch's t-test against the difference of the mean for all variants
	expt <- expt %>% filter(nchar(kmer) == kmer_len & count >= count_cutoff) %>% mutate(tl_std = sqrt(tl_sd / (count-1))) %>% select(kmer, tl_avg, tl_std, count)
	ctrl <- ctrl %>% filter(nchar(kmer) == kmer_len & count >= count_cutoff) %>% mutate(tl_std = sqrt(tl_sd / (count-1))) %>% select(kmer, tl_avg, tl_std, count)
	res.df <- inner_join(expt, ctrl, by = 'kmer', suffix = c('_e', '_c')) %>%
		mutate(diff = tl_avg_e - tl_avg_c, count = (count_e + count_c)/2, diff_sem = (tl_std_e ** 2 / count_e + tl_std_c ** 2 / count_c) ** 0.5 ) %>%
		rowwise() %>% mutate(p_increase = my_t_test(tl_avg_e, tl_std_e, count_e, tl_avg_c, tl_std_c, count_c, diff = mean_diff, lower.tail = F),
												 p_decrease = my_t_test(tl_avg_e, tl_std_e, count_e, tl_avg_c, tl_std_c, count_c, diff = mean_diff, lower.tail = T)) %>% ungroup() %>%
		mutate(log10p = ifelse(diff > mean_diff, -log10(pmax(p_increase, .Machine$double.xmin)), -log10(pmax(p_decrease, .Machine$double.xmin))))
	diff_avg <- mean(res.df$diff)
	diff_sd <- sd(res.df$diff)
	res.df$z_score <- (res.df$diff - diff_avg) / diff_sd
	res.df <- res.df %>% arrange(desc(z_score)) %>% mutate(rank = 1:nrow(.))
	
	# get the adjust p value cutoff for multiple hypothesis testing
	padj_cutoff <- -log10(p_cutoff / nrow(res.df))
	res.df <- res.df %>% mutate(sig = ifelse(log10p > padj_cutoff, TRUE, FALSE)) %>% mutate(sig = factor(sig, levels = c(TRUE, FALSE)))
	
	if (!is.null(kmer_sele_lst)) {
		if (!is.list(kmer_sele_lst)) stop('"kmer_sele_lst" must be a list!')
		
		# get the colors
		if (is.null(color_sele)){
			mycolor = c(alpha('black', 0.15), brewer.pal(9, 'Set1')[1:length(kmer_sele_lst)])
		} else {
			mycolor = c(alpha('black', 0.15), color_sele)
		}
		names(mycolor) <- mycolor # must be named vector, otherwise it won't work in ggplot2 (scale_fill/color_manual)
		res.df <- res.df %>% mutate(color = mycolor[1])
		if (is.null(color_labels)) color_labels = mycolor[2:length(mycolor)]
		
		# get the label if it's not provided
		if (is.null(kmer_labels)){
			kmer_labels = unlist(kmer_sele_lst)
		}
		res.df <- res.df %>% mutate(label = ifelse(kmer %in% kmer_labels, gsub('T', 'U', kmer), ''))
		
		# add color and group
		for (i in 1:length(kmer_sele_lst)){
			kmer_count <- sapply(kmer_sele_lst[[i]], function(x){vcountPattern(x, DNAStringSet(res.df$kmer), fixed = F)}) %>% apply(., 1, sum)
			res.df <- res.df %>% 
				mutate(color = ifelse(kmer_count != 0, mycolor[i+1], color))
		}
		
		res.df <- res.df %>% mutate(color = factor(color, levels = unique(mycolor))) %>% arrange(color)
	} else {
		plot_label = F
	}
	
	
	# x-axis
	xmin = min(res.df$count)
	xmax = max(res.df$count)
	xbreaks = 10 ** seq(floor(log10(xmin)), ceiling(log10(xmax)),1)
	xlabels = paste0('10^', log10(xbreaks))
	if (is.null(xlim)){
		xlim = c(min(xbreaks), max(xbreaks))
	}
	
	# y-axis
	if(is.null(ybreaks)){
		ymin = min(res.df$diff)
		ymax = max(res.df$diff)
		ybreaks = custom_breaks(ymin, ymax, digits = 1, length.out = 7)
	} else {
		ymin = min(ybreaks, min(res.df$diff))
		ymax = max(ybreaks, max(res.df$diff))
	}
	
	z_min = min(min(res.df$z_score), min((ybreaks - diff_avg)/diff_sd))
	z_max = max(max(res.df$z_score), max((ybreaks - diff_avg)/diff_sd))
	z_label = unique(c(rev(seq(0, z_min, by=-3)), 0, seq(0, z_max, by=3)))
	n_kmer = nrow(res.df)
	
	# only draw legend if no labels for kmers are drawn on the plot
	if (plot_label | is.null(kmer_sele_lst)){ 
		legend_pos = 'none'
	} else {
		legend_pos = c(0,1)
	}
	
	# show differnt shapes for significant kmers
	if (show_sig){
		myshape <- c(21, 22)
	} else myshape <- c(21, 21)
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = ifelse(length(kmer_labels) < 6, 0.5, 1)
	
	gplot <- ggplot(res.df, aes(x = count, y = diff)) +
		geom_point(aes(fill = color, shape = sig, size = sig), color = alpha('black', 0.3), stroke = 0.15) + 
		labs(x = 'Mean number of variants', y = "Difference in mean poly(A)-tail length (nt)") + 
		coord_cartesian(xlim = xlim, ylim = c(ymin, ymax)) +
		scale_shape_manual(values = myshape, guide = 'none') +
		scale_size_manual(values = c(1, 1), guide = 'none') +
		scale_x_continuous(trans = 'log10', breaks = xbreaks, labels = do.call(expression, rlang::parse_exprs(xlabels)), expand = expansion(mult=xpad)) +
		scale_y_continuous(breaks = ybreaks, labels = signs(ybreaks, accuracy = .1), expand = expansion(mult=0.1),
											 sec.axis = sec_axis(~ ((. - diff_avg) / diff_sd), name = 'Z-score', 
											 										breaks = z_label, labels = signs(z_label, accuracy = 1))) +
		scale_fill_manual(values = mycolor, breaks = mycolor[2:length(mycolor)], labels = color_labels) + 
		geom_hline(yintercept = diff_avg, linetype='dashed', linewidth = line_w, color = t_col('black', 0.5)) +
		annotate('text', x = xlim[1], y = ymin, hjust = 0, vjust = 1, label = paste0('n = ', signs(n_kmer, format=scales::comma)), size = (font_size - font_size_offset) * 5 /14) + 
		guides(fill = guide_legend(override.aes = list(fill = mycolor[2:length(mycolor)], shape = 21, size = 0.7))) +
		theme(
			legend.position=legend_pos,
			legend.justification = c(0,1),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.title = element_blank(),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			
			axis.title.y = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.y = element_text(size=font_size,color='black', hjust = 1, margin = margin(0,1,0,0)),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			axis.text.y.right = element_text(hjust = 1, margin = margin(0,0,0,1)),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5, margin = margin(0,0,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin = margin(0,0,0,0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
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
	if(plot_label){
		gplot <- gplot + geom_text_repel(data = res.df, aes(x = count, y = diff, label = label), max.overlaps = Inf,
																		 color = 'black', size = (font_size - font_size_offset) * 5 / 14, seed = 6,
																		 min.segment.length = 0, segment.size = 0.5/.pt, box.padding = 0.8)
	}
	if (ft == 'pdf') {
		CairoPDF(paste0(fn, '.pdf'), width = 2, height = 1.6, family = 'Arial')
	} else {
		png(file=paste0(fn, '.png'), width=2, height=1.6,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
	invisible(res.df)
}
