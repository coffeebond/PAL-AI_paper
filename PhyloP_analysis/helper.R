library('signs')
library('Cairo')
library('viridis')
library('RcppAlgos')
library('stats')
library('RColorBrewer')

# a function to adjust exon coordinates for each mRNA
adjust_exons <- function(exons_group, l_max = 2000) {
	mRNA_name <- exons_group$name[1]
	total_length <- sum(width(exons_group))
	
	if (total_length <= l_max) return(exons_group)
	
	# Sort exons correctly based on strand, from 3' to 5'
	if (as.character(strand(exons_group)[1]) == "+") {
		exons_group <- exons_group[order(start(exons_group), decreasing = TRUE)]
	} else {
		exons_group <- exons_group[order(end(exons_group))]
	}
	
	# Calculate cumulative lengths from 3' end
	cum_lengths <- cumsum(width(exons_group))
	
	# Find exons to keep (simplified logic)
	keep_idx <- which(cum_lengths <= l_max)
	if (length(keep_idx) == 0 || max(cum_lengths[keep_idx]) < l_max) {
		keep_idx <- c(keep_idx, max(keep_idx, 0) + 1)
	}
	keep_idx <- keep_idx[keep_idx <= length(exons_group)]  # Ensure within bounds
	
	# Get exons to keep
	trimmed_exons <- exons_group[keep_idx]
	current_length <- sum(width(trimmed_exons))
	
	# Trim last exon if needed
	if (current_length > l_max) {
		last_exon <- trimmed_exons[length(trimmed_exons)]
		trim_amount <- current_length - l_max
		
		if (as.character(strand(last_exon)) == "+") {
			start(last_exon) <- start(last_exon) + trim_amount
		} else {
			end(last_exon) <- end(last_exon) - trim_amount
		}
		
		trimmed_exons[length(trimmed_exons)] <- last_exon
	}
	
	# Restore original order (5' to 3')
	if (as.character(strand(exons_group)[1]) == "+") {
		trimmed_exons <- trimmed_exons[order(start(trimmed_exons))]
	} else {
		trimmed_exons <- trimmed_exons[order(end(trimmed_exons), decreasing = TRUE)]
	}
	
	return(trimmed_exons)
}

# a function to sort exons from 3' to 5' and compute relative positions
compute_relative_positions <- function(exons_group) {
	# sort exons based on strand
	if (as.character(strand(exons_group)[1]) == "+") {
		exons_sorted <- exons_group[order(end(exons_group), decreasing = TRUE)]
	} else {
		exons_sorted <- exons_group[order(start(exons_group), decreasing = FALSE)]
	}
	
	# compute cumulative lengths from the 3' end
	exon_lengths <- width(exons_sorted)
	cumulative_lengths <- cumsum(c(0, exon_lengths[-length(exon_lengths)]))
	
	# create a per-base GRanges object with relative positions
	pos_bed <- unlist(tile(exons_sorted, width = 1))
	pos_offset <- rep(cumulative_lengths, exon_lengths)
	
	# calculate relative positions
	if (as.character(strand(exons_sorted)[1]) == "+") {
		exon_ends <- rep(end(exons_sorted), exon_lengths)
		pos_bed$dis2end <- exon_ends - start(pos_bed) + pos_offset
	} else {
		exon_starts <- rep(start(exons_sorted), exon_lengths)
		pos_bed$dis2end <- start(pos_bed) - exon_starts + pos_offset
	}
	
	if ('pa_id' %in% colnames(mcols(exons_group))){
		pos_bed$pa_id <- exons_group$pa_id[1]
	}
	if ('name' %in% colnames(mcols(exons_group))){
		pos_bed$name <- exons_group$name[1]
	}
	return(pos_bed)
}

# a function to extract and tile positions of a motif in a DNAStringSet object
tile_motif_positions <- function(seqs, motif, fixed = TRUE){
	motif_pos_lst <- vmatchPattern(motif, seqs, fixed = fixed)
	pbmclapply(1:length(seqs), function(i){
		if (length(motif_pos_lst[[i]]) == 0){
			return(list())
		}
		list(
			'pa_id' = names(seqs)[i],
			'pos' = start(unlist(tile(motif_pos_lst[[i]], width = 1))),
			'flag' = TRUE
		)
	}, mc.cores = 40) %>% bind_rows 
} 

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

phylop_box_plot <- function(df_in, group_col = 'group', y_col = 'diff_tl', group_label = NULL, color_sele = NULL, ybreaks = NULL, fn = 'plot', ft = 'pdf'){
	###--- A box plot for phylop scores of select groups ---###
	df.plot <- df_in %>% mutate(y = df_in %>% pull(y_col),
															group = df_in %>% pull(group_col)) 
	n_group <- df.plot %>% group_by(group) %>% summarize(n = n()) %>% pull(n)
	
	# statistical tests
	group_combo <- comboGeneral(levels(df.plot$group), m = 2) 
	pvals <- pbmclapply(1:nrow(group_combo), function(x){
		x_series <- df.plot %>% filter(group == group_combo[x, 1]) %>% pull(y)
		y_series <- df.plot %>% filter(group == group_combo[x, 2]) %>% pull(y)
		bind_rows(
			list('group_x' = group_combo[x, 1], 
					 'group_y' = group_combo[x, 2], 
					 'n_x' = length(x_series),
					 'n_y' = length(y_series),
					 'pval' = wilcox.test(x = x_series, y = y_series, alternative = 'less')$p.value
			),
			list('group_x' = group_combo[x, 2], 
					 'group_y' = group_combo[x, 1], 
					 'n_x' = length(y_series),
					 'n_y' = length(x_series),
					 'pval' = wilcox.test(x = y_series, y = x_series, alternative = 'less')$p.value
			)
		) %>% return(.)
	}, mc.cores = 40) %>% bind_rows()
	
	# get the labels
	if (is.null(group_label)){
		mylabels <- paste0(levels(df.plot$group), ' (', signs(n_group, format = scales::comma), ')')
	} else {
		mylabels <- paste0(group_label, ' (', signs(n_group, format = scales::comma), ')')
	}
	
	# get the y breaks
	if (is.null(ybreaks)){
		ymin <- df.plot %>% pull(y) %>% quantile(0.05)
		ymax <- df.plot %>% pull(y) %>% quantile(0.95) 
		ybreaks <- custom_breaks(ymin, ymax, length.out = 8, digits = 1)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}
	
	# get the colors
	if (is.null(color_sele)){
		color_sele = brewer.pal(9, 'GnBu')[floor(seq(9,1, length.out = length(n_group)))]
	}
	
	line_w = 0.75/.pt
	font_size = 6
	gplot <- ggplot(df.plot, aes(x=group, y=y, fill = group)) +
		stat_summary(geom = 'boxplot', fun.data = custom_quantile, width = 0.4, linewidth = line_w) + 
		labs(y = "phyloP score") +
		coord_cartesian(ylim = c(ymin, ymax), clip = 'on') +
		scale_y_continuous(breaks = ybreaks, labels = ybreaks) +
		scale_fill_manual(values = color_sele, labels = mylabels)+
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_blank(),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black', angle = 45, hjust = 1),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			
			legend.title = element_blank(),
			legend.position='none',
			legend.justification=c(1,1),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size-0.5,color='black', margin = margin(l = 2,b = 0, unit = "points")),
			legend.key.height = unit(font_size, 'points'),
			legend.key.width = unit(font_size, 'points'),
			legend.key = element_blank(),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0),
			
			aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
	invisible(pvals)
}

pval_heatmap <- function(df_in, col_x, col_y, col_val, fn = 'plot', ft = 'pdf'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y),
															val = df_in %>% pull(col_val))
	df.plot <- df.plot %>% mutate(log_val = -log10(pmax(val,.Machine$double.xmin))) %>% 
		mutate(val_bin = cut(.$log_val,
												 breaks = c(0, 1, 2, 3, Inf),
												 labels = c("0-1", "1-2", "2-3", ">3"),
												 right = FALSE))  # Left-closed bins [a, b))
	
	# Create heatmap
	line_w = 0.75/.pt
	font_size = 5
	font_size_offset = 0.5
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = val_bin)) +
		geom_tile(color = 'white', linewidth = line_w/3) +
		scale_fill_manual(values = brewer.pal(11, 'BrBG')[c(2,7,9,11)]) +
		labs(fill = "-log10(p-value)") + 
		scale_x_discrete(position = 'top', guide = guide_axis(angle = 90)) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_text(size=font_size, face = "bold",color='black', hjust = 0.5, margin=margin(0,0,3,0)),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank(),
			
			axis.text.y = element_text(size=font_size, color='black', hjust = 0, margin=margin(t=0,r=0.5,b=0,l=0)),
			axis.text.x = element_text(size=font_size, color='black', hjust = 0, vjust = 0, margin=margin(t=0,r=0,b=0,l=0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			
			legend.position = 'right',
			#legend.justification = c(1,0),
			legend.title = element_text(size=font_size, color='black'),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size - font_size_offset,color='black',margin=margin(0,0,0,0)),
			legend.key = element_blank(),
			legend.key.size = unit(font_size, 'points'),
			legend.spacing.y = unit(font_size/2, 'points'),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0),
			
			aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn,'.pdf'), width = 3, height = 3, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=3, height=3,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
}