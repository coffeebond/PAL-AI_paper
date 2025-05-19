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
		scale_y_continuous(position = ypos, limits = c(ymin, ymax), breaks = mybreaks, labels = signs(mybreaks, accuracy = 1), expand = expansion(mult = 0)) +
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

motif_insertion_line_plot <- function(df_in, motifs = NULL, ybreaks = NULL, xbreak = 200, pos_lims = c(0, 1), col_sele = NULL, fn = 'my'){
	if (!is.null(motifs)){
		df_in <- df_in %>% filter(motif %in% motifs)
	}
	
	df_in <- df_in %>% mutate(label = gsub('T', 'U', motif)) %>%
		mutate(x = max(.$dis2end) - dis2end)
	
	if (is.null(ybreaks)){
		ymin <- min(df_in$avg - df_in$sem, na.rm = T)
		ymax <- max(df_in$avg + df_in$sem, na.rm = T)
		ybreaks = custom_breaks(ymin, ymax, length.out = 8)
	} else {
		ymin <- min(ybreaks, min(df_in$avg - df_in$sem, na.rm = T))
		ymax <- max(ybreaks, max(df_in$avg + df_in$sem, na.rm = T))
	}
	
	xbreaks <- seq(max(df_in$dis2end), min(df_in$dis2end), -xbreak)
	xlabels <- max(df_in$dis2end) - xbreaks
	
	xlimits <- pos_lims * max(df_in$x)
	
	if (is.null(col_sele)){
		col_sele = brewer.pal(9, 'Set1')[1:length(levels(df_in$motif))]
	}
	
	font_size <- 5.5
	font_size_offset = 0.5
	line_w <- 0.75/.pt
	
	gplot <- ggplot(df_in, aes(x = x, y = avg, color = motif, label = label)) +
		geom_line(size = 0.75/.pt) +
		geom_ribbon(aes(ymin=avg-sem, ymax=avg+sem, fill = motif), alpha = 0.3, color = NA) +
		scale_color_manual(name = NULL, values = col_sele) +
		scale_fill_manual(name = NULL, values = col_sele) +
		scale_y_continuous(name = 'Mean predicted tail-length change (nt)', limits = c(ymin, ymax), 
											 breaks = ybreaks, labels = signs(ybreaks), expand = expansion(mult = c(0.05, 0.05))) +
		scale_x_continuous(name = paste0('Distance to PAS (nt)'), 
											 breaks = xbreaks, 
											 limits = xlimits,
											 labels = xlabels, 
											 expand = expansion(mult = c(0.05,0.05))) +
		guides(color = guide_legend(byrow = TRUE)) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black', angle = 45, hjust = 1),
			axis.line.x = element_line(color='black', linewidth=line_w),
			axis.line.y = element_line(color='black', linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black'),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.2),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length = unit(2, 'pt'),
			
			legend.position=c(0,1),
			legend.justification = c(0,1),
			legend.text = element_text(size=font_size,color='black', margin = margin(0,0,0,0)),
			legend.title = element_text(size=font_size,color='black'),
			legend.background = element_blank(),
			#legend.key.size = unit(2,'pt'),
			legend.key.height = unit(font_size,'pt'),
			legend.key.width = unit(font_size,'pt'),
			legend.key = element_blank(),	
			legend.spacing.y = unit(2, 'pt'),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(1,0,0,3),
			
			aspect.ratio = 1
		) 
	#png(paste0(fn, '.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	CairoPDF(paste0(fn, '_line_plot.pdf'), width = 2, height = 2)
	print(gplot)
	dev.off()	
	invisible(df_in)
}
