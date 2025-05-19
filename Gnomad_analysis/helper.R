library('signs')
library('Cairo')
library('viridis')
library('RColorBrewer')

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
tl_diff_vs_position_scatter_plot <- function(df_in, x_col = 'pos', y_col = 'diff_tl', fill_col = 'af', xinterval = NULL, ybreaks = NULL, pos_lims = c(-100, -1), fn = 'plot', ft = 'pdf'){
	df.plot <- df_in
	df.plot$x <- as.data.frame(df_in)[,x_col]
	df.plot$y <- as.data.frame(df_in)[,y_col]
	df.plot$f <- log10(as.data.frame(df_in)[,fill_col])
	
	# get the y breaks
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y, ybreaks)
		ymax <- max(df.plot$y, ybreaks)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 8, digits = 1)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}
	
	n_sub <- df.plot %>% filter(x >= pos_lims[1] & x <= pos_lims[-1]) %>% nrow(.)
	font_size <- 6
	line_w <- 0.75/.pt
	font_size_offset <- 0.5
	gplot <- ggplot(df.plot, aes(x = x, y = y, color = f)) +
		geom_point(size = 0.3, shape = 1, stroke = 0.1) + 
		labs(x = 'Relative position to 3\' ends (nt)', y = "Predicted tail-length difference\nrelative to refenence alleles (nt)") + 
		coord_cartesian(ylim = c(ymin, ymax)) +
		scale_y_continuous(limits = c(ymin, ymax), breaks = ybreaks, 
											 labels = signs(ybreaks, format = scales::comma), expand = expansion(mult = c(0.05, 0.05))) +
		scale_x_continuous(breaks = c(-1, rev(seq(-xinterval, pos_lims[1], -xinterval))), 
											 limits = pos_lims,
											 labels = signs(c(-1, rev(seq(-xinterval, pos_lims[1], -xinterval))), accuracy = 1, format = scales::comma), 
											 expand = expansion(mult = c(0.05,0.05))) +
		scale_color_viridis(name = 'Allele frequency',option = 'magma', breaks = seq(-5,0,1), 
												labels = do.call(expression, rlang::parse_exprs(paste0('10^', seq(-5,0,1))))) + 
		annotate(geom = 'text', x = pos_lims[1], y = ymin, hjust = 0, vjust = 0, size = 5/14*(font_size-font_size_offset),
						 label = paste0('n = ', signs(n_sub, format = scales::comma))) +
		theme(
			legend.position=c(0,1),
			legend.justification = c(0,1),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.title = element_text(size=font_size,color='black', hjust = 0.5, face = 'bold'),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			
			axis.title.y = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.y = element_text(size=font_size,color='black', hjust = 1, margin = margin(0,1,0,0)),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			axis.text.y.right = element_text(hjust = 1, margin = margin(0,0,0,1)),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5, margin = margin(0,0,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin = margin(1,0,0,0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(2, "pt"),
			
			legend.title = element_text(size=font_size- font_size_offset,color='black'),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(font_size - font_size_offset, 'pt'),
			legend.text = element_text(size = font_size - font_size_offset, margin = margin(0,0,0,0)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,3),
			aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn, '.pdf'), width = 2.2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2.2, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
	invisible(df.plot)
}

tl_diff_cdf_plot <- function(df_in, groups, group_col = 'group', val_col = 'diff_tl', group_label = NULL, color_sele = NULL, xbreaks = NULL, ybreaks = NULL, fn = 'plot', ft = 'pdf'){
	###--- A cdf plot for predicted tail length changes of select groups ---###
	df.plot <- df_in
	df.plot$group <- as.data.frame(df_in)[,group_col]
	df.plot$xval <- as.data.frame(df_in)[,val_col]
	
	n_group <- df.plot %>% group_by(group) %>% summarize(n = n()) %>% pull(n)
	
	# statistical test if there are only two groups
	pval <- wilcox.test(x = df.plot %>% filter(group == groups[1]) %>% pull(xval), 
											y = df.plot %>% filter(group == groups[2]) %>% pull(xval),
											alternative = 'less')$p.value
	if (pval < 0.01 && pval != 0) {
		#pv.label <- bquote(italic(P)<10^.(ceiling(log10(pval))))
		pv.label <- bquote(italic(P) == .(format(pval*10^ceiling(-log10(pval)), digits = 2))*'\u00D7'*10^.(-floor(-log10(pval))))
	} else {
		pv.label <- bquote(italic(P)==.(format(pval, digits = 2)))
	}
	
	# get the colors
	if (is.null(group_label)){
		mylabels <- paste0(groups, ' (', signs(n_group, format = scales::comma), ')')
	} else {
		mylabels <- paste0(group_label, ' (', signs(n_group, format = scales::comma), ')')
	}
	
	# get the x breaks
	if (is.null(xbreaks)){
		xmin <- df.plot %>% pull(xval) %>% quantile(0.05)
		xmax <- df.plot %>% pull(xval) %>% quantile(0.95) 
		xbreaks <- custom_breaks(xmin, xmax, length.out = 8, digits = 1)
	} else {
		xmin <- min(xbreaks)
		xmax <- max(xbreaks)
	}
	
	# get the y breaks
	if (is.null(ybreaks)){
		ymin <- 0
		ymax <- 1
		xbreaks <- seq(0,1,0.25)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}
	
	# get the colors
	if (is.null(color_sele)){
		color_sele = brewer.pal(9, 'Set1')[1:length(groups)]
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	gplot <- ggplot(df.plot, aes(x=xval, group=group, color = group)) +
		stat_ecdf(geom = "step", size = 1/.pt )+
		labs(x = 'Predicted tail-length difference (nt)', y = "Cumulative fraction") +
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax), clip = 'on') +
		scale_x_continuous(breaks = xbreaks, labels = signs(xbreaks, format = scales::comma)) +
		scale_y_continuous(breaks = ybreaks, labels = ybreaks) +
		scale_color_manual(values = color_sele, labels = mylabels)+
		annotate('text', x = xmin + 0.05 *(xmax - xmin), y = 0.8, hjust = 0, label = deparse(pv.label), parse = T, size = font_size * 5 /14) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_blank(),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black'), #angle=45,hjust = 1),
			axis.line.x = element_line(color='black',size=line_w),
			axis.line.y = element_line(color='black',size=line_w),
			axis.title.x = element_text(size=font_size,color='black', vjust = 0.5, margin=margin(0.5,0,0,0)),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0),
			axis.ticks.x = element_line(color = "black", size = line_w),
			axis.ticks.y = element_line(color = "black", size = line_w),
			
			legend.title = element_blank(),
			legend.position=c(0.05,1),
			legend.justification=c(0,1),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size-0.5,color='black', margin = margin(l = -2,b = 3, unit = "points")),
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
	invisible(df.plot)
	
}

fraction_barplot <- function(df_in, y_col = 'frac_nor', x_col = 'label', xlabel = NULL, ybreaks = NULL, ylabel = 'Relative frequency (\u03994L < -25 nt)',color_sele = NULL, fn = 'plot'){
	df.plot <- df_in
	df.plot$y <- as.data.frame(df_in)[,y_col]
	df.plot$x <- as.data.frame(df_in)[,x_col]
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 6)
	} else {
		ymin <- min(ybreaks, df.plot$y)
		ymax <- max(ybreaks, df.plot$y)
	}
	
	if (is.null(color_sele)){
		mycolor <- rev(brewer.pal(9, 'GnBu'))[1:nrow(df.plot)]
	} else {
		mycolor <- color_sele
	}
	
	font_size <- 5.5
	line_w <- 0.75/.pt
	font_size_offset <- 0.5
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = x)) +
		geom_bar(stat='identity', width = 0.6) +
		coord_cartesian(ylim = c(ymin, ymax)) +
		scale_fill_manual(values = mycolor) +
		scale_y_continuous(name = ylabel, breaks = ybreaks, expand = expansion(mult = c(0, 0.05)),
											 labels = ybreaks) +
		#guides(fill = guide_legend(nrow = 1)) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_text(size = font_size, hjust =0.5, margin=margin(t=0,r=0,b=1,l=0)),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black',angle = 45, hjust = 1, margin=margin(t=0,r=0,b=0,l=0)),
			axis.line.x = element_line(color='black',size=line_w),
			axis.line.y = element_line(color='black',size=line_w),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.2),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_line(color = "black", size = line_w),
			axis.ticks.length = unit(2, 'pt'),
			
			legend.position='none'
		) 
	#png(paste0(fn, '.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	CairoPDF(paste0(fn, '.pdf'), width = 1.8, height = 2)
	print(gplot)
	dev.off()	
	invisible(df.plot)
}

fraction_barplot_by_bin <- function(df_in, col_x, col_y, col_dodge, col_anno = NULL, ybreaks = NULL, x_title = 'x', y_title = 'y', dodge_color = NULL, dodge_width = 0.7, fn = 'plot'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y),
															dodge = df_in %>% pull(col_dodge))
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 6)
	} else {
		ymin <- min(ybreaks, df.plot$y)
		ymax <- max(ybreaks, df.plot$y)
	}
	
	if (is.null(dodge_color)){
		dodge_color <- rev(brewer.pal(9, 'GnBu'))[1:length(df.plot %>% pull(dodge) %>% unique)]
	} 
	
	# add p-value
	if (!is.null(col_anno)){
		df.anno <- df.plot %>% mutate(y = ymax * 1.05, anno = df.plot %>% pull(col_anno)) %>%
			mutate(x_adj = as.numeric(x) + (as.numeric(dodge) - (length(unique(dodge)) + 1) / 2) * (dodge_width / length(levels(dodge))))
	}
	
	font_size <- 5.5
	line_w <- 0.75/.pt
	font_size_offset <- 0.5
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = dodge)) +
		geom_bar(stat='identity', position = position_dodge(width = dodge_width), width = dodge_width-0.1)
	
	if (!is.null(col_anno)){
		gplot <- gplot + 
			#geom_tile(data = df.anno, aes(fill = anno, x = x_adj, y = y), width = dodge_width/length(unique(df.tile$dodge)), color = 'white', linewidth = line_w/2, height = 0.05)
			geom_point(data = df.anno, aes(fill = anno, x = x_adj, y = y), shape = 21, size = 1, color = 'transparent')
	}
	
	gplot <- gplot +
		coord_cartesian(ylim = c(ymin, ymax)) +
		scale_fill_manual(values = c(dodge_color, brewer.pal(11, 'BrBG')[c(2,7,9,11)])) +
		labs(x = x_title, y = y_title) +
		scale_y_continuous(breaks = ybreaks, expand = expansion(mult = c(0, 0.1)),
											 labels = ybreaks) +
		#guides(fill = guide_legend(nrow = 1)) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_text(size = font_size, hjust =0.5, margin=margin(t=0,r=0,b=1,l=0)),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black', margin=margin(t=0,r=0,b=0,l=0)),
			axis.line.x = element_line(color='black', linewidth = line_w),
			axis.line.y = element_line(color='black', linewidth = line_w),
			axis.title.x = element_text(size=font_size,color='black', margin=margin(t=1,r=0,b=0,l=0)),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.2),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length = unit(2, 'pt'),
			
			legend.title = element_blank(),
			legend.position='right',
			legend.justification=c(0,1),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size,color='black', margin = margin(l = 0,b = 0, unit = "points")),
			legend.key.height = unit(font_size, 'points'),
			legend.key.width = unit(font_size, 'points'),
			legend.key = element_blank(),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0)
			
		) 
	#png(paste0(fn, '.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	CairoPDF(paste0(fn, '.pdf'), width = 4, height = 2)
	print(gplot)
	dev.off()	
	invisible(df.plot)
}

