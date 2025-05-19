library('tidyverse')
library('ggplot2')
library('Cairo')
library('RcppAlgos')
library('RColorBrewer')
library('viridis')
library('signs', quietly = T)

custom.breaks <- function(bmin, bmax, digits = 0, length.out = 8, zero = TRUE) {
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

violin_plot <- function(df_in, col_x, col_y, label_x = NULL, ybreaks = NULL, col_fill = NULL, fn = 'plot', ft = 'pdf'){
	# make a violin plot
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y)) %>%
		mutate(x = factor(x))
	n_group <- df.plot %>% group_by(x) %>% summarise(n =n ()) %>%
		mutate(label = paste0('n = ', n))
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom.breaks(ymin, ymax, length.out = 6)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}
	
	if (is.null(label_x)){
		label_x <- 'x'
	}
	
	if (is.null(col_fill)){
		col_fill <- brewer.pal(9, 'Set1')
	}
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	gplot <- ggplot(data = df.plot, aes(x = x, y = y, fill = x)) + 
		geom_violin(color = NA, width = 0.9, alpha = 0.75) +
		stat_summary(fun.data = custom_quantile, geom = 'boxplot', lwd = line_w, width = 0.5, fill = NA) +
		coord_cartesian(ylim = c(ymin, ymax)) +
		labs(y = 'Pearson R', x = label_x) +
		scale_y_continuous(breaks = ybreaks) +
		scale_x_discrete(expand = expansion(add = c(0.75, 1))) +
		scale_fill_manual(values = col_fill) +
		geom_text(data = n_group, aes(x = x, label = label), y = ymax - 0.12*(ymax - ymin), size = font_size * 5/14, color = 'black', angle = 45, hjust = 0, vjust = 0) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_text(size=font_size, face = "bold",color='black', hjust = 0.5, margin=margin(0,0,3,0)),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank(),
			panel.grid.major.y = element_line(color = 'gray50', linewidth=0.25/.pt, linetype = 'dotted'),
			panel.grid.minor.y = element_blank(),
			
			axis.text.x = element_text(size=font_size,color='black', margin=margin(t=0.5,r=0,b=0,l=0)),
			axis.text.y = element_text(size=font_size,color='black', margin=margin(t=0,r=0.5,b=0,l=0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.y = element_text(size=font_size, color='black', margin=margin(0,1,0,0)),
			axis.title.x = element_text(size=font_size, color='black', margin=margin(1,0,0,0)),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			
			legend.position = 'none',
			legend.justification = c(0,1),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size - font_size_offset,color='black',margin=margin(0,0,0,-3)),
			legend.key = element_blank(),
			legend.key.size = unit(font_size, 'points'),
			legend.spacing.y = unit(font_size/2, 'points'),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0)
			
			#aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn,'.pdf'), width = 1, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=1, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
}

scatter_plot <- function(df_in, col_x, col_y, label_x = NULL, xbreaks = NULL, xlim = NULL, ybreaks = NULL, fn = 'plot', ft = 'pdf'){
	# make a violin plot
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y))
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom.breaks(ymin, ymax, length.out = 6)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x)
		xmax <- max(df.plot$x)
		xbreaks = 10 ** seq(floor(log10(xmin)), ceiling(log10(xmax)),1)
	} 
	xlabels = paste0('10^', log10(xbreaks))
	if (is.null(xlim)){
		xlim = c(min(xbreaks), max(xbreaks))
	}
	
	if (is.null(label_x)){
		label_x <- 'x'
	}

	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	gplot <- ggplot(data = df.plot, aes(x = x, y = y)) + 
		geom_point(size = 0.3, pch = 21, color = alpha('black', 0.2), fill = alpha('black', 0.1), stroke = 0.2) +
		coord_cartesian(ylim = c(ymin, ymax)) +
		labs(y = 'Pearson R', x = label_x) +
		scale_y_continuous(breaks = ybreaks) +
		scale_x_continuous(trans = 'log10', breaks = xbreaks, labels = do.call(expression, rlang::parse_exprs(xlabels)), expand = expansion(mult=0.2)) +
		annotate('text', x = xlim[1], y = ymax, hjust = 0, vjust = 1, label = paste0('n = ', signs(nrow(df.plot), format=scales::comma)), size = font_size * 5 /14) + 
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_text(size=font_size, face = "bold",color='black', hjust = 0.5, margin=margin(0,0,3,0)),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank(),
			panel.grid.major.y = element_line(color = 'gray50', linewidth=0.25/.pt, linetype = 'dotted'),
			panel.grid.minor.y = element_blank(),
			
			axis.text.x = element_text(size=font_size,color='black', margin=margin(t=0.5,r=0,b=0,l=0)),
			axis.text.y = element_text(size=font_size,color='black', margin=margin(t=0,r=0.5,b=0,l=0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.y = element_text(size=font_size, color='black', margin=margin(0,1,0,0)),
			axis.title.x = element_text(size=font_size, color='black', margin=margin(1,0,0,0)),
			axis.ticks.y = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			
			legend.position = 'none',
			legend.justification = c(0,1),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size - font_size_offset,color='black',margin=margin(0,0,0,-3)),
			legend.key = element_blank(),
			legend.key.size = unit(font_size, 'points'),
			legend.spacing.y = unit(font_size/2, 'points'),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0)
			
			#aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn,'.pdf'), width = 1.2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=1.2, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
}