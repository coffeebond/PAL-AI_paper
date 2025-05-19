library('Cairo')
library('ggrepel')
library('RColorBrewer')
library('viridis', quietly = T)
library('ggpointdensity', quietly = T)
library('signs', quietly = T)
library('circlize')
library('cowplot')
library('ggseqlogo')

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

multi_sample_r_value_barplot <- function(df_in, x = 'gene', y = 'r_value', dodge = 'r_type', dodge_label = NULL, mycolor = NULL, ybreaks = NULL, xlab = NULL, ylab = NULL, ft = 'pdf', fn = 'plot'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(x), y = df_in %>% pull(y), z = df_in %>% pull(dodge))
	
	if (is.null(ybreaks)){
		ymin = min(df.plot$y)
		ymax = max(df.plot$y)
		ybreaks = custom_breaks(ymin, ymax, length.out = 7)
	} else{
		ymin = min(ybreaks)
		ymax = max(ybreaks)
	}
	
	if (is.null(dodge_label)){
		dodge_label = df.plot %>% pull(z) %>% unique()
	}
	
	if (is.null(mycolor)){
		mycolor = brewer.pal(8, 'Set2')[1:length(unique(df.plot$z))]
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = z)) +
		geom_bar(stat = 'identity', position = position_dodge(width = 0.55), width = 0.5) +
		coord_cartesian(ylim = c(ymin, ymax)) +
		labs(x = xlab, y = ylab) +
		scale_y_continuous(breaks = ybreaks, expand = expansion(mult = c(0, 0.1))) +
		scale_fill_manual(values = mycolor, labels = dodge_label) +
		theme(
			legend.position=c(1,1),
			legend.justification = c(1,1),
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
			axis.text.x = element_text(size=font_size,color='black', margin = margin(1,0,0,0), angle = 45, hjust = 1),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(2, "pt"),
			
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(font_size - font_size_offset, 'pt'),
			legend.text = element_text(size = font_size - font_size_offset, margin = margin(0,0,0,0)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0)
		)
	
	if (ft =='pdf'){
		CairoPDF(file=paste0(fn,'.pdf'), width=4, height=2)
	} else {
		png(file=paste0(fn,'.png'), width=4, height=2, unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
}

multi_sample_r_value_jitter_barplot <- function(df_in, x = 'gene', y = 'r_value', dodge = 'r_type', dodge_label = NULL, mycolor = NULL, ybreaks = NULL, xlab = NULL, ylab = NULL, ft = 'pdf', fn = 'plot'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(x), y = df_in %>% pull(y), z = df_in %>% pull(dodge))
	
	if (is.null(ybreaks)){
		ymin = min(df.plot$y)
		ymax = max(df.plot$y)
		ybreaks = custom_breaks(ymin, ymax, length.out = 7)
	} else{
		ymin = min(ybreaks)
		ymax = max(ybreaks)
	}
	
	if (is.null(dodge_label)){
		dodge_label = df.plot %>% pull(z) %>% unique()
	}
	
	if (is.null(mycolor)){
		mycolor = brewer.pal(8, 'Set2')[1:length(unique(df.plot$z))]
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = z)) +
		geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.55), width = 0.5) +
		geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), position = position_dodge(width = 0.55), width = 0.3, linewidth = 0.3) +
		geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.55, dodge.width = 0.55), size = 0.2, stroke = 0.2, color = 'black') +
		coord_cartesian(ylim = c(ymin, ymax)) +
		labs(x = xlab, y = ylab) +
		scale_y_continuous(breaks = ybreaks, expand = expansion(mult = c(0, 0.1))) +
		scale_fill_manual(values = mycolor, labels = dodge_label) +
		theme(
			legend.position=c(1,1),
			legend.justification = c(1,1),
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
			axis.text.x = element_text(size=font_size,color='black', margin = margin(1,0,0,0), angle = 45, hjust = 1),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(2, "pt"),
			
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(font_size - font_size_offset, 'pt'),
			legend.text = element_text(size = font_size - font_size_offset, margin = margin(0,0,0,0)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0)
		)
	
	if (ft =='pdf'){
		CairoPDF(file=paste0(fn,'.pdf'), width=4, height=2)
	} else {
		png(file=paste0(fn,'.png'), width=4, height=2, unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
}

multi_sample_r_value_jitter_boxplot <- function(df_in, x = 'gene', y = 'r_value', dodge = 'r_type', dodge_label = NULL, dodge_width = 0.6, mycolor = NULL, ybreaks = NULL, xlab = NULL, ylab = NULL, ft = 'pdf', fn = 'plot'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(x), y = df_in %>% pull(y), z = df_in %>% pull(dodge))
	
	if (is.null(ybreaks)){
		ymin = min(df.plot$y)
		ymax = max(df.plot$y)
		ybreaks = custom_breaks(ymin, ymax, length.out = 7)
	} else{
		ymin = min(ybreaks)
		ymax = max(ybreaks)
	}
	
	if (is.null(dodge_label)){
		dodge_label = df.plot %>% pull(z) %>% unique()
	}
	
	if (is.null(mycolor)){
		mycolor = brewer.pal(8, 'Set2')[1:length(unique(df.plot$z))]
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = z)) +
		#geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.55), width = 0.5) +
		#geom_violin(position = position_dodge(dodge_width), color = NA, width = dodge_width * 1.5, alpha = 0.75) +
		stat_summary(fun.data = custom_quantile, geom = 'boxplot', position = position_dodge(dodge_width), lwd = line_w, width = dodge_width * 0.9) +
		#geom_errorbar(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), position = position_dodge(width = 0.55), width = 0.3, linewidth = 0.3) +
		geom_point(pch = 21, position = position_jitterdodge(jitter.width = dodge_width * 0.8, dodge.width = dodge_width), size = 0.2, stroke = 0.2, color = 'black') +
		coord_cartesian(ylim = c(ymin, ymax)) +
		labs(x = xlab, y = ylab) +
		scale_y_continuous(breaks = ybreaks, expand = expansion(mult = c(0, 0.1))) +
		scale_fill_manual(values = mycolor, labels = dodge_label) +
		theme(
			legend.position=c(1,1),
			legend.justification = c(1,1),
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
			panel.grid.major.y = element_line(color = 'gray50', linewidth=0.25/.pt, linetype = 'dotted'),
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5, margin = margin(0,0,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin = margin(1,0,0,0), angle = 45, hjust = 1),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length = unit(2, "pt"),
			
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.key.size = unit(font_size - font_size_offset, 'pt'),
			legend.text = element_text(size = font_size - font_size_offset, margin = margin(0,0,0,0)), #, family = 'mono'
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0)
		)
	
	if (ft =='pdf'){
		CairoPDF(file=paste0(fn,'.pdf'), width=4, height=2)
	} else {
		png(file=paste0(fn,'.png'), width=4, height=2, unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
}

plot_scatter_w_side_decor <- function(df_in, col_x, col_y, col_label = NULL, col_fill = NULL, label_sele = NULL, color_sele = NULL, xbreaks = NULL, ybreaks = NULL, xy_equal = FALSE, legend_position = 'none', stat_legend_position = 'topright', labs = c('x', 'y'), pad = 0.05, fn = 'plot', ft = 'pdf') {
	###--- A scatter plot with density plots at top and right ---###
	df.plot <- df_in %>% mutate(x = .data[[col_x]], 
															y = .data[[col_y]])

	# labels 
	if (!is.null(label_sele)){
		# if "label_sele" is provided, the colors will used for selected points in "col_label"
		if (is.null(col_label)){
			stop('When "label_sele" is specified, "col_label" cannot be NULL!')
		}
		df.plot <- df.plot %>% mutate(label = as.character(.data[[col_label]])) %>%
			mutate(label = ifelse(label %in% label_sele, label, '')) %>%
			mutate(label = gsub('T', 'U', label)) %>% mutate(label = factor(label)) %>%
			arrange(desc(label))
	} 
	
	# point fill colors
	if (!is.null(col_fill)){
		df.plot <- df.plot %>% 
			mutate(fill_group = factor(.data[[col_fill]]))
		if (!is.null(color_sele)){
			mycolor <- color_sele
		} else {
			mycolor <- brewer.pal(9, 'Set1')[1:length(levels(df.plot$fill_group))]
		}
	} else {
		df.plot <- df.plot %>% 
			mutate(fill_group = 'others')
		mycolor <- alpha('gray50',0.3)
	}
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x)
		xmax <- max(df.plot$x)
		xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
	} else {
		xmin <- min(xbreaks, df.plot$x)
		xmax <- max(xbreaks, df.plot$x)
	}
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 7)
	} else {
		ymin <- min(ybreaks, df.plot$y)
		ymax <- max(ybreaks, df.plot$y)
	}
	
	if (xy_equal){ # in case x and y must be the same scale
		if (xmin <= ymin & xmax >= ymax){
			ymin <- xmin
			ymax <- xmax
			ybreaks <- xbreaks
		} else if (xmin >= ymin & xmax <= ymax){
			xmax <- ymax
			xmin <- ymin
			xbreaks <- ybreaks
		} else {
			xmin <- min(xmin, ymin)
			ymin <- xmin
			xmax <- max(xmax, ymax)
			ymax <- xmax
			xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
			ybreaks <- xbreaks
		}
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	rs <- round(cor(df.plot %>% select(x, y), method = 's')[1,2],2)
	rp <- round(cor(df.plot %>% select(x, y), method = 'p')[1,2],2)
	n <- nrow(df.plot)
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = fill_group)) + 
		geom_point(size = 2/.pt, shape = 21, stroke = 0.1, color = alpha('black', 0.5)) +
		scale_fill_manual(values = mycolor) 	
	
	if(xy_equal){
		gplot <- gplot + geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = line_w, col = alpha('black', 0.4)) 
	}
	
	if (grepl('bottom', stat_legend_position)){
		y_slp_1 <- ymin + 0.02*(ymax - ymin)
		y_slp_2 <- ymin + 0.1*(ymax - ymin)
		y_slp_3 <- ymin + 0.18*(ymax - ymin)
	} else {
		y_slp_1 <- ymax - 0.18*(ymax - ymin)
		y_slp_2 <- ymax - 0.1*(ymax - ymin)
		y_slp_3 <- ymax - 0.02*(ymax - ymin)
	}
	if (grepl('left', stat_legend_position)){
		x_slp <- xmin + 0.01*(xmax - xmin)
	} else {
		x_slp <- xmax - 0.01*(xmax - xmin)
	}
	
	gplot <- gplot +
		labs(x = labs[1], y = labs[2]) +
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_x_continuous(breaks = xbreaks, label = signs(xbreaks), expand = expansion(mult = pad)) +
		scale_y_continuous(breaks = ybreaks, label = signs(ybreaks), expand = expansion(mult = pad)) +
		annotate('text', x = x_slp, y = y_slp_1, label = deparse(bquote(n==.(signs(n, format = scales::comma)))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_2, label = deparse(bquote(italic(R)[s]==.(rs))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_3, label = deparse(bquote(italic(R)[p]==.(rp))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		theme(
			plot.margin = unit(c(0,0,0.01,0.01), "inch"),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			
			legend.position=legend_position,
			legend.justification = c(0,1),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key.size = unit(font_size,'points'),
			legend.key = element_blank(),
			legend.box.margin = margin(0,0,0,0),
			legend.text = element_text(size=font_size-font_size_offset, color='black', margin = margin(0,0,0,0)),
			legend.margin = margin(-6,0,0,0),
			
			axis.text.y = element_text(size=font_size,color='black', margin=margin(0,0.5,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin=margin(0.5,0,0,0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black', margin=margin(0.5,0,0,0)),
			axis.title.y = element_text(size=font_size,color='black', margin=margin(0,0.5,0,0)),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length=unit(1.5/.pt, "mm")
			
			#aspect.ratio=1
		)
	
	if(!is.null(label_sele)){
		gplot <- gplot + geom_text_repel(data = df.plot, aes(x = x, y = y, label = label), max.overlaps = Inf,
																		 color = 'black', size = (font_size - font_size_offset) * 5 / 14, seed = 57,
																		 min.segment.length = 0, segment.size = 0.5/.pt, box.padding = 0.8)
	}
	
	if (length(levels(df.plot$fill_group)) <= 1){
		xplot <- ggplot(df.plot, aes(x = x)) +
			geom_density(color = 'black', linewidth = line_w)
	} else {
		xplot <- ggplot(df.plot, aes(x = x, color = fill_group)) +
			geom_density(linewidth = line_w) +
			scale_color_manual(values = mycolor)
	}
	xplot <- xplot +
		coord_cartesian(xlim = c(xmin, xmax)) +
		scale_x_continuous(breaks = xbreaks, expand = expansion(mult = pad)) +
		theme(
			plot.margin = unit(c(0.01,0.01,0,0.01), "inch"),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			#line = element_blank(),
			text = element_blank(),
			axis.ticks.y = element_blank(),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			axis.line.y = element_blank(),
			axis.line.x = element_line(color='black',linewidth=line_w),
			#text = element_text(size=font_size,color='black'),
			title = element_blank(),
			legend.position = 'none'
		)

	if (length(levels(df.plot$fill_group)) <= 1){
		yplot <- ggplot(df.plot, aes(x = y)) +
			geom_density(color = 'black', linewidth = line_w)
	} else {
		yplot <- ggplot(df.plot, aes(x = y, color = fill_group)) +
			geom_density(linewidth = line_w) +
			scale_color_manual(values = mycolor)
	}
	yplot <- yplot +
		coord_flip(xlim = c(ymin, ymax)) +
		scale_x_continuous(breaks = ybreaks, expand = expansion(mult = pad)) +
		theme(
			plot.margin = unit(c(0.01,0.01,0.01,0), "inch"),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			#line = element_blank(),
			text = element_blank(),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			axis.line.x = element_blank(),
			axis.line.y = element_line(color='black',linewidth=line_w),
			#text = element_text(size=font_size,color='black'),
			title = element_blank(),
			legend.position = 'none'
		)
	
	# align plots separately
	aligned_xplot <- cowplot::align_plots(xplot, gplot, align = 'v', axis = 'l')[[1]]
	aligned_yplot <- cowplot::align_plots(yplot, gplot, align = 'h', axis = 'bl')[[1]]
	plots <- plot_grid(aligned_xplot, NULL, gplot, aligned_yplot, ncol = 2, rel_widths = c(6, 1), rel_heights = c(1, 6))
	
	design = "
		111111#
		2222223
		2222223
		2222223
		2222223
		2222223
		2222223
	"
	#plots <- xplot + gplot + yplot + plot_layout(design = design)
	#plots <- (((xplot | grid::textGrob('')) + plot_layout(ncol = 2, widths = c(1, 1))) / ((gplot | yplot) + plot_layout(ncol = 2, widths = c(1, 1)))) + plot_layout(heights = c(1, 1))
	if (ft == 'pdf'){
		CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	}
	print(plots)
	dev.off()	
	invisible(df.plot)
}

scatter_plot_with_groups <- function(df_in, x, y, group, xbreaks = NULL, ybreaks = NULL, xy_equal = FALSE, legend_position = 'none', stat_legend_position = 'topright', labs = c('x', 'y'), mycolor = NULL, fn ='plot'){
	df.plot <- df_in
	df.plot$x <- df_in %>% pull(x)
	df.plot$y <- df_in %>% pull(y)
	df.plot$group <- df_in %>% pull(group)
	n <- nrow(df.plot)
	rp <- cor(df.plot$x, df.plot$y, method = 'p') %>% round(., digits = 2)
	rs <- cor(df.plot$x, df.plot$y, method = 's') %>% round(., digits = 2)
	n_group <- length(unique(df.plot %>% pull(group)))
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x)
		xmax <- max(df.plot$x)
		xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
	} else {
		xmin <- min(xbreaks, df.plot$x)
		xmax <- max(xbreaks, df.plot$x)
	}
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 7)
	} else {
		ymin <- min(ybreaks, df.plot$y)
		ymax <- max(ybreaks, df.plot$y)
	}
	
	if (xy_equal){ # in case x and y must be the same scale
		if (xmin <= ymin & xmax >= ymax){
			ymin <- xmin
			ymax <- xmax
			ybreaks <- xbreaks
		} else if (xmin >= ymin & xmax <= ymax){
			xmax <- ymax
			xmin <- ymin
			xbreaks <- ybreaks
		} else {
			xmin <- min(xmin, ymin)
			ymin <- xmin
			xmax <- max(xmax, ymax)
			ymax <- xmax
			xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
			ybreaks <- xbreaks
		}
	}
	
	if (grepl('bottom', stat_legend_position)){
		y_slp_1 <- ymin + 0.02*(ymax - ymin)
		y_slp_2 <- ymin + 0.1*(ymax - ymin)
		y_slp_3 <- ymin + 0.18*(ymax - ymin)
	} else {
		y_slp_1 <- ymax - 0.18*(ymax - ymin)
		y_slp_2 <- ymax - 0.1*(ymax - ymin)
		y_slp_3 <- ymax - 0.02*(ymax - ymin)
	}
	if (grepl('left', stat_legend_position)){
		x_slp <- xmin + 0.01*(xmax - xmin)
	} else {
		x_slp <- xmax - 0.01*(xmax - xmin)
	}
	
	if (is.null(mycolor)){
		mycolor = brewer.pal(12, 'Set3')[1:n_group]
	}
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, group = group)) +
		geom_point(size = 1, aes(fill = group), shape = 21, color = t_col('black', 0.5), stroke = 0.2) + 
		labs(x = 'Measured poly(A) tail-length change (nt)', y = "Predicted poly(A) tail-length change (nt)") + 
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_x_continuous(breaks = xbreaks, labels = signs(xbreaks)) +
		scale_y_continuous(breaks = ybreaks, labels = signs(ybreaks)) +
		scale_fill_manual(values = mycolor) + 
		annotate('text', x = x_slp, y = y_slp_1, label = deparse(bquote(n==.(signs(n, format = scales::comma)))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_2, label = deparse(bquote(italic(R)[s]==.(rs))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_3, label = deparse(bquote(italic(R)[p]==.(rp))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		theme(
			legend.position='right',
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
			axis.text.x = element_text(size=font_size,color='black', margin = margin(1,0,0,0)),
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
	
	CairoPDF(paste0(fn, '.pdf'), width = 2.8, height = 2, family = 'Arial')
	print(gplot)
	dev.off()
}

plot_scatter <- function(df_in, col_x, col_y, col_label = NULL, col_fill = NULL, label_sele = NULL, color_sele = NULL, xbreaks = NULL, ybreaks = NULL, xy_equal = FALSE, legend_position = 'none', stat_legend_position = 'topright', labs = c('x', 'y'), pad = 0.05, fn = 'plot', ft = 'pdf') {
	###--- A scatter plot ---###
	df.plot <- df_in %>% mutate(x = .data[[col_x]], 
															y = .data[[col_y]])
	
	# labels 
	if (!is.null(label_sele)){
		# if "label_sele" is provided, the colors will be used for selected points in "col_label"
		if (is.null(col_label)){
			stop('When "label_sele" is specified, "col_label" cannot be NULL!')
		}
		df.plot <- df.plot %>% mutate(label = as.character(.data[[col_label]])) %>%
			mutate(label = ifelse(label %in% label_sele, label, '')) %>%
			mutate(label = gsub('T', 'U', label)) %>% mutate(label = factor(label)) %>%
			arrange(desc(label))
	} 
	
	# point fill colors
	if (!is.null(col_fill)){
		df.plot <- df.plot %>% 
			mutate(fill_group = factor(.data[[col_fill]]))
		if (!is.null(color_sele)){
			mycolor <- color_sele
		} else {
			mycolor <- brewer.pal(9, 'Set1')[1:length(levels(df.plot$fill_group))]
		}
	} else {
		df.plot <- df.plot %>% 
			mutate(fill_group = 'others')
		mycolor <- alpha('gray50',0.3)
	}
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x)
		xmax <- max(df.plot$x)
		xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
	} else {
		xmin <- min(xbreaks, df.plot$x)
		xmax <- max(xbreaks, df.plot$x)
	}
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 7)
	} else {
		ymin <- min(ybreaks, df.plot$y)
		ymax <- max(ybreaks, df.plot$y)
	}
	
	if (xy_equal){ # in case x and y must be the same scale
		if (xmin <= ymin & xmax >= ymax){
			ymin <- xmin
			ymax <- xmax
			ybreaks <- xbreaks
		} else if (xmin >= ymin & xmax <= ymax){
			xmax <- ymax
			xmin <- ymin
			xbreaks <- ybreaks
		} else {
			xmin <- min(xmin, ymin)
			ymin <- xmin
			xmax <- max(xmax, ymax)
			ymax <- xmax
			xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
			ybreaks <- xbreaks
		}
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	rs <- round(cor(df.plot %>% select(x, y), method = 's')[1,2],2)
	rp <- round(cor(df.plot %>% select(x, y), method = 'p')[1,2],2)
	n <- nrow(df.plot)
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = fill_group)) + 
		geom_point(size = 2/.pt, shape = 21, stroke = 0.1, color = alpha('black', 0.5)) +
		scale_fill_manual(values = mycolor) 	
	
	if(xy_equal){
		gplot <- gplot + geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = line_w, col = alpha('black', 0.4)) 
	}
	
	if (grepl('bottom', stat_legend_position)){
		y_slp_1 <- ymin + 0.02*(ymax - ymin)
		y_slp_2 <- ymin + 0.1*(ymax - ymin)
		y_slp_3 <- ymin + 0.18*(ymax - ymin)
	} else {
		y_slp_1 <- ymax - 0.18*(ymax - ymin)
		y_slp_2 <- ymax - 0.1*(ymax - ymin)
		y_slp_3 <- ymax - 0.02*(ymax - ymin)
	}
	if (grepl('left', stat_legend_position)){
		x_slp <- xmin + 0.01*(xmax - xmin)
	} else {
		x_slp <- xmax - 0.01*(xmax - xmin)
	}
	
	gplot <- gplot +
		labs(x = labs[1], y = labs[2]) +
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_x_continuous(breaks = xbreaks, label = signs(xbreaks), expand = expansion(mult = pad)) +
		scale_y_continuous(breaks = ybreaks, label = signs(ybreaks), expand = expansion(mult = pad)) +
		annotate('text', x = x_slp, y = y_slp_1, label = deparse(bquote(n==.(signs(n, format = scales::comma)))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_2, label = deparse(bquote(italic(R)[s]==.(rs))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_3, label = deparse(bquote(italic(R)[p]==.(rp))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		theme(
			plot.margin = unit(c(0,0,0.01,0.01), "inch"),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			
			legend.position=legend_position,
			legend.justification = c(0,1),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key.size = unit(font_size,'points'),
			legend.key = element_blank(),
			legend.box.margin = margin(0,0,0,0),
			legend.text = element_text(size=font_size-font_size_offset, color='black', margin = margin(0,0,0,0)),
			legend.margin = margin(-6,0,0,0),
			
			axis.text.y = element_text(size=font_size,color='black', margin=margin(0,0.5,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin=margin(0.5,0,0,0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black', margin=margin(0.5,0,0,0)),
			axis.title.y = element_text(size=font_size,color='black', margin=margin(0,0.5,0,0)),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length=unit(1.5/.pt, "mm")
			
			#aspect.ratio=1
		)
	
	if(!is.null(label_sele)){
		gplot <- gplot + geom_text_repel(data = df.plot, aes(x = x, y = y, label = label), max.overlaps = Inf,
																		 color = 'black', size = (font_size - font_size_offset) * 5 / 14, seed = 57,
																		 min.segment.length = 0, segment.size = 0.5/.pt, box.padding = 0.8)
	}
	
	if (ft == 'pdf'){
		CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()	
	invisible(df.plot)
}

plot_scatter_w_sem <- function(df_in, col_x, col_y, col_x_sem, col_y_sem, col_label = NULL, col_fill = NULL, label_sele = NULL, color_sele = NULL, xbreaks = NULL, ybreaks = NULL, xy_equal = FALSE, legend_position = 'none', stat_legend_position = 'topright', labs = c('x', 'y'), pad = 0.05, fn = 'plot', ft = 'pdf') {
	###--- A scatter plot with density plots at top and right ---###
	df.plot <- df_in %>% mutate(x = .data[[col_x]], 
															y = .data[[col_y]],
															x_sem = .data[[col_x_sem]],
															y_sem = .data[[col_y_sem]])
	# labels 
	if (!is.null(label_sele)){
		# if "label_sele" is provided, the colors will used for selected points in "col_label"
		if (is.null(col_label)){
			stop('When "label_sele" is specified, "col_label" cannot be NULL!')
		}
		df.plot <- df.plot %>% mutate(label = as.character(.data[[col_label]])) %>%
			mutate(label = ifelse(label %in% label_sele, label, '')) %>%
				mutate(label = gsub('T', 'U', label)) %>% mutate(label = factor(label)) %>%
					arrange(desc(label))
	} 

	# point fill colors
	if (!is.null(col_fill)){
		df.plot <- df.plot %>% 
			mutate(fill_group = factor(.data[[col_fill]]))
		if (!is.null(color_sele)){
			mycolor <- color_sele
		} else {
			mycolor <- brewer.pal(9, 'Set1')[1:length(levels(df.plot$fill_group))]
		}
	} else {
		df.plot <- df.plot %>% 
			mutate(fill_group = 'others')
		mycolor <- alpha('gray50',0.3)
	}

	# x and y axes
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x - df.plot$x_sem)
		xmax <- max(df.plot$x + df.plot$x_sem)
		xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
	} else {
		xmin <- min(xbreaks, df.plot$x - df.plot$x_sem)
		xmax <- max(xbreaks, df.plot$x + df.plot$x_sem)
	}
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y - df.plot$y_sem)
		ymax <- max(df.plot$y + df.plot$y_sem)
		ybreaks <- custom_breaks(ymin, ymax, length.out = 7)
	} else {
		ymin <- min(ybreaks, df.plot$y - df.plot$y_sem)
		ymax <- max(ybreaks, df.plot$y + df.plot$y_sem)
	}
	
	if (xy_equal){ # in case x and y must be the same scale
		if (xmin <= ymin & xmax >= ymax){
			ymin <- xmin
			ymax <- xmax
			ybreaks <- xbreaks
		} else if (xmin >= ymin & xmax <= ymax){
			xmax <- ymax
			xmin <- ymin
			xbreaks <- ybreaks
		} else {
			xmin <- min(xmin, ymin)
			ymin <- xmin
			xmax <- max(xmax, ymax)
			ymax <- xmax
			xbreaks <- custom_breaks(xmin, xmax, length.out = 7)
			ybreaks <- xbreaks
		}
	}
	
	# aesthetics
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	# correlation coefficient and count
	rs <- round(cor(df.plot %>% select(x, y), method = 's')[1,2],2)
	rp <- round(cor(df.plot %>% select(x, y), method = 'p')[1,2],2)
	n <- nrow(df.plot)
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, fill = fill_group)) + 
		#geom_point(size = 1/.pt, shape = 21, stroke = 0.1, color = alpha('black', 0.5)) +
		geom_pointrange(aes(xmin = x - x_sem, xmax = x + x_sem), linewidth = 0.5/.pt, size = 0.33, shape = 21, stroke = 0.2, color = alpha('black', 0.5)) + 
		geom_pointrange(aes(ymin = y - y_sem, ymax = y + y_sem), linewidth = 0.5/.pt, size = 0.33, shape = 21, stroke = 0.2, color = alpha('black', 0.5)) + 
		scale_fill_manual(values = mycolor) 	
	
	if(xy_equal){
		gplot <- gplot + geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = line_w, col = alpha('black', 0.4)) 
	}
	
	if (grepl('bottom', stat_legend_position)){
		y_slp_1 <- ymin + 0.02*(ymax - ymin)
		y_slp_2 <- ymin + 0.1*(ymax - ymin)
		y_slp_3 <- ymin + 0.18*(ymax - ymin)
	} else {
		y_slp_1 <- ymax - 0.18*(ymax - ymin)
		y_slp_2 <- ymax - 0.1*(ymax - ymin)
		y_slp_3 <- ymax - 0.02*(ymax - ymin)
	}
	if (grepl('left', stat_legend_position)){
		x_slp <- xmin + 0.01*(xmax - xmin)
	} else {
		x_slp <- xmax - 0.01*(xmax - xmin)
	}
	
	gplot <- gplot +
		labs(x = labs[1], y = labs[2]) +
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_x_continuous(breaks = xbreaks, label = signs(xbreaks), expand = expansion(mult = pad)) +
		scale_y_continuous(breaks = ybreaks, label = signs(ybreaks), expand = expansion(mult = pad)) +
		annotate('text', x = x_slp, y = y_slp_1, label = deparse(bquote(n==.(signs(n, format = scales::comma)))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_2, label = deparse(bquote(italic(R)[s]==.(rs))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = x_slp, y = y_slp_3, label = deparse(bquote(italic(R)[p]==.(rp))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		theme(
			plot.margin = unit(c(0,0,0.01,0.01), "inch"),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			
			legend.position=legend_position,
			legend.justification = c(0,1),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key.size = unit(font_size,'points'),
			legend.key = element_blank(),
			legend.box.margin = margin(0,0,0,0),
			legend.text = element_text(size=font_size-font_size_offset, color='black', margin = margin(0,0,0,0)),
			legend.margin = margin(-6,0,0,0),
			
			axis.text.y = element_text(size=font_size,color='black', margin=margin(0,0.5,0,0)),
			axis.text.x = element_text(size=font_size,color='black', margin=margin(0.5,0,0,0)),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black', margin=margin(0.5,0,0,0)),
			axis.title.y = element_text(size=font_size,color='black', margin=margin(0,0.5,0,0)),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length=unit(1.5/.pt, "mm")
			
			#aspect.ratio=1
		)
	
	if(!is.null(label_sele)){
		gplot <- gplot + geom_text_repel(data = df.plot, aes(x = x, y = y, label = label), max.overlaps = Inf,
																		 color = 'black', size = (font_size - font_size_offset) * 5 / 14, seed = 57,
																		 min.segment.length = 0, segment.size = 0.5/.pt, box.padding = 0.8)
	}
	
	if (ft == 'pdf'){
		CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()	
	invisible(df.plot)
}

select_gene_ism_score_heatmap_w_logo <- function(mat_in, utr_label = NULL, pos_range = NULL, xbreak = 10, g_fsize = 5.5, block_size = 0.2, nt_label_scaler = 1.1, line_w = 0.75/.pt, fn = 'result', offset = 0, max_min = NULL){
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
	
	if (is.null(max_min)){
		max_min = max(max(mat_in), abs(min(mat_in)))
	}
	
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
												 guide = guide_colorbar(title.position = 'left', title = 'DTLC (nt)')) +
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
	
	CairoPDF(file=paste0(fn,'_ISM_results_logo_lineplot_heatmap.pdf'), width=block_size*(ncol(mat_in) + 5), height=block_size*15)
	plots <- plot_grid(logoplot, line_plot, htmap, ncol = 1, rel_heights = c(1,1.5,2.5), align = 'v', axis = 'lr')
	print(plots)
	dev.off()
	invisible(plots)
}
