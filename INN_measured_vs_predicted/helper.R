library('Cairo')
library('ggrepel')
library('RColorBrewer')
library('viridis', quietly = T)
library('ggpointdensity', quietly = T)
library('signs', quietly = T)
library('circlize')

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

prediction_w_sele_scatter_plot <- function(df_in, gene_id = 'idx', xaxis = 'y', yaxis = 'y_pred', xy_breaks = NULL, label_sele = NULL, color_sele = NULL, fn = 'plot', ft = 'pdf'){
	df.plot <- df_in %>% mutate(x = !!as.symbol(xaxis), y = !!as.symbol(yaxis), gene_id = !!as.symbol(gene_id)) %>%
		mutate(gene_id = sapply(gene_id, function(z){unlist(str_split(z, '__'))[1]})) 
	
	if (!is.null(label_sele)){
		df.plot <- df.plot %>% 
			mutate(label = ifelse(gene_id %in% label_sele, gene_id, ''),
						 group = ifelse(gene_id %in% label_sele, gene_id, 'others')) %>%
			mutate(group = factor(group, levels = c(label_sele, 'others'))) %>%
			arrange(desc(group))
		if (is.null(color_sele)){
			mycolor <- c(brewer.pal(9, 'Set1')[1:length(label_sele)], alpha('black',0.1))
		} else {
			mycolor <- c(color_sele, alpha('black',0.1))
		}
		names(mycolor) <- c(label_sele, 'others')
	} else {
		df.plot <- df.plot %>% 
			mutate(group = 'others')
		mycolor <- alpha('black',0.1)
	}
	
	rs = cor(df.plot %>% select(x, y), method = 's')[1,2] 
	rp = cor(df.plot %>% select(x, y), method = 'p')[1,2] 
	
	if (is.null(xy_breaks)){
		xy_min <- min(df.plot %>% select(x, y))
		xy_max <- max(df.plot %>% select(x, y))
		xy_breaks = custom_breaks(xy_min, xy_max)
	} else {
		xy_min <- min(df.plot %>% select(x, y), min(xy_breaks))
		xy_max <- max(df.plot %>% select(x, y), max(xy_breaks))
	}
	
	mag = 1
	line_w = 0.75 / .pt * mag
	font_size = 5.5 * mag
	font_size_offset = 0.5
	gplot <- ggplot(df.plot, aes(x = x, y = y))
	
	if (is.null(label_sele)){
		gplot <- gplot + 
			geom_pointdensity(adjust = 5, shape = 21, size = 0.2, stroke = 0.1) +
			scale_color_viridis(option="A") 
	} else {
		gplot <- gplot +
			geom_point(size = 0.6, shape = 21, color = t_col('black', 0.2), stroke = 0.2) + 
			scale_fill_manual(values = mycolor) 
	}
	gplot <- gplot +
		geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'gray20', linewidth = line_w) +
		labs(x = 'Measured tail-length change (nt)', y = 'Predicted tail-length change (nt)') +
		coord_cartesian(xlim = c(xy_min, xy_max), ylim = c(xy_min, xy_max), clip = 'off') +
		scale_x_continuous(breaks = xy_breaks, labels = signs(xy_breaks, accuracy = 1), expand = expansion(mult = c(0.05, 0.05))) +
		scale_y_continuous(breaks = xy_breaks, labels = signs(xy_breaks, accuracy = 1), expand = expansion(mult = c(0.05, 0.05))) + 
		annotate('text', x = xy_min, y = xy_max, label = deparse(bquote('n ='~.(signs(nrow(df.plot), format = scales::comma)))), parse = T, hjust = 0, size = (font_size - font_size_offset) * 5 / 14) +
		annotate('text', x = xy_min, y = xy_max - (xy_max - xy_min) * 0.07, label = deparse(bquote(italic('R')[s]~'='~.(signs(rs, accuracy = 0.01)))), parse = T, hjust = 0, size = (font_size - font_size_offset) * 5 / 14) +
		annotate('text', x = xy_min, y = xy_max - (xy_max - xy_min) * 0.14, label = deparse(bquote(italic('R')[p]~'='~.(signs(rp, accuracy = 0.01)))), parse = T, hjust = 0, size = (font_size - font_size_offset) * 5 / 14) +
		theme(	
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_blank(),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(hjust = 0.5, size=font_size, color='black'),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black', vjust = 1, margin=margin(1,0,0,0)),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.2),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			
			legend.position = 'none',
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size-font_size_offset,color='black'),
			legend.key.height = unit(0.1, 'inch'),
			legend.key.width = unit(0.1, 'inch'),
			legend.key = element_blank(),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0),
			
			aspect.ratio = 1
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

scatter_plot <- function(df_in, col_x, col_y, xbreaks = NULL, ybreaks = NULL, xlabel = 'x', ylabel = 'y', fn = 'plot'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y))
	pval <- cor.test(df.plot$x, df.plot$y, method = 'p', alternative = 'greater')$p.value
	pv.label <- bquote(italic(P)==.(format(pval, scientific = TRUE, digits = 2)))
	rp <- round(cor(df.plot$x, df.plot$y, method = 'p'),2)
	rs <- round(cor(df.plot$x, df.plot$y, method = 's'),2)
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x)
		xmax <- max(df.plot$x)
		xbreaks <- custom_breaks(xmin, xmax, 6)
	} else {
		xmin <- min(xbreaks, df.plot$x)
		xmax <- max(xbreaks, df.plot$x)
	}
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		xbreaks <- custom_breaks(ymin, ymax, 6)
	} else {
		ymin <- min(ybreaks, df.plot$y)
		ymax <- max(ybreaks, df.plot$y)
	}
	
	font_size <- 5.5
	font_size_offset <- 0.5
	line_w <- 0.75/.pt
	gplot <- ggplot(df.plot, aes(x = x, y = y)) +
		geom_point(size = 1, shape = 21, color = alpha('black', 1), stroke = 0.3) + 
		labs(x = xlabel, y = ylabel) + 
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_x_continuous(breaks = xbreaks, labels = signs(xbreaks)) +
		scale_y_continuous(breaks = ybreaks, labels = signs(ybreaks)) +
		annotate('text', x = xmax - 0.01*(xmax - xmin), y = ymin + 0.02*(ymax - ymin), label = deparse(pv.label), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = xmax - 0.01*(xmax - xmin), y = ymin + 0.1*(ymax - ymin), label = deparse(bquote(italic(R)[s]==.(rs))), parse = T, size = (font_size-font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		annotate('text', x = xmax - 0.01*(xmax - xmin), y = ymin + 0.18*(ymax - ymin), label = deparse(bquote(italic(R)[p]==.(rp))), parse = T, size = (font_size-font_size_offset) * 5 /14, hjust = 1, vjust = 1) +
		theme(
			plot.title = element_blank(),
			plot.background = element_rect(fill = 'transparent',color=NA),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			
			axis.text.y.left = element_text(size=font_size, color='black', margin=margin(t=0,r=1,b=0,l=0)),
			axis.text.y.right = element_text(size=font_size, color='black', margin=margin(t=0,r=0,b=0,l=1)),
			axis.text.x.top = element_text(size=font_size, color='black', margin=margin(t=0,r=0,b=1,l=0)),
			axis.text.x.bottom = element_text(size=font_size, color='black', margin=margin(t=1,r=0,b=0,l=0)),
			axis.line.x = element_line(color='black',size=line_w),
			axis.line.y = element_line(color='black',size=line_w),
			axis.title.x = element_text(size=font_size,color='black', vjust = 1, margin=margin(0.5,0,0,0)),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.5, margin=margin(0,0.5,0,0)),
			axis.ticks.x = element_line(color = "black", size = line_w),
			axis.ticks.y = element_line(color = "black", size = line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			
			legend.position=c(0.05,1),
			legend.justification=c(0,1),
			legend.title = element_blank(),
			legend.text = element_text(size=font_size - font_size_offset,color='black', margin=margin(0.5,0,0.5,-3)),
			legend.background = element_blank(),
			legend.key.height = unit(0.05, 'inch'),
			legend.key.width = unit(0.05, 'inch'),
			legend.key = element_blank(),
			
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0),
			
			aspect.ratio = 1
		)
	CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2, family = 'Arial')
	print(gplot)
	dev.off()
	invisible(df.plot)
}