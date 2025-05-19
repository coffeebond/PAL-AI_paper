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

motif_density_line_plot <- function(df_in, col_x, col_y, col_group, xbreaks = NULL, xlabels = NULL, ybreaks = NULL, x_lim = NULL, x_title = 'x', y_title = 'y', group_color = NULL, fn = 'my'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y),
															group = df_in %>% pull(col_group)) %>%
		mutate(group = factor(group))
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y, na.rm = T)
		ymax <- max(df.plot$y, na.rm = T)
		ybreaks = custom_breaks(ymin, ymax, length.out = 8)
	} else {
		ymin <- min(ybreaks, min(df.plot$y, na.rm = T))
		ymax <- max(ybreaks, max(df.plot$y, na.rm = T))
	}
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x, na.rm = T)
		xmax <- max(df.plot$x, na.rm = T)
		ybreaks = custom_breaks(xmin, xmax, length.out = 8)
	} else {
		xmin <- min(xbreaks, min(df.plot$x, na.rm = T))
		xmax <- max(xbreaks, max(df.plot$x, na.rm = T))
	}
	
	if (is.null(x_lim)){
		x_lim = c(xmin, xmax)
	}
	
	if (is.null(group_color)){
		group_color <- brewer.pal(9, 'Set1')[1:length(levels(df.plot$group))]
	}

	font_size <- 5.5
	font_size_offset = 0.5
	line_w <- 0.75/.pt
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, color = group)) +
		geom_line(linewidth = 0.75/.pt) +
		scale_color_manual(name = NULL, values = group_color) +
		scale_y_continuous(name = y_title, 
											 limits = c(ymin, ymax), 
											 breaks = ybreaks, 
											 labels = signs(ybreaks), expand = expansion(mult = c(0.05, 0.05))) +
		scale_x_continuous(name = x_title, 
											 breaks = xbreaks, 
											 limits = x_lim,
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
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
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
	CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2)
	print(gplot)
	dev.off()	
	invisible(df.plot)
}
