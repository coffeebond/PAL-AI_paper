library('ggplot2')
library('Cairo')
library('RColorBrewer')
library('signs', quietly = T)

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
line_plot <- function(df_in, col_x, col_y, col_y_delta, col_group, xbreaks = NULL, ybreaks = NULL, y_title = 'y', x_title = 'x', fn = 'my'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y),
															y_delta = df_in %>% pull(col_y_delta),
															group = df_in %>% pull(col_group)) %>%
		mutate(group = factor(group))
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y - df.plot$y_delta, na.rm = T)
		ymax <- max(df.plot$y + df.plot$y_delta, na.rm = T)
		ybreaks = custom_breaks(ymin, ymax, length.out = 8)
	} else {
		ymin <- min(ybreaks, min(df.plot$y - df.plot$y_delta, na.rm = T))
		ymax <- max(ybreaks, max(df.plot$y + df.plot$y_delta, na.rm = T))
	}
	
	if (is.null(xbreaks)){
		xmin <- min(df.plot$x)
		xmax <- max(df.plot$x)
		xbreaks = custom_breaks(xmin, xmax, length.out = 8)
	} else {
		xmin <- min(xbreaks, df.plot$x)
		xmax <- max(xbreaks, df.plot$x)
	}
	
	mycolor = brewer.pal(9, 'Set1')[1:length(levels(df.plot$group))]
	font_size <- 5.5
	font_size_offset = 0.5
	line_w <- 0.75/.pt
	
	gplot <- ggplot(df.plot, aes(x = x, y = y, color = group)) +
		geom_line(size = line_w) +
		geom_ribbon(aes(ymin = y - y_delta, ymax = y + y_delta, fill = group), alpha = 0.3, color = NA) +
		scale_color_manual(name = NULL, values = mycolor) +
		scale_fill_manual(name = NULL, values = mycolor) +
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_y_continuous(name = y_title, breaks = ybreaks, labels = signs(ybreaks), expand = expansion(mult = c(0.05, 0.05))) +
		scale_x_continuous(name = x_title, breaks = xbreaks, labels = signs(xbreaks), expand = expansion(mult = c(0.05,0.05))) +
		guides(color = guide_legend(byrow = TRUE)) +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black'),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black'),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.2),
			axis.ticks.x = element_line(color = "black", linewidth = line_w),
			axis.ticks.y = element_line(color = "black", linewidth = line_w),
			axis.ticks.length = unit(2, 'pt'),
			
			legend.position = c(1,1),
			legend.justification = c(1, 1),
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
