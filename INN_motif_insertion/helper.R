library('signs')
library('Cairo')
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
		scale_x_continuous(name = paste0('Distance to 3\' ends (nt)'), 
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

pas_cpe_line_plot <- function(df_in, ybreaks = NULL, xbreaks = seq(-100,100,50), fn = 'my'){
	df.plot <- df_in %>% filter(dis2pas <= xbreaks[length(xbreaks)] & dis2pas >= xbreaks[1])
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$avg - df.plot$sem, na.rm = T)
		ymax <- max(df.plot$avg + df.plot$sem, na.rm = T)
		ybreaks = custom_breaks(ymin, ymax, length.out = 8)
	} else {
		ymin <- min(ybreaks, min(df.plot$avg - df.plot$sem, na.rm = T))
		ymax <- max(ybreaks, max(df.plot$avg + df.plot$sem, na.rm = T))
	}
	
	xmin <- min(xbreaks)
	xmax <- max(xbreaks)
	
	mycolor = brewer.pal(9, 'Set1')[1]
	font_size <- 5.5
	font_size_offset = 0.5
	line_w <- 0.75/.pt
	
	gplot <- ggplot(df.plot, aes(x = pos, y = avg)) +
		geom_line(size = 0.75/.pt, color = mycolor) +
		geom_ribbon(aes(ymin=avg-sem, ymax=avg+sem), fill = mycolor, alpha = 0.3, color = NA) +
		scale_y_continuous(name = 'Mean predicted tail-length change (nt)', limits = c(ymin, ymax), 
											 breaks = ybreaks, labels = signs(ybreaks), expand = expansion(mult = c(0.05, 0.05))) +
		scale_x_continuous(name = paste0('Relative position to PAS (nt)'), 
											 breaks = xbreaks, 
											 limits = c(xmin, xmax),
											 labels = signs(xbreaks, format = scales::comma), 
											 expand = expansion(mult = c(0.05,0.05))) +
		annotate('rect', xmin = (-4) , xmax = 6, ymin = ymin, ymax = ymax, alpha = 0.3, fill = brewer.pal(9, 'Set1')[2])+
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin =unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			
			axis.text.y = element_text(size=font_size,color='black'),
			axis.text.x = element_text(size=font_size,color='black', angle = 45, hjust = 1),
			axis.line.x = element_line(color='black',size=line_w),
			axis.line.y = element_line(color='black',size=line_w),
			axis.title.x = element_text(size=font_size,color='black'),
			axis.title.y = element_text(size=font_size,color='black', vjust = 0.2),
			axis.ticks.x = element_line(color = "black", size = line_w),
			axis.ticks.y = element_line(color = "black", size = line_w),
			axis.ticks.length = unit(2, 'pt'),
			
			legend.position='none',
			aspect.ratio = 1
		) 
	#png(paste0(fn, '.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	CairoPDF(paste0(fn, '_line_plot.pdf'), width = 2, height = 2)
	print(gplot)
	dev.off()	
}
