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

r_distribution_violin_plot <- function(df_in, groups, ybreaks = NULL, group_colors = NULL, fn = 'plot', ft = 'pdf'){
	# perform t-test for all combinations
	df_in <- df_in %>% mutate(r = sqrt(r2))
	group_combo <- comboGeneral(unique(df_in$group), m = 2)
	r_test <- lapply(1:nrow(group_combo), function(x){
		x_series <- df_in %>% filter(group == group_combo[x, 1]) %>% pull(r)
		y_series <- df_in %>% filter(group == group_combo[x, 2]) %>% pull(r)
		bind_rows(
			list('group_x' = group_combo[x, 1], 
					 'group_y' = group_combo[x, 2], 
					 'pval' = t.test(x = x_series, y = y_series, alternative = 'less')$p.value
					 ),
			list('group_x' = group_combo[x, 2], 
					 'group_y' = group_combo[x, 1], 
					 'pval' = t.test(x = y_series, y = x_series, alternative = 'less')$p.value
					 )
		) %>% return(.)
	}) %>% bind_rows()
	
	# make a violin plot
	df.plot <- df_in %>% filter(group %in% groups)
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$r)
		ymax <- max(df.plot$r)
		ybreaks <- custom.breaks(ymin, ymax, length.out = 6)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}

	if (is.null(group_colors)){
		if (length(levels(df.plot$group)) > 9){
			group_colors <- viridis(length(levels(df.plot$group)))
		} else {
			group_colors <- brewer.pal(9, 'Set1')[1:length(levels(df.plot$group))]
		}
	} 
	
	line_w = 0.75/.pt
	font_size = 6
	font_size_offset = 0.5
	gplot <- ggplot(data = df.plot, aes(x = group, y = r, fill = group)) + 
		geom_violin(color = NA, width = 0.9, alpha = 0.75) +
		stat_summary(fun.data = custom_quantile, geom = 'boxplot', lwd = line_w, width = 0.5, fill = NA) +
		labs(y = 'Pearson R') +
		scale_y_continuous(breaks = ybreaks, limits = c(ymin, ymax))+
		scale_fill_manual(values = group_colors) +
		coord_flip() +
		theme(
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.title = element_text(size=font_size, face = "bold",color='black', hjust = 0.5, margin=margin(0,0,3,0)),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			panel.background = element_blank(),
			panel.grid.major.x = element_line(color = 'gray50', linewidth=0.25/.pt, linetype = 'dotted'),
			panel.grid.minor.x = element_blank(),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank(),
			
			axis.text.x = element_text(size=font_size,color='black', margin=margin(t=0,r=0.5,b=0,l=0)),
			axis.text.y = element_blank(),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.y = element_blank(),
			axis.title.x = element_text(size=font_size, color='black', margin=margin(0,1,0,0)),
			axis.ticks.y = element_blank(),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.ticks.length=unit(1.5/.pt, "mm"),
			
			legend.position = 'right',
			legend.justification = c(0,1),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size - font_size_offset,color='black',margin=margin(0,0,0,-3)),
			legend.key = element_blank(),
			legend.key.size = unit(font_size, 'points'),
			legend.spacing.y = unit(font_size/2, 'points'),
			legend.box.margin = margin(0,0,0,0),
			legend.margin = margin(0,0,0,0),
			
			aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn,'.pdf'), width = 4, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=4, height=2,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
	invisible(r_test)
}

pval_heatmap <- function(df_in, col_x, col_y, col_val, fn = 'plot', ft = 'pdf'){
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y),
															val = df_in %>% pull(col_val))
	df.plot <- df.plot %>% mutate(log_val = -log10(val)) %>% 
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

violin_plot <- function(df_in, col_x, col_y, col_dodge = NULL, dodge_width = 0.6, label_x = NULL, ybreaks = NULL, fill_color = NULL, fn = 'plot', ft = 'pdf'){
	# make a violin plot
	df.plot <- df_in %>% mutate(x = df_in %>% pull(col_x),
															y = df_in %>% pull(col_y)) %>%
		mutate(x = factor(x))
	
	if (!is.null(col_dodge)){
		df.plot <- df.plot %>% mutate(dodge = df.plot %>% pull(col_dodge)) %>%
			mutate(dodge = factor(dodge))
	}
	
	if (is.null(ybreaks)){
		ymin <- min(df.plot$y)
		ymax <- max(df.plot$y)
		ybreaks <- custom.breaks(ymin, ymax, length.out = 6)
	} else {
		ymin <- min(ybreaks)
		ymax <- max(ybreaks)
	}
	
	if (is.null(fill_color)){
		fill_color <- brewer.pal(9, 'Set1')
	}
	line_w = 0.75/.pt
	font_size = 5.5
	font_size_offset = 0.5
	
	if (is.null(col_dodge)){
		gplot <- ggplot(data = df.plot, aes(x = x, y = y, fill = x)) +
			geom_violin(color = 'NA', width = 0.9, alpha = 0.75) +
			stat_summary(fun.data = custom_quantile, geom = 'boxplot', lwd = line_w, width = 0.5, fill = NA) 
	} else {
		gplot <- ggplot(data = df.plot, aes(x = x, y = y, fill = dodge)) +
			geom_violin(color = NA, width = dodge_width, alpha = 0.75, position = position_dodge(dodge_width)) #+
			#stat_summary(fun.data = custom_quantile, geom = 'boxplot', lwd = line_w, width = dodge_width - 0.1, fill = NA, position = position_dodge(dodge_width)) 
	}
	
	gplot <- gplot +	
		coord_cartesian(ylim = c(ymin, ymax)) +
		labs(y = 'Pearson R', x = label_x) +
		scale_y_continuous(breaks = ybreaks) +
		scale_x_discrete(expand = expansion(add = c(0.75, 1))) +
		scale_fill_manual(values = fill_color) +
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
			legend.justification = c(1,0),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.text = element_text(size=font_size - font_size_offset,color='black',margin=margin(0,0,0,0)),
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
