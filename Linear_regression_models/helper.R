library('tidyverse')
library('ggplot2')
library('Cairo')
library('RcppAlgos')
library('RColorBrewer')
library('viridis')
library('signs', quietly = T)
library('ggpointdensity', quietly = T)

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
			axis.text.y = element_text(size=font_size,color='black', margin=margin(t=0,r=0.5,b=0,l=0)),
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
			legend.margin = margin(0,0,0,0)
			
			#aspect.ratio = 1
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn,'.pdf'), width = 2.5, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2.5, height=2,unit='in',res=600,bg='transparent')
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


# make a scatter plot for comparing the predicted values and measured values
prediction_w_sele_scatter_plot <- function(df_in, gene_id = 'id', xaxis = 'y', yaxis = 'y_pred', xy_breaks = NULL, label_sele = NULL, color_sele = NULL, fn = 'plot', ft = 'pdf'){
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
		xy_breaks = custom.breaks(xy_min, xy_max)
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
		geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'gray20', size = line_w) +
		labs(x = 'Measured tail-length change (nt)', y = 'Predicted tail-length change (nt)') +
		coord_cartesian(xlim = c(xy_min, xy_max), ylim = c(xy_min, xy_max), clip = 'off') +
		scale_x_continuous(breaks = xy_breaks, labels = signs(xy_breaks, accuracy = 1), expand = expansion(mult = c(0.05, 0.05))) +
		scale_y_continuous(breaks = xy_breaks, labels = signs(xy_breaks, accuracy = 1), expand = expansion(mult = c(0.05, 0.05))) + 
		annotate('text', x = xy_min, y = xy_max, label = deparse(bquote('n ='~.(signs(nrow(df.plot), format = scales::comma)))), parse = T, hjust = 0, size = (font_size - font_size_offset) * 5 / 14) +
		annotate('text', x = xy_min, y = xy_max - (xy_max - xy_min) * 0.07, label = deparse(bquote(italic('R')[s]~'='~.(format(rs, digits = 2)))), parse = T, hjust = 0, size = (font_size - font_size_offset) * 5 / 14) +
		annotate('text', x = xy_min, y = xy_max - (xy_max - xy_min) * 0.14, label = deparse(bquote(italic('R')[p]~'='~.(format(rp, digits = 2)))), parse = T, hjust = 0, size = (font_size - font_size_offset) * 5 / 14) +
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

# make a bar plot for features
coefficient_bar_plot <- function(df_in, x = 'kmer', y = 'coefficient', y_std = 'std', top_n = 20, ybreaks = NULL, color_sele = NULL, fn = 'plot', ft = 'pdf'){
	df_in$y_axis <- df_in %>% pull(!!as.symbol(y))
	df_in$x_axis <- df_in %>% pull(!!as.symbol(x))
	df_in$y_axis_std <- df_in %>% pull(!!as.symbol(y_std))
	
	df.plot <- df_in %>% filter(x_axis != 'intercept') %>% 
		arrange(desc(abs(y_axis))) %>% mutate(group = ifelse(y_axis > 0, 'g1', 'g2')) %>%
		dplyr::slice(1:min(top_n, nrow(.))) %>% arrange(abs(y_axis)) %>%
		mutate(x_axis = gsub('T', 'U', x_axis)) %>%
		mutate(x_axis = factor(x_axis, levels = .$x_axis)) 
	
	#my_label <- paste0(df.plot %>% pull(x_axis) %>% str_split_fixed(., '_', 2) %>% {.[,2]}, ' (',
	#									 df.plot %>% pull(x_axis) %>% str_split_fixed(., '_', 2) %>% {.[,1]}, ')') 
	
	ymax = df.plot %>% pull(y_axis) %>% max(., 0)
	ymin = df.plot %>% pull(y_axis) %>% min(., 0)
	
	if(is.null(ybreaks)){
		ybreaks = custom.breaks(-ymax, ymax, length.out = 7, digits = 2)
	} else{
		ymax <- max(ymax, ybreaks)
		ymin <- min(ymin, ybreaks)
	}
	
	if(is.null(color_sele)){
		mycolor = c(brewer.pal(11, 'PuOr')[9], brewer.pal(8, 'Set2')[1])
	} else{
		mycolor = color_sele
	}
	
	line_w = 0.75/.pt
	font_size = 5.5
	gplot <- ggplot(df.plot, aes(x = x_axis, y = y_axis, fill = group)) +
		geom_bar(stat = 'identity', position = position_dodge(), width = 0.5) + 
		geom_errorbar(aes(ymin = y_axis - y_axis_std, ymax = y_axis + y_axis_std), linewidth = 0.2, width = 0.4) +
		labs(y = 'Model coefficients') +
		scale_y_continuous(position = 'right', breaks = ybreaks, labels = signs(ybreaks, format = scales::comma), expand = expansion(mult = 0.05)) +
		scale_x_discrete(expand = expansion(add = 0.6)) +
		scale_fill_manual(values = mycolor) +
		geom_hline(yintercept = 0, linewidth = line_w) +
		coord_flip(expand = T, clip = 'off', ylim = c(ymin, ymax)) +
		theme(
			legend.position='none',
			
			axis.title.x = element_text(size=font_size,color='black',vjust=0.5),
			axis.text.x = element_text(size=font_size,color='black'),
			axis.ticks.x = element_line(color='black',linewidth=line_w),
			axis.line.x = element_line(color='black',linewidth=line_w),
			
			axis.title.y = element_blank(),
			axis.text.y = element_text(size=font_size,color='black'),
			axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			
			plot.title = element_blank(), 
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), 'inch'),
			panel.background = element_blank(),
			panel.grid = element_blank()
		)
	if (ft == 'pdf'){
		CairoPDF(paste0(fn,'.pdf'), width = 2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2, height=3,unit='in',res=600,bg='transparent')
	}
	print(gplot)
	dev.off()
	invisible(df.plot)
}
