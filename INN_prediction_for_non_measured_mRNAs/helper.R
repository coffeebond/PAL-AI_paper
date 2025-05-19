library('ggplot2')
library('Cairo')
library('ggrepel')
library('RColorBrewer')
library('viridis', quietly = T)
library('ggpointdensity', quietly = T)
library('signs', quietly = T)
library('circlize')
library('ComplexHeatmap')
library('gprofiler2')

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

plot_heatmaps <- function(mat_in, info_df, pred_mat = NULL, legend_title = "Tail-length change", fn = 'heatmap', bound = Inf, main_fill = NULL, mybreaks = NULL, row_split = NULL, show_row_names = F, show_column_names = F, column_split = NULL, column_labels = NULL, g_fsize = 7, anno_width = 0.08, f_width = 2, f_height = 3, f_type = 'pdf') {	
	
	row_order = info_df %>% pull(gene_name)
	
	# get the breaks for tail-length heatmap
	
	# get the row categories and order
	if (!is.null(row_split)){
		row_split <- info_df %>% pull(!!as.symbol(row_split))
		info_df$row_cat <- info_df %>% pull(!!as.symbol(row_split))
		n_cat_labels <- paste0('n = ', info_df %>% group_by(row_cat) %>% summarize(n=n(), .groups = 'drop') %>% pull(n))
	} else{
		n_cat_labels <- NULL
	}
	
	# set the value limits and breaks
	max_min = max(min(max(mat_in, na.rm = T), bound), min(abs(min(mat_in, na.rm = T)), bound))
	if (is.null(mybreaks)){
		my_min <- min(mat_in, na.rm = T)
		my_max <- max(mat_in, na.rm = T)
		mybreaks <- custom_breaks(my_min, my_max, length.out = 7)
	} 
	mylabels <- signs(mybreaks, accuracy = 1)
	if ((bound <= max(mat_in, na.rm = T)) | (-bound >= min(mat_in, na.rm = T))){
		mybreaks <- mybreaks[mybreaks <= bound & mybreaks >= (-bound)]
		mylabels <- signs(mybreaks, accuracy = 1)
		if (bound <= max(mat_in, na.rm = T)){
			mybreaks <- unique(c(mybreaks, bound))
			mylabels[length(mylabels)] <- paste0('\u2265', mylabels[length(mylabels)])
		}
		if (-bound >= min(mat_in, na.rm = T)) {
			mybreaks <- unique(c(-bound, mybreaks))
			mylabels[1] <- paste0('\u2264', mylabels[1])
		}
	} 
	
	# add prediction layer
	if (!is.null(pred_mat)){
		my_layer_func = function(j, i, x, y, width, height, fill) {
			v = pindex(pred_mat, i, j)
			grid.points(x[v], y[v], pch = 16, gp = gpar(col = 'black'), size = unit(g_fsize/2, "pt"))
			#grid.rect(x[v], y[v], width = unit(1/ncol(mat_in), "npc"), height = unit(1/nrow(mat_in), "npc"), gp = gpar(lwd = 1/.pt, fill = "transparent"))
		}
	} else {
		my_layer_func = NULL
	}
	
	# fill color
	if (is.null(main_fill)){
		main_fill <- colorRamp2(c(-max_min, 0, max_min), c(brewer.pal(9, 'Set1')[2], 'white', brewer.pal(9, 'Set1')[1]))
	} else {
		main_fill <- colorRamp2(c(-max_min, 0, max_min), c(brewer.pal(11, 'PiYG')[11], 'white', brewer.pal(11, 'PiYG')[1]))
	}
	
	# column labels
	if (is.null(column_labels)){
		column_labels = colnames(mat_in)
	}
	ht = Heatmap(mat_in %>% {.[row_order,]}, 
							 col = main_fill,
							 column_title = NULL,
							 column_title_side = c("bottom"),
							 column_title_gp = gpar(fontsize = g_fsize),
							 column_order = colnames(mat_in),
							 column_names_gp = gpar(fontsize = g_fsize, fontface = "italic"),
							 column_names_rot = 0,
							 column_names_centered = TRUE,
							 column_labels = column_labels,
							 column_split = column_split, 
							 show_column_names = show_column_names,
							 show_row_names = show_row_names,
							 row_names_gp = gpar(fontsize = g_fsize),
							 row_names_side = c("left"),
							 row_order = row_order,
							 row_split = row_split,
							 row_title = n_cat_labels,
							 row_title_gp = gpar(fontsize = g_fsize),
							 row_title_rot = 90,
							 row_gap = unit(2, "pt"),
							 border = TRUE,
							 border_gp = gpar(lwd = 1/.pt),
							 rect_gp = gpar(col = "white", lwd = 1/.pt), # add white lines among blocks
							 width = unit(anno_width*1.5*ncol(mat_in), 'inch'),
							 use_raster = F,
							 na_col = "gray50",
							 layer_fun = my_layer_func,
							 heatmap_legend_param = list(
							 	title = legend_title,
							 	title_gp = gpar(fontsize = g_fsize),
							 	title_position = "leftcenter-rot",
							 	#col_fun = colorRamp2(c(min(mat), 0, max(mat)), c(brewer.pal(9, 'Set1')[2], 'white', brewer.pal(9, 'Set1')[1])),
							 	at = mybreaks,
							 	labels = mylabels,
							 	labels_gp = gpar(fontsize = g_fsize),
							 	grid_width = unit(anno_width/1.25, "inch"),
							 	legend_height = unit(1, "inch")
							 )
	)
	
	if (f_type =='pdf'){
		CairoPDF(file=paste0(fn,'.pdf'), width=f_width, height=f_height)
	} else {
		png(file=paste0(fn,'.png'), width=f_width, height=f_height, unit='in',res=600,bg='transparent')
	}
	pushViewport(viewport(gp = gpar(lwd = 1/.pt)))
	draw(ht, 
			 legend_gap = unit(0.45, "inch"), ht_gap = unit(0.05, "inch"), 
			 merge_legend = TRUE,#annotation_legend_list = list(ht_mm_lgd),
			 background = 'transparent', newpage = FALSE)
	dev.off()
	invisible(info_df)
}

gost_bar_plot <- function(df_in, padj_cutoff = 0.05, mybreaks = NULL, mycolor = alpha(brewer.pal(9,'Set1')[1], 0.3), ft = 'pdf', fn = 'plot'){
	# a function to make a bar plot from gprofiler results
	dplot <- df_in %>% filter(p_value < padj_cutoff) %>% 
		mutate(log10p = -log10(p_value)) %>% arrange(desc(log10p))
	if (nrow(dplot) == 0){
		warning('No significant GO categories!')
		return(-1)
	}
	dplot$term_name <- factor(dplot$term_name, levels = rev(dplot$term_name))
	
	ymin <- 0
	ymax <- max(dplot$log10p)
	if(is.null(mybreaks)){
		mybreaks = custom_breaks(ymin, ymax, length.out = 6)
	}
	
	font_size = 5.5
	font_size_offset = 1
	line_w = 0.75/.pt
	gplot <- ggplot(dplot, aes(x = term_name, y = log10p)) +
		geom_bar(stat = 'identity', fill = mycolor) +
		labs(y = bquote('\u2212log'[10]*'('*italic(P)[adjust]*')')) +
		geom_text(aes(label = term_name), y = ymin + (ymax - ymin)*0.05, hjust = 0, size = (font_size-font_size_offset)*5/14) +
		scale_y_continuous(breaks = mybreaks) +
		coord_flip() +
		theme(
			legend.position = 'none',	
			plot.background = element_rect(fill = 'transparent',color=NA),
			plot.margin = unit(c(0.01,0.01,0.01,0.01), "inch"),
			plot.title = element_blank(),
			panel.grid = element_blank(),
			panel.background = element_blank(),
			axis.text.y = element_blank(),#element_text(size=2.5,color='black'),
			axis.text.x = element_text(hjust = 0.5, size=font_size, color='black'),
			axis.line.x = element_line(color='black',linewidth=line_w),
			axis.line.y = element_line(color='black',linewidth=line_w),
			axis.title.x = element_text(size=font_size,color='black', vjust = 1, margin=margin(1,0,0,0)),
			axis.title.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.ticks.x = element_line(color = "black", linewidth = line_w)
		)
	
	#myheight<- nrow(dplot) * 0.3 + 0.3
	if (ft == 'pdf'){
		CairoPDF(paste0(fn, '.pdf'), width = 2, height = 2, family = 'Arial')
	} else {
		png(file=paste0(fn,'.png'), width=2, height=2,unit='in',res=600,bg='transparent')
	}
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
		vjust_slp <- 0
	} else {
		y_slp_1 <- ymax - 0.18*(ymax - ymin)
		y_slp_2 <- ymax - 0.1*(ymax - ymin)
		y_slp_3 <- ymax - 0.02*(ymax - ymin)
		vjust_slp <- 1
	}
	if (grepl('left', stat_legend_position)){
		x_slp <- xmin + 0.01*(xmax - xmin)
		hjust_slp <- 0
	} else {
		x_slp <- xmax - 0.01*(xmax - xmin)
		hjust_slp <- 1
	}
	
	gplot <- gplot +
		labs(x = labs[1], y = labs[2]) +
		coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
		scale_x_continuous(breaks = xbreaks, label = signs(xbreaks), expand = expansion(mult = pad)) +
		scale_y_continuous(breaks = ybreaks, label = signs(ybreaks), expand = expansion(mult = pad)) +
		annotate('text', x = x_slp, y = y_slp_1, label = deparse(bquote(n==.(signs(n, format = scales::comma)))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = hjust_slp, vjust = vjust_slp) +
		annotate('text', x = x_slp, y = y_slp_2, label = deparse(bquote(italic(R)[s]==.(rs))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = hjust_slp, vjust = vjust_slp) +
		annotate('text', x = x_slp, y = y_slp_3, label = deparse(bquote(italic(R)[p]==.(rp))), parse = T, size = (font_size - font_size_offset) * 5 /14, hjust = hjust_slp, vjust = vjust_slp) +
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

