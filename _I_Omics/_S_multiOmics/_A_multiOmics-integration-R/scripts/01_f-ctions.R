####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####
# f-ctions
####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

# .rs.restartR()

# install.packages("installr")
# library(installr)

# updateR()

# update.packages(checkBuilt=TRUE)

# packageVersion("rlang")
# remove.packages("rlang")
# devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
# packageVersion("rlang")
# update.packages("rlang")


####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

# rm(list = ls(all = TRUE))
gc()



library(scales)
library(ggplot2)
# packageVersion("ggplot2")
library(magrittr)
library(crayon)

library(rstatix)
library(ggpubr)

library(Compositional)
library(equalCovs)

library(psych)

library(corrplot) 

library(pdist)

library(RColorBrewer)

library(plotly)

library(heatmaply)

library(Hmisc)

library(caret)

library(htmlwidgets)

# library(RcmdrMisc)
# packageVersion("car")

####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####
par(mfrow = c(2,3))
n = 9
# (intervals = cut(c(-1,1), n*2))
intervals = cut(c(-1,1), 
                breaks = setdiff(seq(-1, 1, 0.1), 0), 
                include.lowest = TRUE)

# my.heatmap.col1 = c(rev(brewer.pal(n, 'Reds')), brewer.pal(n, 'Blues'))
my.heatmap.col1 = c(rev(brewer.pal(n, 'Reds')), 'white', brewer.pal(n, 'Blues'))
pie(rep(1, length(my.heatmap.col1)), 
    col = my.heatmap.col1, 
    labels = levels(intervals), 
    main = 'heatmap palette 1')

my.heatmap.col2 = c(my.heatmap.col1[1:3], 
                    rep('grey90', length(4:(n*2-2))),
                    my.heatmap.col1[(n*2-1): (n*2 + 1)])
pie(rep(1, length(my.heatmap.col1)), 
    col = rev(my.heatmap.col2), 
    labels = levels(intervals), 
    main = 'heatmap palette 2')



my.ggplot.palette = brewer.pal(9, "Set1")[c(3, 1, 5, 8, 2, 4)]
pie(rep(1, length(my.ggplot.palette)), 
    col = my.ggplot.palette, 
    main = 'ggplot palette',
    labels = c('C', 'H', 'D', 'HD', 'HDW', 'W'))

my.dist.col = brewer.pal(11, "RdYlGn")
pie(rep(1, length(my.dist.col)), 
    col = rev(my.dist.col), 
    main = 'abs cor dist col',
    labels = seq(0, 2, 0.2))
my.dist.col2 = my.dist.col
my.dist.col2[9:11] = 'grey90'
pie(rep(1, length(my.dist.col)), 
    col = rev(my.dist.col2), 
    main = 'abs cor dist col2',
    labels = seq(0, 2, 0.2))

time.levels = c(1, 7, 8, 14, 15, 21, 28)
plot(time.levels, 
     c(rep(0, 2), rep(1, 2), rep(2, 2), 3), 
     pch = c(1, 16, 0, 15, 2, 17, 8), 
     col = 'black', cex = 2,
     yaxt='n', xaxt='n', ann=FALSE, bty='n')
axis(side = 1, at = time.levels, labels = time.levels)

par(mfrow = c(1,1))

NMDS.palette = brewer.pal(9, 'Paired')[c(3,4,8,6,1,2,9)]
show_col(brewer_pal(palette = "Paired")(9)[c(3,4,8,6,1,2,9)])

####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

## https://stackoverflow.com/questions/38068774/rstudio-suddenly-stopped-showing-plots-in-the-plot-pane

dev.cur()
i= as.numeric(dev.cur())
dev.off(i) #where i = index of device to be switched off

# get the device back with

# getOption("device")

# or

dev.set(which = dev.next())


####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

myplot.ggplot <- function(df, x, y, ncol, scales.y, title, xlab, ylab, palette){

  ggplot(df, aes(x = x,
                   y = y,
                   colour = x,
                   pch = x)) +
    geom_point(size=2, shape=16) +
    facet_wrap( ~ variable, ncol = ncol, #6, 
                drop = FALSE, scales = scales.y #"free_y"
                ) +
    theme_bw() +
    theme(legend.position = "right",
          legend.title = element_blank(),
                         # element_text(color = "black",
                         #             size = 10,
                         #            face = 2),
          legend.background = element_rect(fill = "grey90", 
                                           colour = 1)) +
    scale_shape_manual(name = "Point",
                       values = c(15, 19)) +
    labs(x = xlab, y = ylab) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = palette) +
    ggtitle(title) 
}



####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

# https://r-coder.com/correlation-plot-r/

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  #on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- (cor(x, y, use = 'na.or.complete')) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * 0.75) #Cor) # Resize the text by level of correlation
}

# Scatterplot Matrices
my.pairs <- function(df, main, cex.labels, legend, col, ncol) {
  cat(blue('Scatterplot Matrices\n\n'))
  pairs(df,                          # Data frame of variables
        upper.panel = panel.cor,    # Correlation panel
        lower.panel = panel.smooth, # Smoothed regression lines
        pch = 21, 
        bg = col,                   # Background color of the symbol (pch 21 to 25)
        col = col,                  # Border color of the symbol
        main = main,               # Title of the plot
        gap = 1,                    # Distance between subplots
        cex.labels = cex.labels,          # Size of the diagonal text
        font.labels = 1,            # Font style of the diagonal text
        oma = c(10,5,5,5))
  par(xpd = TRUE)
  legend("bottomright", 
         ncol = n, 
         fill = col, 
         legend = legend, 
         bg = 'grey90',
         cex = 0.75, bty="n", box.lwd = 0,
         xpd=TRUE, 
         text.width = rep(0.015, n),
         inset = c(0.25, -0.25)  
         #x.intersp=0,
         #xjust=0, yjust=0
  )
}




# library(psych)
# corr.test {psych}	R Documentation
# Find the correlations, sample sizes, and probability values between elements of a matrix or data.frame.

my.corr.test <- function(df, method, main, n, cex){
  
  cat(blue('Find the correlations, sample sizes, and probability values 
  between elements of a matrix or data.frame\n\n'))
  
  cts = corr.test(x = df, 
                  y = NULL, 
                  use = "pairwise",
                  method = method,
                  adjust = "holm", 
                  alpha = .05,
                  ci = TRUE,
                  minlength = 15, # What is the minimum length for abbreviations. Defaults to 5.
                  normal = FALSE)
  # cat('names(object)', names(cts), '\n\n')
  # print(corr.p(cts$r, n = ncol(cts$r)),short=FALSE)
  # print(cts$stars, quote=FALSE)
  
  # cex = abs(cts$r)
  # cex[upper.tri(cex, diag=TRUE)] = 0# NA
  # cex = as.vector(cex[!is.na(cex)])
  # not working as expected
  
  cor.plot(cts$r, 
           cex = cex, 
           upper = FALSE, 
           diag = FALSE, 
           n = n*2, # The number of levels of shading to use. Defaults to 51
           adjust = "holm", 
           p = 0.05,  
           alpha = 0.75,
           gr = colorRampPalette(c("darkblue", "white", "darkred")),
           stars = TRUE,
           las = 2,
           main = main,
           scale = TRUE)
  
  
}


# SPLOM, histograms and correlations for a data matrix
# pairs.panels(d,
#              smooth = TRUE,      # If TRUE, draws loess smooths
#              scale = TRUE,      # If TRUE, scales the correlation text font
#              density = TRUE,     # If TRUE, adds density plots and histograms
#              ellipses = TRUE,    # If TRUE, draws ellipses
#              method = "pearson", # Correlation method (also "spearman" or "kendall")
#              pch = 21,           # pch symbol
#              lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
#              cor = TRUE,         # If TRUE, reports correlations
#              jiggle = FALSE,     # If TRUE, data points are jittered
#              factor = 2,         # Jittering factor
#              hist.col = 4,       # Histograms color
#              stars = TRUE,       # If TRUE, adds significance level with stars
#              ci = TRUE,          # If TRUE, adds confidence intervals
#              main = 'just proof of drawing concept')

my.pairs.panels <- function(df, method, palette, scale, main, cex, cex.labels){
  cat(blue('SPLOM, histograms and correlations for a data matrix\n\n'))
  pairs.panels(df,
               smooth = TRUE,      # If TRUE, draws loess smooths
               scale = scale, #TRUE,      # If TRUE, scales the elation text font
               density = TRUE,     # If TRUE, adds density plots and histograms
               ellipses = FALSE,    # If TRUE, draws ellipses
               method = method, # elation method (also "spearman" or "kendall")
               pch = 21,           # pch symbol
               lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
               cor = TRUE,         # If TRUE, reports elations
               jiggle = FALSE,     # If TRUE, data points are jittered
               factor = 2,         # Jittering factor
               hist.col = 4,       # Histograms color
               stars = TRUE,       # If TRUE, adds significance level with stars
               ci = TRUE,          # If TRUE, adds confidence intervals
               bg = palette,
               col = palette,
               cex.labels = cex.labels, # 1, # diagonal
               cex = cex, #4,
               main = main) 
  
}


# corCi {psych}	R Documentation
# Bootstrapped and normal confidence intervals for raw and composite correlations




# rcorr {Hmisc}	R Documentation
# Matrix of Correlations and P-values
# rcorr Computes a matrix of Pearson's r or Spearman's rho rank correlation coefficients 
# for all possible pairs of columns of a matrix. Missing values are deleted 
# in pairs rather than deleting all rows of x having any missing variables. 
# Ranks are computed using efficient algorithms (see reference 2), using midranks 
# for ties.

# heatmaply {heatmaply}
# Cluster heatmap based on plotly

my.heatmaply.cor <- function(df, type, dist_method, hclust_method, main, colors){
  cat(blue('Matrix of Correlations and P-values with Cluster heatmap based on plotly\n\n'))
  # r = cor(d, use = 'na.or.complete')
  # Hmisc::rcorr		Matrix of Correlations and P-values
  r.rcorr = rcorr(as.matrix(df), type = type)
  r = r.rcorr$r
  # range(r - r.rcorr$r, na.rm = TRUE)
  p = r.rcorr$P
  
  ind = which(p==0, arr.ind=TRUE)
  zeroes = min(p[p != 0], na.rm = TRUE)/2
  # cbind(colnames(p)[ind[,1]], colnames(p)[ind[,2]])
  # p[ind]
  p[ind] = zeroes
  p[row(p)==col(p)] = zeroes
  # p[ind]
  
  # Cluster heatmap based on plotly
  heatmaply_cor(
    r,
    node_type = "scatter",
    point_size_mat = -log10(p), 
    point_size_name = "-log10(p-value)",
    label_names = c("x", "y", "Correlation"),
    dist_method = dist_method, 
    hclust_method = hclust_method, 
    main = main,
    colors = colors
  )
  
}



# corrplot {corrplot}
# A visualization of a correlation matrix.

my.cor.plot <- function(mymat, main, col, order, type, choose){
  cat(blue('A visualization of a correlation matrix\n\n'))
  M = NULL
  # https://stackoverflow.com/questions/40340081/return-a-matrix-with-ifelse
  M = ifelse(choose == 1, 
             list(cor(mymat, use = 'na.or.complete', method = type)),
             list(rcorr(mymat, type = type)$r)
             )
  # colnames(M) 
  corrplot(M[[1]], 
           method = "color", 
           main = main, 
           mar = c(6, 0, 2, 0), 
           hclust.method = 'ward.D2',
           tl.cex = 0.75,
           tl.col = "grey45", 
           col = col, 
           order = order)
  
}



# pdist {pdist}	R Documentation
# Partitioned Distances
# Description
# Compute distance matrix between two matrices of observations, or two subsets of one matrix




####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####
# pheatmap {pheatmap}	R Documentation
# A function to draw clustered heatmaps
my.pheatmap <- function(df, main, col, breaks, silent){
  pheatmap::pheatmap(df, 
                     method = "color", 
                     main = main, 
                     mar = c(1, 1, 1, 1),
                     tl.cex = 0.75,
                     tl.col = "grey45", 
                     col = col,
                     cluster_cols=FALSE, 
                     cluster_rows=FALSE,
                     cellwidth = 10, 
                     cellheight = 10,
                     breaks = breaks,
                     silent = silent
  )
}


####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

my.customised.t.test <- function(data, var.levels, stress.levels, 
                                 plot.violin, plot.box, plot.dot,
                                 y.lab, single,
                                 p.cex.labels, p.cex, p.palette){
  
  mydata = dplyr::as_tibble(data.table::data.table(data))
  mydata.long <- mydata %>%
    tidyr::pivot_longer(-Stress, names_to = "variables", values_to = "value")
  
  mydata.long$variables = factor(mydata.long$variables, levels = var.levels)
  mydata.long$Stress = factor(mydata.long$Stress, levels = stress.levels)
  
  mydata.long$value = as.numeric(mydata.long$value)
  
  mydata.long %>%
    dplyr::group_by(variables, Stress) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(value),
      sd = sd(value)
    ) %>%
    dplyr::ungroup()
  
  stat.test <- mydata.long %>%
    dplyr::group_by(variables) %>%
    rstatix::t_test(value ~ Stress, p.adjust.method = "holm")
  # Remove unnecessary columns and display the outputs
  stat.test %>% dplyr::select(-.y., -statistic, -df)
  stat.test <- stat.test %>% rstatix::add_xy_position(x = "Stress", 
                                                      fun = "max", # "mean_sd"
                                                      step.increase = 0)
  
  if(plot.violin){
    # https://mode.com/blog/violin-plot-examples/
    myplot <- ggpubr::ggdotplot(
      mydata.long, x = "Stress", y = "value",
      binwidth=0.1, 
      stackgroups = FALSE, # add = "mean_se",
      fill = "Stress", 
      palette = p.palette, #"npg", 
      legend = "none",
      ylab = y.lab, 
      add = c("violin", "jitter", "mean_sd"), #"dotplot", 
      title = 'Input data',
      # breaks = seq(0,4,0.5), 
      ggtheme = ggpubr::theme_pubr(border = TRUE)+ 
        theme( axis.text = element_text( angle = 90 ))## +
      # geom_point(pch = 5, position = position_jitterdodge())
    ) + 
      facet_wrap(~variables, ncol = length(unique(mydata.long$variables))) + #, scales="free_y") + 
      scale_y_continuous(trans = log2_trans(),
                         labels =  trans_format("log2", math_format(2^.x)),
                         n.breaks = 4) + 
      theme( axis.text = element_text( angle = 90 ))
    
    print(myplot)
  }
  
  if(plot.box){
    myplot <- ggpubr::ggboxplot( 
      mydata.long, x = "Stress", y = "value",
      cex = 0.5, 
      stackgroups = FALSE, 
      add = c("mean_sd", "jitter"),
      title = 'Input data',
      fill = "Stress", 
      palette = p.palette, #"npg", 
      legend = "none",
      ylab = y.lab, 
      ggtheme = ggpubr::theme_pubr(border = TRUE)+ 
        theme( axis.text = element_text( angle = 90 ))## +
      # geom_point(pch = 5, position = position_jitterdodge())
    ) + 
      facet_wrap(~variables, 
                 ncol = ceiling(length(unique(mydata.long$variables))/2), 
                 scales="free_y") + 
      scale_y_continuous(#trans = log2_trans(),
                         #labels =  trans_format("log2", math_format(2^.x)),
                         n.breaks = 4) +
      theme( axis.text = element_text( angle = 90 ))
    
    print(myplot) 
  }
  
  if(plot.dot){
    myplot <- ggpubr::ggdotplot(
      mydata.long, x = "Stress", y = "value",
      stackgroups = FALSE, # add = "mean_se",
      fill = "Stress", 
      palette = p.palette, #"npg", 
      legend = "none",
      ylab = y.lab, 
      add = c("jitter", "mean_sd"), #"dotplot", 
      title = 'Input data',
      # breaks = seq(0,4,0.5), 
      ggtheme = ggpubr::theme_pubr(border = TRUE)+ 
        theme( axis.text = element_text( angle = 90 )),
        position = position_jitter(0.05)## +
      # geom_point(pch = 5, position = position_jitterdodge())
    ) + 
      facet_wrap(~variables, 
                 scales = "free", 
                 ncol = ceiling(length(unique(mydata.long$variables))/2)) + 
      scale_y_continuous(#trans = log2_trans(),
                         #labels =  trans_format("log2", math_format(2^.x)),
                         n.breaks = 4) +
      theme( axis.text = element_text( angle = 90 ))
    
    print(myplot)
  }
  
  # if(nrow(stat.test.sig) > 0)
  
  stat.test.sig = stat.test[stat.test$p.adj.signif != 'ns', ]
  # print(stat.test.sig %>% dplyr::arrange(group2, group1))
  
  graphs <- mydata.long %>%
    dplyr::group_by(variables) %>%
    rstatix::doo(
      ~ggpubr::ggdotplot(
        data =., 
        x = "Stress", 
        y = "value",
        fill = "Stress", 
        palette = p.palette, #"npg", 
        legend = "none",
        add = c("jitter", "mean_sd"),
        # ggtheme = ggpubr::theme_pubr(border = TRUE) + 
        #   theme(axis.text = element_text( angle = 90 )),
        position = position_jitter(0.05),
        ggtheme = ggpubr::theme_pubr()
      ) # + 
      # scale_y_continuous(# trans = log2_trans(),
      #                    labels =  trans_format("log2", math_format(2^.x)))
      , result = "plots"
    )
  # graphs
  
  variables <- levels(graphs$variables)
  plots <- graphs$plots %>% set_names(variables)
  
  dl = droplevels(stat.test.sig$variables)
  dlX = droplevels(stat.test.sig[intersect(grep('W', stat.test.sig$group1, invert = TRUE), 
                                           grep('W', stat.test.sig$group2, invert = TRUE)), ])$variables
  
  
  if(length(levels(dl)) > 0) {
  
    # cat(blue('sig\n\n'))
    bp <- ggdotplot(
      mydata.long[mydata.long$variables %in% levels(dl),], 
      x = "Stress", 
      y = "value", 
      fill = "Stress", #add = "mean_se",
      facet.by = c("variables"), 
      scales = "free", 
      palette = p.palette, #"npg", 
      legend = "none",
      add = c("jitter", 
              'mean_sd'),
      title = 'Sig diff',
      position = position_jitter(0.05)
    ) 
    bp = bp +
      stat_pvalue_manual(stat.test.sig, 
                         hide.ns = TRUE, 
                         tip.length = 0.02, 
                         step.increase = 0.05) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme(axis.text.x = element_text (angle = 45, 
                                        vjust = 1, 
                                        hjust = 1, 
                                        size = 7.5),
            axis.text.y = element_text (angle = 45, 
                                        vjust = 1, 
                                        hjust = 1, 
                                        size = 7.5))
    plot(bp)
  }
    
  # cat(blue('Without W\n\n'))
  # bp <- ggdotplot(
  #   mydata.long[mydata.long$variables %in% levels(dlX), ], 
  #   x = "Stress", 
  #   y = "value", 
  #   fill = "Stress", #add = "mean_se",
  #   facet.by = c("variables"), 
  #   scales = "free", 
  #   palette = p.palette, #"npg", 
  #   legend = "none",
  #   add = c("jitter", 
  #           'mean_sd',
  #           'mean_var'),
  #   title = 'Without W-specific differences',
  #   position = position_jitter(0.05)
  # ) 
  # bp +
  #   stat_pvalue_manual(stat.test.sig[stat.test.sig$variables %in% levels(dlX), ], 
  #                      hide.ns = TRUE, 
  #                      tip.length = 0.02, 
  #                      step.increase = 0.05) +
  #   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  #   theme(axis.text.x = element_text (angle = 45, 
  #                                     vjust = 1, 
  #                                     hjust=1, 
  #                                     size = 7.5))

  # plot(bp)
  
  if(length(levels(dlX)) > 0){
    
    cat(blue('Without W-only specific difference\n\n'))
    if(single){
      for(variable in levels(dlX)){
        stat.test.i <- dplyr::filter(stat.test.sig, variables == variable)
        stat.test.i$y.position = stat.test.i$y.position
        graph.i <- plots[[variable]] +
          labs(title = variable) +
          ggpubr::stat_pvalue_manual(stat.test.i, 
                                     label = "p.adj.signif",
                                     tip.length = 0.05, 
                                     step.increase = 0.2) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) #+
        # theme(axis.text.x = element_text (angle = 45, 
        #                                   vjust = 1, 
        #                                   hjust = 1, 
        #                                   size = 7.5),
        #       axis.text.y = element_text (angle = 45, 
        #                                   vjust = 1, 
        #                                   hjust = 1, 
        #                                   size = 7.5))
        print(graph.i)
      }  
    }
    
    
    ind = which(colnames(data) %in% levels(dlX))
    if(length(ind) > 1) {
      
      my.pairs.panels(df = data[, ind],
                      method = 'pearson',
                      palette = p.palette,
                      scale = FALSE,
                      main = 'Pairs of without W-specific differences',
                      cex.labels = p.cex.labels,
                      cex = p.cex)
    }
  }
  
  
  
  return(stat.test)
  
}


####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

my.logFC <- function(Ctrl, Ctrl.H, Ctrl.D, H, D, HD, W, tp, title){
  
  n = 11

  if(tp == 1){
    

    Ctrl.log = log2(Ctrl)
    H.log = log2(H)
    D.log = log2(D)
    HD.log = log2(HD)
    W.log = log2(W)
    
    Ctrl.log.mean = apply(Ctrl.log, 2, mean)
    H.log.mean = apply(H.log, 2, mean)
    HD.log.mean = apply(HD.log, 2, mean)
    D.log.mean = apply(D.log, 2, mean)
    W.log.mean = apply(W.log, 2, mean)
    
    H.logFC = H.log.mean - Ctrl.log.mean
    HD.logFC = HD.log.mean - Ctrl.log.mean
    D.logFC = D.log.mean - Ctrl.log.mean
    W.logFC = W.log.mean - Ctrl.log.mean
    
    log2FC = rbind(H.logFC, HD.logFC, D.logFC, W.logFC)
    # rownames(log2FC)
    # range(log2FC)
    
    

    m = max(abs(min(log2FC)), 
            abs(max(log2FC)))
    par(mfrow = c(2,2))
    pie(rep(1,n),
        col = rev(brewer.pal(n, 'RdBu')),
        labels = round(seq(-m, m, length.out = n), 2))
    par(mfrow = c(1,1))
    
    return(log2FC)
    
  } else if(tp == 2){
    
    
    keep1 = grep('Plant|Time|Stress', colnames(Ctrl.H), invert = TRUE)
    
    Ctrl.H.log = log2(Ctrl.H[, keep1])
    Ctrl.D.log = log2(Ctrl.D[, keep1])
    H.log = log2(H[,keep1])
    D.log = log2(D[,keep1])
    HD.log = log2(HD[,keep1])
    W.log = log2(W[,keep1])
    
    t1 = which(Ctrl.H$Time == 1)
    t2 = which(Ctrl.H$Time == 2)
    # https://9to5answer.com/using-apply-function-on-a-matrix-with-na-entries
    Ctrl.H.log.mean.1 = apply(Ctrl.H.log[t1, ], 2, mean, na.rm=TRUE)
    Ctrl.H.log.mean.2 = apply(Ctrl.H.log[t2, ], 2, mean, na.rm=TRUE)
    
    t1 = which(Ctrl.D$Time == 1)
    t2 = which(Ctrl.D$Time == 2)
    Ctrl.D.log.mean.1 = apply(Ctrl.D.log[t1, ], 2, mean, na.rm=TRUE)
    Ctrl.D.log.mean.2 = apply(Ctrl.D.log[t2, ], 2, mean, na.rm=TRUE)
    
    t1 = which(H$Time == 1)
    t2 = which(H$Time == 2)
    H.log.mean.1 = apply(H.log[t1, ], 2, mean, na.rm=TRUE)
    H.log.mean.2 = apply(H.log[t2, ], 2, mean, na.rm=TRUE)
    
    t1 = which(HD$Time == 1)
    t2 = which(HD$Time == 2)
    HD.log.mean.1 = apply(HD.log[t1, ], 2, mean, na.rm=TRUE)
    HD.log.mean.2 = apply(HD.log[t2, ], 2, mean, na.rm=TRUE)
    
    t1 = which(D$Time == 1)
    t2 = which(D$Time == 2)
    D.log.mean.1 = apply(D.log[t1, ], 2, mean, na.rm=TRUE)
    D.log.mean.2 = apply(D.log[t2, ], 2, mean, na.rm=TRUE)
    
    t1 = which(W$Time == 1)
    t2 = which(W$Time == 2)
    W.log.mean.1 = apply(W.log[t1, ], 2, mean, na.rm=TRUE)
    W.log.mean.2 = apply(W.log[t2, ], 2, mean, na.rm=TRUE)
    
    
    
    H.logFC.1 = H.log.mean.1 - Ctrl.H.log.mean.1
    HD.logFC.1 = HD.log.mean.1 - Ctrl.D.log.mean.1
    D.logFC.1 = D.log.mean.1 - Ctrl.D.log.mean.1
    W.logFC.1 = W.log.mean.1 - Ctrl.H.log.mean.1
    
    H.logFC.2 = H.log.mean.2 - Ctrl.H.log.mean.2
    HD.logFC.2 = HD.log.mean.2 - Ctrl.D.log.mean.2
    D.logFC.2 = D.log.mean.2 - Ctrl.D.log.mean.2
    W.logFC.2 = W.log.mean.2 - Ctrl.H.log.mean.2
    
    log2FC = rbind(H.logFC.1, H.logFC.2,
                   HD.logFC.1, HD.logFC.2,
                   D.logFC.1, D.logFC.2,
                   W.logFC.1, W.logFC.2)
    # rownames(log2FC)
    # range(log2FC, na.rm = TRUE)
    
    

    m = max(abs(min(log2FC, na.rm = TRUE)), abs(max(log2FC, na.rm = TRUE)))
    par(mfrow = c(2,2))
    pie(rep(1,n), col = rev(brewer.pal(n, 'RdBu')),
        labels = round(seq(-m, m, length.out = n), 2))
    par(mfrow = c(1,1))

    
    return(log2FC)

    
    
  } else {
    cat(red('ERROR'))
  }
  
}


####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####

my.NMDS <- function(mymat, stress.levels, time.levels, title, k, subset){
  
  set.seed(123456)
  
  if(k == 3){
    

    nmds = vegan::metaMDS(mymat, 
                          distance = "euclidean",
                          k = 3, try = 100, trymax = 1000,
                          trace = FALSE)
    
    data.scores = as.data.frame(vegan::scores(nmds)$sites)
    data.scores$Treatment = factor(subset$Stress, levels = stress.levels)
    data.scores$SamplingDay = factor(subset$Time, levels = time.levels)
    data.scores$PlantNo = as.factor(subset$Plant)
    
    g12 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
      geom_point(size = 3, aes( shape = Treatment, colour = SamplingDay )) +
      theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10),
            legend.text = element_text(size = 6, 
                                       #face ="bold", 
                                       colour ="black"
            ),
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 8, 
                                        colour = "black"#, 
                                        #face = "bold"
            ),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", 
                                                                            fill = NA, 
                                                                            size = 1.2),
            legend.key=element_blank()) +
      labs(x = "MDS1", colour = "SamplingDay", y = "MDS2", shape = "Treatment") + 
      ggtitle(title) +
      scale_color_manual(values = NMDS.palette)
    
    g13 = ggplot(data.scores, aes(x = NMDS1, y = NMDS3)) +
      geom_point(size = 3, aes( shape = Treatment, colour = SamplingDay )) +
      theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10),
            legend.text = element_text(size = 6, 
                                       #face ="bold", 
                                       colour ="black"
                                       ),
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 8, 
                                        colour = "black"#, 
                                        #face = "bold"
                                        ),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", 
                                                                            fill = NA, 
                                                                            size = 1.2),
            legend.key=element_blank()) +
      labs(x = "MDS1", colour = "SamplingDay", y = "MDS3", shape = "Treatment") + 
      ggtitle(title) +
      scale_color_manual(values = NMDS.palette)
    
    g23 = ggplot(data.scores, aes(x = NMDS2, y = NMDS3)) +
      geom_point(size = 3, aes( shape = Treatment, colour = SamplingDay )) +
      theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10),
            legend.text = element_text(size = 6, 
                                       #face ="bold", 
                                       colour ="black"
            ),
            legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 8, 
                                        colour = "black"#, 
                                        #face = "bold"
            ),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", 
                                                                            fill = NA, 
                                                                            size = 1.2),
            legend.key=element_blank()) +
      labs(x = "MDS2", colour = "SamplingDay", y = "MDS3", shape = "Treatment") + 
      ggtitle(title) +
      scale_color_manual(values = NMDS.palette)
    
    
    ds = data.scores
    ds$trt = paste(ds$Treatment, ds$SamplingDay, sep = '_')
    cent <- aggregate(cbind(NMDS1, NMDS2) ~ trt, data = ds, FUN = mean)
    cent$time = factor(as.numeric(gsub('.*_', '', cent$trt)), levels = time.levels)
    cent$stress = factor(gsub('_.*', '', cent$trt), levels = stress.levels)
    segs <- merge(ds, setNames(cent, c('trt','oNMDS1','oNMDS2')),
                  by = 'trt', sort = FALSE)
    ds$time = factor(as.numeric(gsub('.*_', '', ds$trt)), levels = time.levels)
    ds$stress = factor(gsub('_.*', '', ds$trt), levels = stress.levels)
    
    gc12 = ggplot(ds, aes(x = NMDS1, y = NMDS2, colour = ds$stress)) +
      geom_segment(data = segs,
                   mapping = aes(xend = oNMDS1, yend = oNMDS2)) + 
      geom_point(data = cent, 
                 size = 3, 
                 aes(colour = stress,
                 shape = time)) +
      geom_point(size = 1) +
      coord_fixed()     + 
      scale_color_manual(values = my.ggplot.palette) +
      ggtitle('Centroids') + 
      scale_shape_manual(name = "Day",
                         breaks = time.levels,
                         values =  c(1, 16, 0, 15, 2, 17, 8)) +
      theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10),
            legend.text = element_text(size = 6, 
                                       #face ="bold", 
                                       colour ="black"
            ),
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 8, 
                                        colour = "black"#, 
                                        #face = "bold"
            ),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", 
                                                                            fill = NA, 
                                                                            size = 1.2),
            legend.key=element_blank())

    
    ds = data.scores
    ds$trt = paste(ds$Treatment, ds$SamplingDay, sep = '_')
    cent <- aggregate(cbind(NMDS1, NMDS3) ~ trt, data = ds, FUN = mean)
    cent$time = factor(as.numeric(gsub('.*_', '', cent$trt)), levels = time.levels)
    cent$stress = factor(gsub('_.*', '', cent$trt), levels = stress.levels)
    segs <- merge(ds, setNames(cent, c('trt','oNMDS1','oNMDS3')),
                  by = 'trt', sort = FALSE)
    ds$time = factor(as.numeric(gsub('.*_', '', ds$trt)), levels = time.levels)
    ds$stress = factor(gsub('_.*', '', ds$trt), levels = stress.levels)
    
    gc13 = ggplot(ds, aes(x = NMDS1, y = NMDS3, colour = ds$stress)) +
      geom_segment(data = segs,
                   mapping = aes(xend = oNMDS1, yend = oNMDS3)) +
      geom_point(data = cent, 
                 size = 3, 
                 aes(colour = stress,
                     shape = time)) +
      geom_point(size = 1) +
      coord_fixed()     + 
      scale_color_manual(name = "Stress", values = my.ggplot.palette) +
      ggtitle('Centroids') + 
      scale_shape_manual(name = "Day",
                         breaks = time.levels,
                         values =  c(1, 16, 0, 15, 2, 17, 8)) +
      theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10),
            legend.text = element_text(size = 6, 
                                       #face ="bold", 
                                       colour ="black"
            ),
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 8, 
                                        colour = "black"#, 
                                        #face = "bold"
            ),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", 
                                                                            fill = NA, 
                                                                            size = 1.2),
            legend.key=element_blank())
    
    ds = data.scores
    ds$trt = paste(ds$Treatment, ds$SamplingDay, sep = '_')
    cent <- aggregate(cbind(NMDS2, NMDS3) ~ trt, data = ds, FUN = mean)
    cent$time = factor(as.numeric(gsub('.*_', '', cent$trt)), levels = time.levels)
    cent$stress = factor(gsub('_.*', '', cent$trt), levels = stress.levels)
    segs <- merge(ds, setNames(cent, c('trt','oNMDS2','oNMDS3')),
                  by = 'trt', sort = FALSE)
    ds$time = factor(as.numeric(gsub('.*_', '', ds$trt)), levels = time.levels)
    ds$stress = factor(gsub('_.*', '', ds$trt), levels = stress.levels)
    
    gc23 = ggplot(ds, aes(x = NMDS2, y = NMDS3, colour = ds$stress)) +
      geom_segment(data = segs,
                   mapping = aes(xend = oNMDS2, yend = oNMDS3)) +
      geom_point(data = cent, 
                 size = 3, 
                 aes(colour = stress,
                     shape = time)) +
      geom_point(size = 1) +
      coord_fixed()     + 
      scale_color_manual(values = my.ggplot.palette) +
      ggtitle('Centroids') + 
      scale_shape_manual(name = "Day",
                         breaks = time.levels,
                         values =  c(1, 16, 0, 15, 2, 17, 8)) +
      theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10),
            legend.text = element_text(size = 6, 
                                       #face ="bold", 
                                       colour ="black"
            ),
            legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 8, 
                                        colour = "black"#, 
                                        #face = "bold"
            ),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", 
                                                                            fill = NA, 
                                                                            size = 1.2),
            legend.key=element_blank())
    
    
    return(list(g12, g13, g23, gc12, gc13, gc23))
  
  } else {
    
    nmds = vegan::metaMDS(mymat, distance = "euclidean",
                          k = 2, try = 100, trymax = 1000,
                          trace = FALSE)
    
    data.scores = as.data.frame(vegan::scores(nmds)$sites)
    data.scores$Treatment = factor(subset$Stress, levels = stress.levels)
    data.scores$SamplingDay = factor(subset$Time, levels = time.levels)
    data.scores$PlantNo = as.factor(subset$Plant)
    
    g12 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2, label =  SamplingDay)) +
      geom_point(size = 4, aes( shape = SamplingDay, colour = Treatment)) +
      theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
            axis.text.x = element_text(colour = "black", face = "bold", size = 12),
            legend.text = element_text(size = 12, face ="bold", colour ="black"),
            legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
            legend.title = element_text(size = 14, colour = "black", face = "bold"),
            panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
            legend.key=element_blank()) +
      labs(x = "MDS1", colour = "Treatment", y = "MDS2", shape = "SamplingDay") + 
      ggtitle('Hormonomics') +
      geom_text() +
      scale_color_manual(values = my.ggplot.palette) +
      scale_shape_manual(values=c(19, 19, 15, 15, 18, 18, 18))# c(1, 2, 6, 5, 8))
    
    
    ds = data.scores
    ds$trt = paste(ds$Treatment, ds$SamplingDay, sep = '_')
    cent <- aggregate(cbind(NMDS1, NMDS2) ~ trt, data = ds, FUN = mean)
    cent$t = ordered(as.numeric(gsub('.*_', '', cent$trt)), levels = time.levels)
    segs <- merge(ds, setNames(cent, c('trt','oNMDS1','oNMDS2')),
                    by = 'trt', sort = FALSE)
    
    ggplot(ds, aes(x = NMDS1, y = NMDS2, colour = trt)) +
      geom_segment(data = segs,
                   mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
      geom_point(data = cent, size = 5, pch = cent$t) +# centroids
      geom_point(pch = ordered(as.numeric(gsub('.*_', '', ds$SamplingDay)), levels = time.levels)) + # sample scores
      coord_fixed() + # same axis scaling
      scale_color_manual(values = rep(my.ggplot.palette[c(1, 3, 2, 4, 5)], c(7, 5, 7, 5, 7))) +
      ggtitle('Hormonomics')

    
    return(list(g12))
  }
}





