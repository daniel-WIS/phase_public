# Plot Functions ----------------------------------------------------------
library(scales) 
library(ggplot2)
library(cowplot)

#Plot without foci, takes dataframe with *1* multimer as argument
plotphase.nofoci = function (xmer.data, y_label = "[Tetramer or Hexamer] (nM)",plot_color="red") { ###TO-DO: y_label SHOULD BE INFERRED
  ###Check whether data supplied, if not plot empty coordinate system
  if (missing(xmer.data)) {
    
    dat2add <- geom_blank()
    ylab2add <- ylab(y_label)
    annotation2add <- geom_blank()
    title2add <- geom_blank()
    hline2add <- geom_blank()
    vline2add <- geom_blank()
    
    
  } else {
    dat2add <- geom_point(
      data = xmer.data,
      aes(RFPconc, GFPconc),
      #alpha = 0.5,
      size = 0.5,
      color = plot_color,
      shape = 16
      
    )
    
    diagonal2add <-  geom_segment(aes(
      x = 0,
      xend = 1e4,
      y = 0,
      yend = 1e4
    ),
    color = "darkgrey",
    alpha = .7)
    
    ylab2add <- ylab(y_label)
    
    title2add <-
      ggtitle(paste(as.character(unique(xmer.data$xmer))))
    
    annotation2add <- annotate("text",x=2e4, y=1.5, label=paste(as.character(unique(xmer.data$xmer)),"n =",as.character(nrow(xmer.data))),size=3.5)
    hline2add <- geom_hline(yintercept = 3, linetype="dotted", col="darkgrey")
    vline2add <- geom_vline(xintercept = 3, linetype="dotted", col="darkgrey")
  }
  
  ####Start ggplot build
  ggplot() +
    
    dat2add +
    
    #diagonal2add +
    
    #title2add +
    
    #xlab("[Dimer] (nM)") +
    
    #ylab2add +
    
    annotation2add +
    
    geom_segment(aes(
      x = 0,
      xend = 1e4,
      y = 0,
      yend = 1e4
    ),
    color = "darkgrey",
    alpha = .7) +
    
    hline2add +
    vline2add +
    
    scale_y_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000)
      #label = trans_format(log10, math_format())
    ) +
    
    scale_x_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000)
      #label = trans_format(log10, math_format())
    ) +
    
    theme(
      axis.line = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) 
}

#Plot without foci, takes dataframe with *1* multimer as argument
plotphase.overlay2 = function (bottom_layer,top_layer, y_label = "[Tetramer or Hexamer] (nM)") { ###TO-DO: y_label SHOULD BE INFERRED
  ###Check whether data supplied, if not plot empty coordinate system
  
  ggplot() + geom_point(
    data = bottom_layer,
    aes(RFPconc, GFPconc),
    #alpha = 0.5,
    size = 0.5,
    color = "blue",
    shape = 16
  ) +
    
    geom_point(
      data = top_layer,
      aes(RFPconc, GFPconc),
      #alpha = 0.5,
      size = 0.5,
      color = "orange",
      shape = 16
      
    ) +
    
    geom_segment(aes(
      x = 0,
      xend = 1e4,
      y = 0,
      yend = 1e4
    ),
    color = "darkgrey",
    alpha = .7) +
    
    ylab(y_label) +
    annotate("text",x=2e4, y=1.5, label=paste(as.character(unique(top_layer$xmer)),"n =",as.character(nrow(top_layer))),size=3.5) +
    geom_hline(yintercept = 3, linetype="dotted", col="darkgrey") +
    geom_vline(xintercept = 3, linetype="dotted", col="darkgrey") +
    xlab("[Dimer] (nM)") +
    geom_segment(aes(
      x = 0,
      xend = 1e4,
      y = 0,
      yend = 1e4
    ),
    color = "darkgrey",
    alpha = .7) +
    
    scale_y_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      label = trans_format(log10, math_format())
    ) +
    
    scale_x_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      label = trans_format(log10, math_format())
    ) +
    
    theme(
      axis.line = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
      
    ) 
}


#Plot without foci takes multimer as argument
plotphase.dan = function (g) {
  ggplot () + geom_point(
    data = subset(tabdan, xmer == g  &
                    GFP_foci == 0 &
                    RFP_foci == 0),
    aes(RFPconc, GFPconc),
    alpha = 0.5,
    size = 0.7
  ) +
    
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) +
    ggtitle(paste(as.character(g))) +
    scale_y_continuous(
      trans = log2_trans(),
      limits = c(1, 50000),
      breaks = c(0, 1, 10, 100, 1000, 10000)
    ) +
    scale_x_continuous(
      trans = log2_trans(),
      limits = c(1, 50000),
      breaks = c(0, 1, 10, 100, 1000, 10000)
    ) +
    geom_abline(slope = 1,
                intercept = (-0):(0),
                color = "red")
  
}

#Plot with foci  
plotphase.dan.wfoci = function (g) {
  ggplot () + geom_point(
    data = subset(tabdan, xmer == g  &
                    GFP_foci == 0 &
                    RFP_foci == 0),
    aes(RFPconc, GFPconc),
    alpha = 0.5,
    size = 0.7
  ) +
    geom_point(
      data = subset(tabdan, xmer == g  &
                      GFP_foci != 0 &
                      RFP_foci != 0),
      aes(RFPconc, GFPconc),
      alpha = 0.5,
      color = "red"
    ) +
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) +
    ggtitle(paste(as.character(g))) +
    scale_y_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      label = trans_format(log10,math_format())
    ) +
    scale_x_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      label = trans_format(log10,math_format())
    ) +
    geom_abline(slope = 1,
                intercept = (-0):(0),
                color = "red")
  
}

#Overlay Plot
plotphase.dan.overlap = function (bottom_layer, top_layer) {
  ggplot () +
    
    geom_point(
      data = subset(bottom_layer,  GFP_foci == 0 & RFP_foci == 0),
      aes(RFPconc, GFPconc), #aes(RFPconc, GFPconc, colour=Tetramer) to make legend
      alpha = .5,
      size = .5, 
      color = "#ff8902", #orange
      shape = 16
    ) +
    
    geom_point(
      data = subset(top_layer, GFP_foci == 0 & RFP_foci == 0),
      aes(RFPconc, GFPconc), #aes(RFPconc, GFPconc, colour=Hexamer) to make legend
      alpha = .5,
      size = .5, 
      color = "#1302ff", # blue
      shape = 16
    ) +
    
    geom_segment(
      aes(
        x = 0,
        xend = 1e4,
        y = 0,
        yend = 1e4
      ),
      color = "darkgrey",
      alpha = .7
    ) +
    
    
    scale_x_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      label = trans_format(log10,math_format())
    ) +
    
    scale_y_continuous(
      trans = log2_trans(),
      limits = c(1, 1e5),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      label = trans_format(log10,math_format())
    ) +
    
    
    
    theme(
      axis.line = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #legend.position = c(.9, .9),
      #legend.key = element_rect(fill = NA, color = NA)
      
    ) +
    
    ggtitle(as.character(unique(top_layer$xmer))) +
    
    #guides(colour = guide_legend(override.aes = list(size = 3))) +
    
    #scale_colour_manual(
    #  name = '',
    #  values = c('Tetramer' = 'royalblue2', 'Hexamer' = 'orange')
    #) +
    
    xlab("[Dimer] (nM)") +
    
    ylab("[Tetramer or Hexamer] (nM)")+
    
    annotate("text",x=2e4, y=1.5, label=paste("n =",as.character(nrow(top_layer))),size=3.5)
  
}

plotphase.dan.overlap(bottom_layer = dat.p53_reduced, top_layer = dat.601_reduced)
