plotTSNE <- function(ctClust, colorby = c("kmeans.cluster", "patient", "cohort", "marker", "markerType")){
  ctRep <- ctClust
  
  ## remove duplicated rows
  if(anyDuplicated(ctRep[,10:ncol(ctRep)]) > 0){
    ctRep <- ctRep[-which(duplicated(ctRep[,10:ncol(ctRep)])),]
  }
  
  ## remove genes without variance
  vars <- NULL
  for(i in 10:ncol(ctRep)){
    vars <- c(vars, var(ctRep[,i], na.rm=T))
  }
  ctGenesNoVar <- ctRep[which(vars == 0 | is.na(vars))+9]
  ctRep <- ctRep[,c(1:9, (which(!is.na(vars) & vars!=0)+9))]
  
  ## run t-SNE
  set.seed(1)
  tsne_out <- Rtsne(as.matrix(ctRep[, 10:ncol(ctRep)]), perplexity = 30)
  
  tsne_y <- as.data.frame(cbind(tsne_out$Y, 
                                ctRep$cohort,
                                ctRep$kmeans.cluster,
                                ctRep$marker,
                                ctRep$patients,
                                ctRep$markerType,
                                ctRep[, 10:ncol(ctRep)]))
  
  names(tsne_y)[1:7] <- c("y1", "y2", "cohort", "kmeans.cluster", "marker", "patient", "markerType")
  #tsne_y$kmeans.cluster <- factor(tsne_y$kmeans.cluster, levels=clusters)
  #relcohorts <- unique(tsne_y$cohort)
  #relcohorts <- relcohorts[match(c("I", "LN", "S", "PB"), relcohorts, nomatch=F)]
  #tsne_y$cohort <- factor(tsne_y$cohort, levels=relcohorts)
  
  for(i in c(1,2,8:ncol(tsne_y))){
    tsne_y[, i] <- as.numeric(tsne_y[, i])
  }
  
  
  
  ## set color for t-SNE plot
  pointSize <- 4
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  clusterPalette <- "Set3"
  shapeVals <- c(19, 17, 15, 18)
  numcohorts <- length(unique(tsne_y$cohort))
  shapeVals <- shapeVals[1:numcohorts]
  
  
  ## label cells by kmeans.cluster
  if("kmeans.cluster" %in% colorby){
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$kmeans.cluster)), size = pointSize, alpha = 1) +
      scale_colour_brewer(palette = clusterPalette) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="cluster", order = 1)) +
      ggtitle("t-SNE colored by kmeans.cluster") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        panel.grid.minor=element_line(color="gray90"),
        panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"))
    print(plTSNE)
  }
  
  
  ## label cells by patient
  if("patient" %in% colorby){
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$patient)), size = pointSize, alpha = 1) +
      scale_color_manual(values = c("gray19","#DD8D29","#E1C408","#84BB78","#859C78","#DB6E07","#B40F20","#798E87",
                                    "#C27D38","#CCC591","#29211F","#F3DF6C","#CEAB07","#D5D5D3","#24281A",
                                    "violetred4", "chartreuse3")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="patient", order = 1)) +
      ggtitle("t-SNE colored by different patients") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        panel.grid.minor=element_line(color="gray90"),
        panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"))
    print(plTSNE)
  }
  
  
  ## label cells by cohort
  if("cohort" %in% colorby){
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$cohort)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("#21908CFF","#FDE725FF","gray48","#440154FF")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="cohort", order = 1)) +
      ggtitle("t-SNE colored by patient cohorts") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        panel.grid.minor=element_line(color="gray90"),
        panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"))
    print(plTSNE)
  }
  
  
  ## label cells by type of markers
  if("markerType" %in% colorby){
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$markerType)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("#9C964A", "#85D4E3", "#F4B5BD")) +
      #scale_colour_manual(values=wes_palette("Moonrise3", 3)) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="markers", order = 1)) +
      ggtitle("t-SNE colored by different types of markers") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        panel.grid.minor=element_line(color="gray90"),
        panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"))
    print(plTSNE)
  }
  
  
  ## label cells by marker
  if("marker" %in% colorby){
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$marker)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("#9A8822","#F5CDB4","#F8AFA8","#FDDDA0","#74A089","#FF0000","#556A5B","#50A45C",
                                   "#F2AD00","#F69100","#C49647","#5BBCD6", "royalblue3","moccasin","#9986A5")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="markers", order = 1)) +
      ggtitle("t-SNE colored by different markers") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        panel.grid.minor=element_line(color="gray90"),
        panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"))
    print(plTSNE)
  }
  
  
}