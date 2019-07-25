plotTSNE <- function(ctClust, colorby = c("kmeans.cluster", "probe", "patient", "age")){
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
                                ctRep$cellSource,
                                ctRep$kmeans.cluster,
                                ctRep$age,
                                ctRep$probe,
                                ctRep[, 10:ncol(ctRep)]))
  
  names(tsne_y)[1:6] <- c("y1", "y2", "tissue", "kmeans.cluster", "age", "probe")
  #tsne_y$kmeans.cluster <- factor(tsne_y$kmeans.cluster, levels=clusters)
  #relTissues <- unique(tsne_y$tissue)
  #relTissues <- relTissues[match(c("I", "LN", "S", "PB"), relTissues, nomatch=F)]
  #tsne_y$tissue <- factor(tsne_y$tissue, levels=relTissues)
  
  for(i in c(1,2,7:ncol(tsne_y))){
    tsne_y[, i] <- as.numeric(tsne_y[, i])
  }
  
  
  
  ## set color for t-SNE plot
  pointSize <- 4
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  clusterPalette <- "Set3"
  shapeVals <- c(19, 17, 15, 18)
  numTissues <- length(unique(tsne_y$tissue))
  shapeVals <- shapeVals[1:numTissues]
  
  
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
      ggtitle("t-SNE between tissues (colored by kmeans.cluster)") +
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
  
  
  ## label cells by probe
  if("probe" %in% colorby){
    tsne_y$probeColor <- NA
    tsne_y[which(tsne_y$probe=="UCM5"), "probeColor"] <- "#A6CEE3"
    tsne_y[which(tsne_y$probe=="UCM6"), "probeColor"] <- "#1F78B4"
    tsne_y[which(tsne_y$probe=="UCM9"), "probeColor"] <- "blue3"
    tsne_y[which(tsne_y$probe=="UCM10"), "probeColor"] <- "turquoise"
    tsne_y[which(tsne_y$probe=="UCM12"), "probeColor"] <- "navy"
    tsne_y[which(tsne_y$probe=="UCM13"), "probeColor"] <- "maroon"
    tsne_y[which(tsne_y$probe=="UCM14"), "probeColor"] <- "darkmagenta"
    tsne_y[which(tsne_y$probe=="UCM15"), "probeColor"] <- "thistle3"
    tsne_y[which(tsne_y$probe=="UCM16"), "probeColor"] <- "peru"
    tsne_y[which(tsne_y$probe=="UCM17"), "probeColor"] <- "purple4"
    tsne_y[which(tsne_y$probe=="NBD1"), "probeColor"] <- "#B2DF8A"
    tsne_y[which(tsne_y$probe=="NBD3"), "probeColor"] <- "#33A02C"
    tsne_y[which(tsne_y$probe=="NBD4"), "probeColor"] <- "yellow2"
    #Probes <- tsne_y$probeColor
    #tsne_y <- subset(tsne_y, select=-c(probeColor))
    #AnnoColors <- cbind(AnnoColors, Probes)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE between tissues (colored by probe)", col=tsne_y$probeColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$probe)), size = pointSize, alpha = 1) +
      scale_color_manual(values=c("#A6CEE3", "#1F78B4", "blue3", "turquoise", "navy", "maroon", "darkmagenta", "darkmagenta", "thistle3", "peru", "purple4", "#B2DF8A", "#33A02C", "yellow2")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="probe", order = 1)) +
      ggtitle("t-SNE between tissues (colored by probe)") +
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
  
  
  ## label cells by tissue
  if("patient" %in% colorby){
    tsne_y$patientColor <- NA
    tsne_y$patient <- NA
    tsne_y[grep("UCM", tsne_y$probe, fixed=TRUE), "patientColor"] <- "turquoise3"
    tsne_y[grep("NBD", tsne_y$probe, fixed=TRUE), "patientColor"] <- "gold"
    tsne_y[grep("UCM", tsne_y$probe, fixed=TRUE), "patient"] <- "UCM"
    tsne_y[grep("NBD", tsne_y$probe, fixed=TRUE), "patient"] <- "NBD"
    #Tissues <- tsne_y$patientColor
    
    #tsne_y <- subset(tsne_y, select=-c(patientColor))
    #AnnoColors <- cbind(AnnoColors, Tissues)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE colored by patient", col=tsne_y$patientColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$patient)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("turquoise3", "gold")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="patient", order = 1)) +
      ggtitle("t-SNE colored by patient") +
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
  
  
  ## label cells by age
  if("age" %in% colorby){
    tsne_y$ageColor <- NA
    tsne_y[which(tsne_y$age=="6"), "ageColor"] <- "deepskyblue2"
    tsne_y[which(tsne_y$age=="12"), "ageColor"] <- "peachpuff2"
    tsne_y[which(tsne_y$age=="infl"), "ageColor"] <- "deepskyblue2"
    tsne_y[which(tsne_y$age=="uninfl"), "ageColor"] <- "navy"
    tsne_y[which(tsne_y$age=="blood"), "ageColor"] <- "orangered"
    #Ages <- tsne_y$ageColor
    #tsne_y <- subset(tsne_y, select=-c(ageColor))
    #AnnoColors <- cbind(AnnoColors, Ages)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE colored by inflamed site", col=tsne_y$ageColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$age)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("orangered", "navy", "deepskyblue2")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="inflamation", order = 1)) +
      ggtitle("t-SNE colored by inflamed site") +
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