compareTopGenes <- function(ctSelect, topGenes){
  ctSelect$crossCompare <- NA
  
  ## loop through the data set to assign
  for (i in 1:nrow(ctSelect)){
    num_express = 0
    for (j in topGenes){
      if (ctSelect[i,j] != 0){
        num_express = num_express+1
      }
    }
    ctSelect[i, "crossCompare"] = num_express
  }
  ctRep <- ctSelect
  
  
  ## remove duplicated rows
  if(anyDuplicated(ctRep[,8:ncol(ctRep)]) > 0){
    ctRep <- ctRep[-which(duplicated(ctRep[,8:ncol(ctRep)])),]
  }
  
  ## remove genes without variance
  vars <- NULL
  for(i in 8:ncol(ctRep)){
    vars <- c(vars, var(ctRep[,i], na.rm=T))
  }
  ctGenesNoVar <- ctRep[which(vars == 0 | is.na(vars))+7]
  ctRep <- ctRep[,c(1:7, (which(!is.na(vars) & vars!=0)+7))]
  
  ## run t-SNE
  set.seed(1)
  tsne_out <- Rtsne(as.matrix(ctRep[, 8:ncol(ctRep)]), perplexity = 30)
  
  tsne_y <- as.data.frame(cbind(tsne_out$Y, 
                                ctRep$cellSource,
                                ctRep$crossCompare,
                                ctRep$age,
                                ctRep$probe,
                                ctRep[, 8:(ncol(ctRep)-1)]))
  
  names(tsne_y)[1:6] <- c("y1", "y2", "cohort", "crossCompare", "marker", "patient")
  
  for(i in c(1,2,7:ncol(tsne_y))){
    tsne_y[, i] <- as.numeric(tsne_y[, i])
  }
  
  
  ## set color for t-SNE plot
  pointSize <- 4
  # myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  # clusterPalette <- "Set3"
  shapeVals <- c(19, 17, 15, 18)
  numcohorts <- length(unique(tsne_y$cohort))
  shapeVals <- shapeVals[1:numcohorts]
  
  
  ## label cells by corssCompare
  plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
    geom_point(aes(color=tsne_y$crossCompare), size = pointSize, alpha = 1) +
    scale_color_viridis(option = "D") +
    #scale_colour_gradient(low = "white", high = "black") +
    scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                       minor_breaks = NULL) +
    scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                       minor_breaks = NULL) +
    scale_shape_manual(values = shapeVals) +
    guides(color=guide_legend(title="number of top genes expressed", order = 1)) +
    ggtitle("t-SNE colored by top expressed genes") +
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
  print(paste("Total number of cells: ", nrow(ctSelect),sep=""))
}