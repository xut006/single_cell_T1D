##### Based on the number of genes expressed plot a histogram for each cell ####
##### with x-axis of number of genes expressed and y-axis of number of cells ####
geneExpressHistogram <- function(ctSelect, num_top = 20){
  
  ## append a histogram column to store number of gene expressed for each cell
  ctSelect$histogram <- NA
  ## loop through the data set to get num_express
  for (i in 1:nrow(ctSelect)){
    num_express = 0
    for (j in 8:(ncol(ctSelect)-1)){
      if (ctSelect[i,j] != 0){
        num_express = num_express+1
      }
    }
    ctSelect[i, "histogram"] = num_express
  }
  
  
  ## use the histogram column to plot histogram
  par(mar=c(50, 50, 50, 50))
  hist <- ggplot(ctSelect, aes(ctSelect$histogram)) + 
    geom_histogram(colour = "wheat2", fill = "grey50", binwidth = 1) + 
    scale_x_continuous(breaks = seq(0, 96, 5), minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(0, 100, 5), minor_breaks = NULL) +
    ggtitle("Gene Expression Histogram") +
    xlab("Number of Gene Expressed") +
    ylab("Number of Cells") +
    theme(text=element_text(size=25),
          axis.title.x=element_text(vjust=-0.5),
          plot.margin=unit(c(15,2,15,2),"cm"))
  print(hist)
  par(mar=c(5, 5, 5, 5))
  
  print(paste("Total number of cells: ", nrow(ctSelect),sep=""))

  
  ## Find the most active genes
  active_genes <- data.frame()
  for (j in 8:(ncol(ctSelect)-1)){
    num_active = 0
    for (i in 1:(nrow(ctSelect)-1)){
      if (ctSelect[i,j] != 0){
        num_active  = num_active+1
      }
    }
    new_gene <- data.frame(c(num_active))
    rownames(new_gene) <- c(colnames(ctSelect)[j])
    active_genes <- rbind(active_genes, new_gene)
  }
  colnames(active_genes) <- "num_active"
  
  ## get the top expressed genes
  top_20 <- row.names(active_genes)[order(active_genes$num_active, decreasing = TRUE)[1:num_top]]
  print(paste("Top ", num_top, " expressing genes: ",sep=""))
  print(top_20)
  
  return(top_20)
}



