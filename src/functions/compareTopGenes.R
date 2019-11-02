compareTopGenes <- function(ctSelect, topGenes){
  ctSelect$crossCompare <- NA
  
  ## loop through the data set to assign
  # for (i in 1:nrow(ctSelect)){
  #   num_express = 0
  #   for (j in topGenes){
  #     if (ctSelect[i,j] != 0){
  #       num_express = num_express+1
  #     }
  #   }
  #   ctSelect[i, "crossCompare"] = num_express
  # }
  
  
  ## Transform datafram for heatmap plot
  heatmap_df <- ctSelect[topGenes]%>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
  names(heatmap_df) <- c("cells", "Top_Genes", "expression")
  
  ## plot Heatmap
  plHeatmap <- ggplot(heatmap_df, aes(x = cells, y = Top_Genes, fill = expression)) +
    geom_tile(colour = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    ggtitle("Heatmap for expresseion level of top expressed genes") +
    theme_minimal() +
    theme(text=element_text(size=28),
      plot.margin=unit(c(14,9,14,9),"cm"))
  print(plHeatmap)
  
  print(paste("Total number of cells: ", nrow(ctSelect),sep=""))
  
}