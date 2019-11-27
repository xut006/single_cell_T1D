SpecificHeatMap <- function(ctSelect, low_Cut, high_cut, num_top = 20, all_PD1ICOS, top_20PD1ICOS){  
  ## append a explevel column to store number of gene expressed for each cell
  ctSelect$explevel <- NA
  ## loop through the data set to get num_express
  for (i in 1:nrow(ctSelect)){
    num_express = 0
    for (j in 8:(ncol(ctSelect)-1)){
      if (ctSelect[i,j] != 0){
        num_express = num_express+1
      }
    }
    ctSelect[i, "explevel"] = num_express
  }
  
  low_expression <- subset(ctSelect, explevel < high_cut & explevel >= low_Cut) %>%
                    subset(select = -c(explevel))
  high_expression <- subset(ctSelect, explevel >= high_cut) %>%
                    subset(select = -c(explevel))
  
  
  ## For cells with low expression level
  low_top_20 <- geneExpressHistogram(low_expression, num_top)
  topT_T_df <- low_expression[low_top_20]
  pheatmap(t(topT_T_df), fontsize = 28)
  print("Specific analysis for cells with low expression level")
  # Cross Comparision with PD1/ICOS genes
  topPI_T_df <- low_expression[top_20PD1ICOS]
  pheatmap(t(topPI_T_df), fontsize = 28)
  print("Cross Comparision: Low Expression Tetramer cells with Top 20 PD1/ICOS genes")
  ## Cross Comparision with low express tetramer genes
  topT_PI_df <- all_PD1ICOS[low_top_20]
  pheatmap(t(topT_PI_df), fontsize = 28)
  print("Cross Comparision: PD1/ICOS cells with Top 20 Low Expression Tetramer genes")
  
  
  ## For cells with high expression level
  high_top_20 <- geneExpressHistogram(high_expression, num_top)
  topT_T_df <- high_expression[high_top_20]
  pheatmap(t(topT_T_df), fontsize = 28)
  print("Specific analysis for cells with high expression level")
  # Cross Comparision with PD1/ICOS
  topPI_T_df <- high_expression[top_20PD1ICOS]
  pheatmap(t(topPI_T_df), fontsize = 28)
  print("Cross Comparision: High Expression Tetramer cells with Top 20 PD1/ICOS genes")
  ## Cross Comparision with high express tetramer genes
  topT_PI_df <- all_PD1ICOS[high_top_20]
  pheatmap(t(topT_PI_df), fontsize = 28)
  print("Cross Comparision: PD1/ICOS cells with Top 20 high Expression Tetramer genes")
  
}
  
  