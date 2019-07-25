#### load UC-MAIT data ####


dataLoadUCMAIT <- function(){
  ## set data directory
  dataDir <- paste(baseDir, "data/", sep="")
  
  ### read in data
  ## probes: NBD1 - UCM5
  ctTable1 <- read.csv(paste(dataDir, "1362292369-NBD1_UCM5-HeatMapResultsUpdated.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable1$cellSource <- "blood"
  ctTable1$probe <- NA
  for(i in 1:nrow(ctTable1)) {
    if (grepl("UCM5", ctTable1[i,2], fixed=TRUE)){
      ctTable1[i, "probe"] <- "UCM5"
    }
    if (grepl("NBD1", ctTable1[i,2], fixed=TRUE)){
      ctTable1[i, "probe"] <- "NBD1"
    }
  }
  ctTable1$age <- "blood"
  
  ### read in data
  ## probes: NBD3 - UCM6
  ctTable2 <- read.csv(paste(dataDir, "1362292371-NBD3-UCM6-HeatMapResults.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable2$cellSource <- "blood"
  ctTable2$probe <- NA
  for(i in 1:nrow(ctTable2)) {
    if (grepl("UCM6", ctTable2[i,2], fixed=TRUE)){
      ctTable2[i, "probe"] <- "UCM6"
    }
    if (grepl("NBD3", ctTable2[i,2], fixed=TRUE)){
      ctTable2[i, "probe"] <- "NBD3"
    }
  }
  ctTable2$age <- "blood"
  
  ### read in data
  ## probes: NBD4
  ctTable3 <- read.csv(paste(dataDir, "1362351425-NBD4.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable3$cellSource <- NA
  ctTable3$probe <- "NBD4"
  ctTable3$age <- NA
  
  ### read in data
  ## probes: UCM14
  ctTable4 <- read.csv(paste(dataDir, "1362351426-UCM14.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable4$cellSource <- NA
  ctTable4$probe <- "UCM14"
  ctTable4$age <- NA
  
  ### read in data
  ## probes: UCM13
  ctTable5 <- read.csv(paste(dataDir, "1362351427-UCM13.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable5$cellSource <- NA
  ctTable5$probe <- "UCM13"
  ctTable5$age <- NA
  
  ### read in data
  ## probes: UCM12
  ctTable6 <- read.csv(paste(dataDir, "1362351428_UCM12.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable6$cellSource <- NA
  ctTable6$probe <- "UCM12"
  ctTable6$age <- NA
  
  ### read in data
  ## probes: UCM10
  ctTable7 <- read.csv(paste(dataDir, "1362351431_UCM10.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable7$cellSource <- NA
  ctTable7$probe <- "UCM10"
  ctTable7$age <- NA
  
  ### read in data
  ## probes: UCM15
  ctTable8 <- read.csv(paste(dataDir, "1362356328-UCM15.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable8$cellSource <- NA
  ctTable8$probe <- "UCM15"
  ctTable8$age <- NA
  
  ### read in data
  ## probes: UCM16
  ctTable9 <- read.csv(paste(dataDir, "1362356329-UCM16.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable9$cellSource <- NA
  ctTable9$probe <- "UCM16"
  ctTable9$age <- NA
  
  ### read in data
  ## probes: UCM17
  ctTable10 <- read.csv(paste(dataDir, "1362351435_Processed-TNET001-UCM17.csv", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable10 <- subset(ctTable10, grepl("UCM17", ctTable10$Name, fixed = TRUE))
  ctTable10$cellSource <- NA
  ctTable10$probe <- "UCM17"
  ctTable10$age <- NA
  
  
  ## everything
  ctTableCombine <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, ctTable9)
  ctTableCombine <- rbind(ctTableCombine, ctTable10)
  
  ## change column names
  names(ctTableCombine)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## rename cellSource column value
  ctTableCombine[grepl("Infl", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"
  ctTableCombine[grepl("Uninf", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"
  ctTableCombine[grepl("Blood", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "blood"
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  
  ## create age (inflamed/uninflamed) column
  ctTableCombine[grepl("Infl", ctTableCombine$cellType, fixed=TRUE), "age"] <- "infl"
  ctTableCombine[grepl("Uninf", ctTableCombine$cellType, fixed=TRUE), "age"] <- "uninfl"
  ctTableCombine[grepl("Blood", ctTableCombine$cellType, fixed=TRUE), "age"] <- "blood"
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "age"] <- "controlIslet"
  
  
  return(ctTableCombine)
}







