#### load UC-MAIT data ####


dataLoadT1D <- function(){
  ## set data directory
  baseDir <- "~/Desktop/Su Lab/single_cell_T1D/"
  dataDir <- paste(baseDir, "data/", sep="")
  
  ### read in Children data
  ctTable1 <- read.csv(paste(dataDir, "1362015470-RAD001-FrozenvsFresh.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable1$probe <- "RAD001"
  
  ctTable2 <- read.csv(paste(dataDir, "RAD003.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable2$probe <- "RAD003"
  
  ctTable3 <- read.csv(paste(dataDir, "1362351438-RAD006.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable3$probe <- "RAD006" 
  
  ctTable4 <- read.csv(paste(dataDir, "1362351437-Processed_RAD005.csv.trimmed", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE)
  ctTable4$probe <- "RAD005"
  
  ctTable5 <- read.csv(paste(dataDir, "1362356331-RAD004-TEY013.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable5 <- subset(ctTable5, grepl("RAD004", ctTable5$Name, fixed = TRUE))
  ctTable5$probe <- "RAD004"
  
  ctTable6 <- read.csv(paste(dataDir, "1362356332-RAD003_Insulin_ICOS_PD1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable6$probe <- "RAD003"
  
  ctTable7 <- read.csv(paste(dataDir, "RAD001andRAD002.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable7$probe <- "RAD001"
  
  ctTable8 <- read.csv(paste(dataDir, "RAD003(p1220p1321).csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable8$probe <- "RAD003"
  
  
  ### read in Adults data
  ctTable9 <- read.csv(paste(dataDir, "1362292377-TEY006-TEY014.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable9$probe <- "TEY006"
  
  ctTable10 <- read.csv(paste(dataDir, "1362356331-RAD004-TEY013.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable10 <- subset(ctTable10, grepl("TEY013", ctTable10$Name, fixed = TRUE))
  ctTable10$probe <- "TEY013"
  
  
  ### read in At_Ristk data
  ctTable11 <- read.csv(paste(dataDir, "1362351434-processedresults-TNET002.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable11$probe <- "TNET002"
  
  ctTable12 <- read.csv(paste(dataDir, "1362351435_Processed-TNET001-UCM17.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable12 <- subset(ctTable12, grepl("TNET001", ctTable12$Name, fixed = TRUE))
  ctTable12$probe <- "TNET001"
  
  ctTable13 <- read.csv(paste(dataDir, "1362351436_Processed-TNET003.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable13$probe <- "TNET003"
  
  ctTable14 <- read.csv(paste(dataDir, "1362351439-TNET004-Plate1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable14$probe <- "TNET004"
  
  ctTable15 <- read.csv(paste(dataDir, "1362351440-TNET004-Plate2.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable16$probe <- "TNET004"
  
  
  
  
  ## everything
  ctTableCombine <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, 
                          ctTable9, ctTable10, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15)
  
  ## change column names
  names(ctTableCombine)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## Creat and name cellSource column value
  ctTableCombine$cellSource <- NA
  ctTableCombine[grepl("RAD", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "Child"
  ctTableCombine[grepl("TEY", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "Adult"
  ctTableCombine[grepl("TNET", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "Risk"
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  
  ## create age (inflamed/uninflamed) column
  # ctTableCombine[grepl("Infl", ctTableCombine$cellType, fixed=TRUE), "age"] <- "infl"
  # ctTableCombine[grepl("Uninf", ctTableCombine$cellType, fixed=TRUE), "age"] <- "uninfl"
  # ctTableCombine[grepl("Blood", ctTableCombine$cellType, fixed=TRUE), "age"] <- "blood"
  # ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "age"] <- "controlIslet"
  
  
  return(ctTableCombine)
}







