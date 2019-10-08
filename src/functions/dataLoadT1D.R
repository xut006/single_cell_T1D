#### load UC-MAIT data ####


dataLoadT1D <- function(){
  ## set data directory
  baseDir <- "~/Desktop/Su Lab/single_cell_T1D/"
  dataDir <- paste(baseDir, "data/", sep="")
  
  ### read in Children data
  ## probe -> patients
  ctTable1 <- read.csv(paste(dataDir, "1362015470-RAD001-FrozenvsFresh.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable1$probe <- "RAD1"
  
  ctTable2 <- read.csv(paste(dataDir, "RAD003.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable2$probe <- "RAD3"
  
  ctTable3 <- read.csv(paste(dataDir, "1362351438-RAD006.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable3$probe <- "RAD6" 
  
  ctTable4 <- read.csv(paste(dataDir, "1362351437-Processed_RAD005.csv.trimmed", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE)
  ctTable4$probe <- "RAD5"
  
  ctTable5 <- read.csv(paste(dataDir, "1362356331-RAD004-TEY013.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable5$probe <- NA
  ctTable5[grepl("RAD4", ctTable5$Name, fixed=TRUE), "probe"] <- "RAD4"
  ctTable5[grepl("RAD001", ctTable5$Name, fixed=TRUE), "probe"] <- "RAD1"
  ctTable5[grepl("TEY008", ctTable5$Name, fixed=TRUE), "probe"] <- "TEY8"
  ctTable5[grepl("TEY003", ctTable5$Name, fixed=TRUE), "probe"] <- "TEY3"
  
  ctTable6 <- read.csv(paste(dataDir, "1362356332-RAD003_Insulin_ICOS_PD1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable6$probe <- "RAD3"
  
  ctTable7 <- read.csv(paste(dataDir, "RAD001andRAD002.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable7$probe <- NA
  ctTable7[grepl("RAD1", ctTable7$Name, fixed=TRUE), "probe"] <- "RAD1"
  ctTable7[grepl("RAD2", ctTable7$Name, fixed=TRUE), "probe"] <- "RAD2"
  
  ctTable8 <- read.csv(paste(dataDir, "RAD003(p1220p1321).csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable8$probe <- "RAD3"
  
  
  ### read in Adults data
  ctTable9 <- read.csv(paste(dataDir, "1362292377-TEY006-TEY014.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable9$probe <- NA
  ctTable9[grepl("TEY-006", ctTable9$Name, fixed=TRUE), "probe"] <- "TEY6"
  ctTable9[grepl("TEY-014", ctTable9$Name, fixed=TRUE), "probe"] <- "TEY14"
  ctTable9[grepl("NBD", ctTable9$Name, fixed=TRUE), "probe"] <- "Healthy"
  
  ctTable10 <- read.csv(paste(dataDir, "1362356331-RAD004-TEY013.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable10 <- subset(ctTable10, grepl("TEY", ctTable10$Name, fixed = TRUE))
  ctTable10$probe <- NA
  ctTable10[grepl("TEY003", ctTable10$Name, fixed=TRUE), "probe"] <- "TEY3"
  ctTable10[grepl("TEY008", ctTable10$Name, fixed=TRUE), "probe"] <- "TEY8"
  
  
  ### read in At_Ristk data
  ctTable11 <- read.csv(paste(dataDir, "1362351434-processedresults-TNET002.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable11$probe <- "TNET2"
  
  ctTable12 <- read.csv(paste(dataDir, "1362351435_Processed-TNET001-UCM17.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable12 <- subset(ctTable12, grepl("TNet001", ctTable12$Name, fixed = TRUE))
  ctTable12$probe <- "TNET1"
  
  ctTable13 <- read.csv(paste(dataDir, "1362351436_Processed-TNET003.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable13$probe <- "TNET3"
  
  ctTable14 <- read.csv(paste(dataDir, "1362351439-TNET004-Plate1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable14$probe <- "TNET4"
  
  ctTable15 <- read.csv(paste(dataDir, "1362351440-TNET004-Plate2.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable15$probe <- "TNET4"
  
  
  
  
  ## everything
  ctTableCombine <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, 
                          ctTable9, ctTable10, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15)
  ## change column names
  names(ctTableCombine)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## Creat and name cellSource (cohort) column value
  ctTableCombine$cellSource <- NA
  ctTableCombine[grepl("RAD", ctTableCombine$probe, fixed=TRUE), "cellSource"] <- "Child"
  ctTableCombine[grepl("TEY", ctTableCombine$probe, fixed=TRUE), "cellSource"] <- "Adult"
  ctTableCombine[grepl("TNET", ctTableCombine$probe, fixed=TRUE), "cellSource"] <- "Risk"
  ctTableCombine[grepl("Healthy", ctTableCombine$probe, fixed=TRUE), "cellSource"] <- "Healthy"
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  ctTableCombine[grepl("EMTPY", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  ctTableCombine[grepl("Empty", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  ctTableCombine[grepl("empty", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  
  
  
  ## create age (different probes that extract cells) column
  ## age -> markers
  ctTableCombine$age <- NA
  # ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "age"] <- "control"
  # ctTableCombine[grepl("Empty", ctTableCombine$cellType, fixed=TRUE), "age"] <- "control"
  # ctTableCombine[grepl("empty", ctTableCombine$cellType, fixed=TRUE), "age"] <- "control"
  # ctTableCombine[grepl("EMTPY", ctTableCombine$cellType, fixed=TRUE), "age"] <- "control"
  
  ## Order of labeling is improtant
  
  ## Surface marker: PD1/CXCR3 double negtive & PD1/CXCR3 double positive
  ctTableCombine[grepl("DN", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/CXCR3-"
  ctTableCombine[grepl("DP", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/CXCR3+"
  ## Surface marker: PD1 positive
  ctTableCombine[grepl("PD1", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/ICOS-"
  ctTableCombine[grepl("PD1+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/ICOS-"
  ## Surface marker: ICOS positive
  ctTableCombine[grepl("ICOS", ctTableCombine$cellType, fixed=TRUE), "age"] <- "ICOS+/PD1-"
  ctTableCombine[grepl("ICOS+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "ICOS+/PD1-"
  ctTableCombine[grepl("Icosp", ctTableCombine$cellType, fixed=TRUE), "age"] <- "ICOS+/PD1-"
  ## Surface marker: CXCR3 positive
  ctTableCombine[grepl("CXCR3", ctTableCombine$cellType, fixed=TRUE), "age"] <- "CXCR3+/PD1-"
  ## Surface marker: PD1/ICOS double negtive
  ctTableCombine[grepl("IPn", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/ICOS-" 
  ctTableCombine[grepl("PD1-ICOS-DN", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/ICOS-"
  ctTableCombine[grepl("PD1-/ICOS-", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/ICOS-"
  ctTableCombine[grepl("PD1/ICOS-/-", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/ICOS-" ## 0
  ctTableCombine[grepl("PD1 ICOS DN", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/ICOS-" 
  ## Surface marker: PD1/ICOS double positive
  ctTableCombine[grepl("IPp", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/ICOS+"
  ctTableCombine[grepl("PD1+ICOS+ DP", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/ICOS+"
  ctTableCombine[grepl("PD1/ICOS +/+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/ICOS+" ## 0
  ctTableCombine[grepl("PD1 ICOS DP", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/ICOS+"
  ## Surface marker: PD1/CXCR3 double negtive
  ctTableCombine[grepl("CXCR3 PD1 DN", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/CXCR3-"
  ctTableCombine[grepl("PD1/CXCR3 -/-", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1-/CXCR3-" ## 0
  ## Surface marker: PD1/CXCR3 double positive
  ctTableCombine[grepl("CXCR3+ PD1+ DP", ctTableCombine$cellType, fixed=TRUE), "age"] <- "PD1+/CXCR3+"
  
  
  
  ## Tetramer markers: 12-20
  ctTableCombine[grepl("12-20+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "12-20"
  ctTableCombine[grepl("12-20", ctTableCombine$cellType, fixed=TRUE), "age"] <- "12-20"
  ctTableCombine[grepl("p12", ctTableCombine$cellType, fixed=TRUE), "age"] <- "12-20"
  ## Tetramer markers: 13-21
  ctTableCombine[grepl("13-21+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "13-21"
  ctTableCombine[grepl("13-21", ctTableCombine$cellType, fixed=TRUE), "age"] <- "13-21"
  ctTableCombine[grepl("p13", ctTableCombine$cellType, fixed=TRUE), "age"] <- "13-21"
  ## Tetramer markers: CP11
  ctTableCombine[grepl("CP11+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "CP11"
  ctTableCombine[grepl("CP11", ctTableCombine$cellType, fixed=TRUE), "age"] <- "CP11"
  ## Tetramer markers: CP13
  ctTableCombine[grepl("CP13+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "CP13"
  ## Tetramer markers: CP18
  ctTableCombine[grepl("CP18+", ctTableCombine$cellType, fixed=TRUE), "age"] <- "CP18"
  
  

  
  
  ## Test if the correct number of cells have been label with markers
  # markers <- c("PD1-/ICOS-", "PD1+/ICOS+", "PD1+/ICOS-", "ICOS+/PD1-", "PD1-/CXCR3-", "PD1+/CXCR3+", "CXCR3+/PD1-",
  #              "12-20", "13-21", "CP11", "CP13", "CP18")
  # sum <- 0
  # for (i in markers){
  #   num <- length(ctTableCombine[grepl(i, ctTableCombine$age, fixed=TRUE), "age"])
  #   result <- cat("marker ", i, " have ", num/96, " cells")
  #   print(result)
  #   sum <- sum + num
  # }
  # print (sum)
  # print (sum/96)
  # 
  # control <- length(ctTableCombine[grepl("control", ctTableCombine$cellSource, fixed=TRUE), "cellSource"])
  # sum <- sum + control
  # print (sum)
  
  
  
  
  return(ctTableCombine)
}







