#### load T1D data ####

## preprocess the datafram by: sed '1,11d' *.csv > *.csv.trimmed
## OR preprocess the datafram by: sh file_trimming.sh 

dataLoadT1D <- function(){
  ## set data directory
  baseDir <- "~/Desktop/Su Lab/single_cell_T1D/"
  dataDir <- paste(baseDir, "data/", sep="")
  
  ### read in Children data
  ctTable1 <- read.csv(paste(dataDir, "1362015470-RAD001-FrozenvsFresh.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable1$patients <- "RAD1"
  
  ctTable2 <- read.csv(paste(dataDir, "RAD003.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable2$patients <- "RAD3"
  
  ctTable3 <- read.csv(paste(dataDir, "1362351438-RAD006.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable3$patients <- "RAD6" 
  
  ctTable4 <- read.csv(paste(dataDir, "1362351437-Processed_RAD005.csv.trimmed", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE)
  ctTable4$patients <- "RAD5"
  
  ctTable5 <- read.csv(paste(dataDir, "1362356331-RAD004-TEY013.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable5$patients <- NA
  ctTable5[grepl("RAD4", ctTable5$Name, fixed=TRUE), "patients"] <- "RAD4"
  ctTable5[grepl("RAD001", ctTable5$Name, fixed=TRUE), "patients"] <- "RAD1"
  ctTable5[grepl("TEY008", ctTable5$Name, fixed=TRUE), "patients"] <- "TEY8"
  ctTable5[grepl("TEY003", ctTable5$Name, fixed=TRUE), "patients"] <- "TEY3"
  
  ctTable6 <- read.csv(paste(dataDir, "1362356332-RAD003_Insulin_ICOS_PD1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable6$patients <- "RAD3"
  
  ctTable7 <- read.csv(paste(dataDir, "RAD001andRAD002.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable7$patients <- NA
  ctTable7[grepl("RAD1", ctTable7$Name, fixed=TRUE), "patients"] <- "RAD1"
  ctTable7[grepl("RAD2", ctTable7$Name, fixed=TRUE), "patients"] <- "RAD2"
  
  ctTable8 <- read.csv(paste(dataDir, "RAD003(p1220p1321).csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable8$patients <- "RAD3"
  
  ctTable16 <- read.csv(paste(dataDir, "1362427424_RAD005.csv.trimmed", sep=""),
                      header=TRUE, stringsAsFactors=FALSE)
  ctTable16$patients <- "RAD5"
  
  
  ### read in Adults data
  ctTable9 <- read.csv(paste(dataDir, "1362292377-TEY006-TEY014.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable9$patients <- NA
  ctTable9[grepl("TEY-006", ctTable9$Name, fixed=TRUE), "patients"] <- "TEY6"
  ctTable9[grepl("TEY-014", ctTable9$Name, fixed=TRUE), "patients"] <- "TEY14"
  ctTable9[grepl("NBD", ctTable9$Name, fixed=TRUE), "patients"] <- "Healthy"
  
  ctTable10 <- read.csv(paste(dataDir, "1362356331-RAD004-TEY013.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable10 <- subset(ctTable10, grepl("TEY", ctTable10$Name, fixed = TRUE))
  ctTable10$patients <- NA
  ctTable10[grepl("TEY003", ctTable10$Name, fixed=TRUE), "patients"] <- "TEY3"
  ctTable10[grepl("TEY008", ctTable10$Name, fixed=TRUE), "patients"] <- "TEY8"
  
  ctTable17 <- read.csv(paste(dataDir, "1362427423_T004tetramers.csv.trimmed", sep=""),
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable17$patients <- "TEY4"
  
  
  ### read in At_Ristk data
  ctTable11 <- read.csv(paste(dataDir, "1362351434-processedresults-TNET002.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable11$patients <- "TNET2"
  
  ctTable12 <- read.csv(paste(dataDir, "1362351435_Processed-TNET001-UCM17.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable12 <- subset(ctTable12, grepl("TNet001", ctTable12$Name, fixed = TRUE))
  ctTable12$patients <- "TNET1"
  
  ctTable13 <- read.csv(paste(dataDir, "1362351436_Processed-TNET003.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable13$patients <- "TNET3"
  
  ctTable14 <- read.csv(paste(dataDir, "1362351439-TNET004-Plate1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable14$patients <- "TNET4"
  
  ctTable15 <- read.csv(paste(dataDir, "1362351440-TNET004-Plate2.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable15$patients <- "TNET4"
  
  ctTable20 <- read.csv(paste(dataDir, "1362438043_TNet005_Plate1.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable20$patients <- "TNET5"
  
  ctTable21 <- read.csv(paste(dataDir, "1362427237_TNET006-Plate1.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable21$patients <- "TNET6"
  
  ctTable22 <- read.csv(paste(dataDir, "1362427296_TNET006-Plate2.csv.trimmed", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable22$patients <- "TNET6"
  
  
  ### read in Normal Blood Donor data
  ctTable18 <- read.csv(paste(dataDir, "1362427412-TC017.csv.trimmed", sep=""),
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable18$patients <- "Healthy"

  ctTable19 <- read.csv(paste(dataDir, "1362427336_T-C012.csv.trimmed", sep=""),
                        header=TRUE, stringsAsFactors=FALSE)
  ctTable19$patients <- "Healthy"
  
  
  
  ## everything
  ctTableCombine <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, 
                          ctTable9, ctTable10, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15,
                          ctTable16, ctTable17, ctTable18, ctTable19, ctTable20, ctTable21, ctTable22)
  ## change column names
  names(ctTableCombine)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## Creat and name cohort (cohort) column value
  ctTableCombine$cohort <- NA
  ctTableCombine[grepl("RAD", ctTableCombine$patients, fixed=TRUE), "cohort"] <- "Child"
  ctTableCombine[grepl("TEY", ctTableCombine$patients, fixed=TRUE), "cohort"] <- "Adult"
  ctTableCombine[grepl("TNET", ctTableCombine$patients, fixed=TRUE), "cohort"] <- "Risk"
  ctTableCombine[grepl("Healthy", ctTableCombine$patients, fixed=TRUE), "cohort"] <- "Healthy"
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "cohort"] <- "control"
  ctTableCombine[grepl("EMTPY", ctTableCombine$cellType, fixed=TRUE), "cohort"] <- "control"
  ctTableCombine[grepl("Empty", ctTableCombine$cellType, fixed=TRUE), "cohort"] <- "control"
  ctTableCombine[grepl("empty", ctTableCombine$cellType, fixed=TRUE), "cohort"] <- "control"
  
  
  
  ## create marker column
  ctTableCombine$marker <- NA
  
  ## create Type of marker column (Tetramer/Surface)
  ctTableCombine$markerType <- NA
  
  ## Order of labeling is improtant!!!
  
  ## Surface marker: PD1/CXCR3 double negtive & PD1/CXCR3 double positive
  ctTableCombine[grepl("DN", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/CXCR3-"
  ctTableCombine[grepl("DP", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/CXCR3+"
  ctTableCombine[grepl("PD1-/CXCR3-", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ctTableCombine[grepl("PD1-/CXCR3+", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: PD1 positive
  ctTableCombine[grepl("PD1", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/ICOS-"
  ctTableCombine[grepl("PD1+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/ICOS-"
  ctTableCombine[grepl("PD1+/ICOS-", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: ICOS positive
  ctTableCombine[grepl("ICOS", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "ICOS+/PD1-"
  ctTableCombine[grepl("ICOS+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "ICOS+/PD1-"
  ctTableCombine[grepl("Icosp", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "ICOS+/PD1-"
  ctTableCombine[grepl("ICOS+/PD1-", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: CXCR3 positive
  ctTableCombine[grepl("CXCR3", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CXCR3+/PD1-"
  ctTableCombine[grepl("CXCR3+/PD1-", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: PD1/ICOS double negtive
  ctTableCombine[grepl("IPn", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/ICOS-" 
  ctTableCombine[grepl("PD1-ICOS-DN", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/ICOS-"
  ctTableCombine[grepl("PD1-/ICOS-", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/ICOS-"
  ctTableCombine[grepl("PD1/ICOS-/-", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/ICOS-" ## 0
  ctTableCombine[grepl("PD1 ICOS DN", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/ICOS-" 
  ctTableCombine[grepl("PD1-/ICOS-", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: PD1/ICOS double positive
  ctTableCombine[grepl("IPp", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/ICOS+"
  ctTableCombine[grepl("PD1+ICOS+ DP", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/ICOS+"
  ctTableCombine[grepl("PD1/ICOS +/+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/ICOS+" ## 0
  ctTableCombine[grepl("PD1 ICOS DP", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/ICOS+"
  ctTableCombine[grepl("PD1+/ICOS+", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: PD1/CXCR3 double negtive
  ctTableCombine[grepl("CXCR3 PD1 DN", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/CXCR3-"
  ctTableCombine[grepl("PD1/CXCR3 -/-", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1-/CXCR3-" ## 0
  ctTableCombine[grepl("PD1-/CXCR3-", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  ## Surface marker: PD1/CXCR3 double positive
  ctTableCombine[grepl("CXCR3+ PD1+ DP", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "PD1+/CXCR3+"
  ctTableCombine[grepl("PD1+/CXCR3+", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Surface"
  
  
  
  ## Tetramer markers: 12-20
  ctTableCombine[grepl("12-20+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "12-20"
  ctTableCombine[grepl("12-20", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "12-20"
  ctTableCombine[grepl("p12", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "12-20"
  ctTableCombine[grepl("12-20", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  ## Tetramer markers: 13-21
  ctTableCombine[grepl("13-21+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "13-21"
  ctTableCombine[grepl("13-21", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "13-21"
  ctTableCombine[grepl("p13", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "13-21"
  ctTableCombine[grepl("13-21", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  ## Tetramer markers: CP11
  ctTableCombine[grepl("CP11+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CP11"
  ctTableCombine[grepl("CP11", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CP11"
  ctTableCombine[grepl("CP11", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  ## Tetramer markers: CP13
  ctTableCombine[grepl("CP13+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CP13"
  ctTableCombine[grepl("CP13", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CP13"
  ctTableCombine[grepl("CP13", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  ## Tetramer markers: CP18
  ctTableCombine[grepl("CP18+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CP18"
  ctTableCombine[grepl("CP18", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CP18"
  ctTableCombine[grepl("CP18", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  ## Tetramer markers: Alt 1-9 
  ctTableCombine[grepl("Alt_1-9", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "Alt1-9"
  ctTableCombine[grepl("Alt1-9", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  ## Tetramer markers: CLIP 
  ctTableCombine[grepl("CLIP", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "CLIP"
  ctTableCombine[grepl("CLIP", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Tetramer"
  
  
  
  ## Bulk_CD4+_T_cell Markers
  ctTableCombine[grepl("Bulk_CD4+", ctTableCombine$cellType, fixed=TRUE), "marker"] <- "Bulk_CD4+"
  ctTableCombine[grepl("Bulk_CD4+", ctTableCombine$marker, fixed=TRUE), "markerType"] <- "Bulk_CD4+"
  

  
  
  ### Check if all new columns have value assignmed
  # sum(is.na(ctTableCombine$marker))
  # na_ctTable<-subset(ctTableCombine,is.na(ctTableCombine$marker))
  # sum(na_ctTable$cohort != "control")
  # na_ctTable<-subset(na_ctTable,na_ctTable$cohort != "control")
  
  
  
  
  return(ctTableCombine)
}







