#### cluster and filter to remove outliers

clusterFilter <- function(ctInput, testK = F, numCenters = 2, plotHeatmap=F, plotClustOnly=F, 
                          heatmapFactor, heatmapColorBy, heatmapcohortLabel,
                          # fisherTests = c("patients", "cohort", "marker", "markerType"),
                          fisherTests = NULL,
                          cumulativeExpHist = T, filterClusters = F, clustersToRemove = NULL,
                          removalText = T){
  
  #### kmeans analysis ####
  ## order by ctSum
  ctInput$ctSum <- rowSums(ctInput[,9:ncol(ctInput)])
  ctInput <- ctInput[order(ctInput$ctSum) ,]
  ctInput$ctSum <- NULL
  
  ctTestK <- ctInput[,9:ncol(ctInput)]
  
  ## set seed
  set.seed(1)
  
  if(testK == T){
    ## run kmeans with 1-15 centers
    wss <- (nrow(ctTestK)-1)*sum(apply(ctTestK,2,var))
    for(i in 2:15){ 
      wss[i] <- sum(kmeans(ctTestK, centers=i)$withinss)
    }
    
    ## create within sum of squares dataframe and plot
    centers <- 1:15
    wssc <- as.data.frame(cbind(centers, wss))
    wssPlot <- ggplot(wssc, aes(centers, wss)) +
      geom_point(colour = "blue", fill = "cyan", size = 5, shape=21) +
      geom_path(size = 2.3) +
      ggtitle("kmeans scree plot") +
      ylab("total within-cluster sum of squares") +
      scale_x_continuous(breaks=1:15) +
      theme(text=element_text(size=25),
            axis.title.x=element_text(vjust=-0.5),
            plot.margin=unit(c(10,5,10,5),"cm"))
    print(wssPlot)
  }
  
  ## run kmeans with chosen number of centers
  set.seed(1)
  normFit <- kmeans(ctTestK, numCenters)
  
  ## attach results to dataframe
  ctClust <- data.frame(ctInput, normFit$cluster)
  
  ## rename clustering results column
  names(ctClust)[which(names(ctClust)=="normFit.cluster")] <- "kmeans.cluster"
  
  ## move kmeans column to column 9 to keep with id variables
  ctClust <- ctClust[, c(1:8, ncol(ctClust), 9:(ncol(ctClust)-1))]
  
  
  #### heatmap and cluster membership ####
  if(plotHeatmap == T){
    par(mar=c(1,1,1,1))
    
    ## account for kmeans.cluster column when in use
    if(heatmapFactor != "kmeans.cluster"){
      idCols <- as.numeric(8)
      ctClust <- subset(ctClust, select = -kmeans.cluster)
    } else{
      idCols <- as.numeric(9)
    }
    
    ## get p-values for differential expression for factors
    pvals <- NULL
    for(i in (idCols+1):ncol(ctClust)){
      pvals <- c(pvals, kruskal.test(ctClust[,i], factor(ctClust[, heatmapFactor]))$p.value)
    }
    ctClust <- ctClust[,c(1:idCols, (order(pvals)+idCols))]
    
    
    AnnoColors <- NULL
    ## label cells by marker
    ## marker -> markers
    if("marker" %in% heatmapColorBy){
      ctClust$markerColor <- NA
      ctClust[which(ctClust$marker=="12-20"), "markerColor"] <- "#9A8822"
      ctClust[which(ctClust$marker=="13-21"), "markerColor"] <- "#F5CDB4"
      ctClust[which(ctClust$marker=="CP11"), "markerColor"] <- "#F8AFA8"
      ctClust[which(ctClust$marker=="CP13"), "markerColor"] <- "#FDDDA0"
      ctClust[which(ctClust$marker=="CP18"), "markerColor"] <- "#74A089"
      ctClust[which(ctClust$marker=="Alt1-9"), "markerColor"] <- "royalblue3"
      ctClust[which(ctClust$marker=="CLIP"), "markerColor"] <- "#9986A5"
      ctClust[which(ctClust$marker=="CXCR3+/PD1-"), "markerColor"] <- "#FF0000"
      ctClust[which(ctClust$marker=="ICOS+/PD1-"), "markerColor"] <- "#556A5B"
      ctClust[which(ctClust$marker=="PD1-/CXCR3-"), "markerColor"] <- "#50A45C"
      ctClust[which(ctClust$marker=="PD1-/ICOS-"), "markerColor"] <- "#F2AD00"
      ctClust[which(ctClust$marker=="PD1+/CXCR3+"), "markerColor"] <- "#F69100"
      ctClust[which(ctClust$marker=="PD1+/ICOS-"), "markerColor"] <- "#C49647"
      ctClust[which(ctClust$marker=="PD1+/ICOS+"), "markerColor"] <- "#5BBCD6"
      ctClust[which(ctClust$marker=="Bulk_CD4+"), "markerColor"] <- "moccasin"
      markers <- ctClust$markerColor
      ctClust <- subset(ctClust, select=-c(markerColor))
      AnnoColors <- cbind(AnnoColors, markers)
    }
    
    ## label cells by cohort: Child, Adult, Risk
    if("cohort" %in% heatmapColorBy){
      ctClust$cohortColor <- NA
      ctClust[which(ctClust$cohort=="Child"), "cohortColor"] <- "#FDE725FF"
      ctClust[which(ctClust$cohort=="Adult"), "cohortColor"] <- "#21908CFF"
      ctClust[which(ctClust$cohort=="Risk"), "cohortColor"] <- "#440154FF"
      ctClust[which(ctClust$cohort=="Healthy"), "cohortColor"] <- "gray48"
      cohorts <- ctClust$cohortColor
      ctClust <- subset(ctClust, select=-c(cohortColor))
      AnnoColors <- cbind(AnnoColors, cohorts)
    }
    
    ## label cells by markerType: Surface, Tetramer, Bulk_CD4+
    if("markerType" %in% heatmapColorBy){
      ctClust$markerTypeColor <- NA
      ctClust[which(ctClust$markerType=="Surface"), "markerTypeColor"] <- "#85D4E3"
      ctClust[which(ctClust$markerType=="Tetramer"), "markerTypeColor"] <- "#F4B5BD"
      ctClust[which(ctClust$markerType=="Bulk_CD4+"), "markerTypeColor"] <- "#9C964A"
      markerTypes <- ctClust$markerTypeColor
      ctClust <- subset(ctClust, select=-c(markerTypeColor))
      AnnoColors <- cbind(AnnoColors, markerTypes)
    }
    
    ## label cells by patients
    ## patients -> patients
    if("patients" %in% heatmapColorBy){
      ctClust$patientsColor <- NA
      ctClust[which(ctClust$patients=="Healthy"), "patientsColor"] <- "gray19"
      ctClust[which(ctClust$patients=="RAD1"), "patientsColor"] <- "#DD8D29"
      ctClust[which(ctClust$patients=="RAD2"), "patientsColor"] <- "#E1C408"
      ctClust[which(ctClust$patients=="RAD3"), "patientsColor"] <- "#84BB78"
      ctClust[which(ctClust$patients=="RAD4"), "patientsColor"] <- "#859C78"
      ctClust[which(ctClust$patients=="RAD5"), "patientsColor"] <- "#DB6E07"
      ctClust[which(ctClust$patients=="RAD6"), "patientsColor"] <- "#B40F20"
      ctClust[which(ctClust$patients=="TEY14"), "patientsColor"] <- "#798E87"
      ctClust[which(ctClust$patients=="TEY3"), "patientsColor"] <- "#C27D38"
      ctClust[which(ctClust$patients=="TEY4"), "patientsColor"] <- "violetred4"
      ctClust[which(ctClust$patients=="TEY6"), "patientsColor"] <- "#CCC591"
      ctClust[which(ctClust$patients=="TEY8"), "patientsColor"] <- "#29211F"
      ctClust[which(ctClust$patients=="TNET1"), "patientsColor"] <- "#F3DF6C"
      ctClust[which(ctClust$patients=="TNET2"), "patientsColor"] <- "#CEAB07"
      ctClust[which(ctClust$patients=="TNET3"), "patientsColor"] <- "#D5D5D3"
      ctClust[which(ctClust$patients=="TNET4"), "patientsColor"] <- "#24281A"
      patientss <- ctClust$patientsColor
      ctClust <- subset(ctClust, select=-c(patientsColor))
      AnnoColors <- cbind(AnnoColors, patientss)
    }
    
    ## label cells by cluster
    if(heatmapFactor == "kmeans.cluster"){
      
      ## for heatmap figure
      # clusterColorKey <- brewer.pal(8, "Dark2")
      # extraColorPalette <- clusterColorKey[c(2,3,4,6)]
      
      extraColorPalette <- brewer.pal(12, "Set3")
      clusterVals <- unique(ctClust$kmeans.cluster)
      clusterVals <- clusterVals[order(clusterVals)]
      
      ctClust$kmeansColor <- NULL
      for(i in 1:length(clusterVals)){
        ctClust[which(ctClust$kmeans.cluster == clusterVals[i]), 
                "kmeansColor"] <- extraColorPalette[i]
      }
      
      Kmeans.clusters <- ctClust$kmeansColor
      ctClust$kmeansColor <- NULL
      AnnoColors <- cbind(AnnoColors, Kmeans.clusters)
    }
    
    markerColorKey <- rbind(c("12-20", "#9A8822"),
                         c("13-21", "#F5CDB4"),
                         c("CP11", "#F8AFA8"),
                         c("CP13", "#FDDDA0"),
                         c("CP18", "#74A089"),
                         c("Alt1-9", "royalblue3"),
                         c("CLIP", "#9986A5"),
                         c("CXCR3+/PD1-", "#FF0000"),
                         c("ICOS+/PD1-", "#556A5B"),
                         c("PD1-/CXCR3-", "#50A45C"),
                         c("PD1-/ICOS-", "#F2AD00"),
                         c("PD1+/CXCR3+", "#F69100"),
                         c("PD1+/ICOS-", "#C49647"),
                         c("PD1+/ICOS+", "#5BBCD6"),
                         c("Bulk_CD4+", "moccasin"))
    markerColorKey <- markerColorKey[which(markerColorKey[,1] %in% ctClust$marker),]
    
    cohortColorKey <- rbind(c("Child", "#FDE725FF"),
                            c("Adult","#21908CFF"),
                            c("Risk","#440154FF"),
                            c("Healthy", "gray48"))
    cohortColorKey <- cohortColorKey[which(cohortColorKey[,1] %in% ctClust$cohort),]
    
    markerTypeColorKey <- rbind(c("Surface","#85D4E3"),
                            c("Tetramer","#F4B5BD"),
                            c("Bulk_CD4+", "#9C964A"))
    markerTypeColorKey <- markerTypeColorKey[which(markerTypeColorKey[,1] %in% ctClust$markerType),]
    
    patientsColorKey <- rbind(c("Healthy", "gray19"),
                           c("RAD1", "#DD8D29"),
                           c("RAD2", "#E1C408"),
                           c("RAD3", "#84BB78"),
                           c("RAD4", "#859C78"),
                           c("RAD5", "#DB6E07"),
                           c("RAD6", "#B40F20"),
                           c("TEY14", "#798E87"),
                           c("TEY3", "#C27D38"),
                           c("TEY4", "violetred4"),
                           c("TEY6", "#CCC591"),
                           c("TEY8", "#29211F"),
                           c("TNET1", "#F3DF6C"),
                           c("TNET2", "#CEAB07"),
                           c("TNET3", "#D5D5D3"),
                           c("TNET4", "#24281A"))
    patientsColorKey <- patientsColorKey[which(patientsColorKey[,1] %in% ctClust$patients),]
    
    insetX <- 0.256
    
      
#       clusterVals <- unique(ctClust$kmeans.cluster)
#       clusterVals <- clusterVals[order(clusterVals)]
#       
#       if(length(clusterVals) == 6){
#         extremeClusterVals <- c(-15, -9, -4, 4, 9, 15)
#       } else{
#         extremeClusterVals <- seq(-15, 15, (30/(length(clusterVals) - 1)))
#       }
#         
#       ctClust$extremeKmeans <- NULL
#       for(i in 1:length(clusterVals)){
#         ctClust[which(ctClust$kmeans.cluster == clusterVals[i]), 
#                 "extremeKmeans"] <- extremeClusterVals[i]
#       }
#       
#       kmeans.cluster <- ctClust$kmeans.cluster
#       ctClust$kmeans.cluster <- ctClust$extremeKmeans
#       ctClust$extremeKmeans <- NULL
#     }
    
    ## plot heatmap
    par(cex.main = 1.5)
#     heatmap.2(as.matrix(ctClust[(idCols + 1):ncol(ctClust)]), Colv=F, 
#               dendrogram="row", trace="none", 
#               col=colorpanel(15, "red", "white", "blue"),
#               RowSideColors = extraColors,
#               colRow=rowColors, adjRow = heatmapRowAdjust,
#               adjCol = c(1, NA), cexCol = 1.3, cexRow=heatmapRowSize, margins=c(10,10), 
#               density.info="density", denscol="black", key.title=NA, 
#               main=paste("heatmap for", heatmapcohortLabel, sep=" "), 
#               xlab="log2(gene expression / Gapdh expression)", 
#               ylab=heatmapYLabel)
    if(plotClustOnly == FALSE){
      heatmap.3(as.matrix(ctClust[(idCols + 1):ncol(ctClust)]), Colv=F, 
                dendrogram="row", trace="none", 
                col=colorpanel(15, "blue", "white", "red"),
                RowSideColors = AnnoColors,
                labRow=NA,
                adjCol = c(1, NA), cexCol = 1.3, margins=c(10,10), 
                density.info="density", denscol="black", key.title=NA, 
                main=paste("heatmap for", heatmapcohortLabel, sep=" "), 
                xlab="log2(gene expression / Gapdh expression)", 
                ylab="Individual Cells\n\n\n")
      if("marker" %in% heatmapColorBy){
        legend("topleft",
               legend=markerColorKey[,1],
               fill=markerColorKey[,2], 
               title="marker",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
        insetX <- insetX + 0.12
      }
      if("cohort" %in% heatmapColorBy){
        legend("topleft",
               legend=cohortColorKey[,1],
               fill=cohortColorKey[,2], 
               title="cohort",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
        insetX <- insetX + 0.12
      }
      if("markerType" %in% heatmapColorBy){
        legend("topleft",
               legend=markerTypeColorKey[,1],
               fill=markerTypeColorKey[,2], 
               title="markerType",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
        insetX <- insetX + 0.12
      }
      if("patients" %in% heatmapColorBy){
        legend("topleft",
               legend=patientsColorKey[,1],
               fill=patientsColorKey[,2], 
               title="patients",
             border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
      insetX <- insetX + 0.075
        #         border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.1))
        # insetX <- insetX + 0.08
      }
      if(heatmapFactor == "kmeans.cluster"){
        legend("topleft",
               legend=clusterVals,
               fill=extraColorPalette[1:length(clusterVals)], 
               title="cluster",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
      }
    }
      
      if(plotClustOnly == TRUE){
        heatmap.2(as.matrix(ctClust[(idCols + 1):ncol(ctClust)]), Colv=F, 
                  dendrogram="row", trace="none", 
                  col=colorpanel(15, "blue", "white", "red"),
                  RowSideColors = Kmeans.clusters,
                  labRow=NA,
                  adjCol = c(1, NA), cexCol = 1.3, margins=c(10,10), 
                  density.info="density", denscol="black", key.title=NA, 
                  main=paste("heatmap for", heatmapcohortLabel, sep=" "), 
                  xlab="log2(gene expression / Gapdh expression)", 
                  ylab="Individual Cells\n\n\n")
        legend("topleft",
               legend=clusterVals,
               fill=extraColorPalette[1:length(clusterVals)], 
               title="cluster",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(0.26,0.15))
      }
    }
    par(cex.main = 1.2)
    
    
    ## remap clusters values to original
#     if(heatmapFactor == "kmeans.cluster"){
#       ctClust$kmeans.cluster <- kmeans.cluster
#     }
    
    
    #### fisher's exact test for all covariates ####
    if(heatmapFactor == "kmeans.cluster"){
      ## get cluster numbers
      clusterVals <- unique(ctClust$kmeans.cluster)
      clusterVals <- clusterVals[order(clusterVals)]
      
      ## patients
      if("patients" %in% fisherTests){
        patientss <- unique(ctClust$patients)
        patientss <- patientss[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY4","TEY6",
                                       "TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), 
                                     patientss, nomatch=F)]
        patientsTable <- matrix(nrow=length(clusterVals), ncol=length(patientss))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(patientss)){
            patientsTable[i, j] <- nrow(subset(ctClust, patients==patientss[j] & 
                                                      kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(patientsTable) <- paste("patients_", patientss, sep="")
        rownames(patientsTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient vs. Cluster", quote=F)
        # print(patientsTable)
        print(patientsTable, quote = F)
        # print(fisher.test(factor(ctClust$patients), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(patientsTable, simulate.p.value=T))
        
        ## Bar Plots
        patientsTablePlots <- tableBarPlots(plotFactor = "patients", plotFactorVals = patientss, clusterVals = clusterVals, factorTable = patientsTable, colorKey = patientsColorKey)
      }
      
      
      ## patient cohorts
      if("cohort" %in% fisherTests){
        cohorts <- unique(ctClust$cohort)
        cohorts <- cohorts[match(c("Child","Adult","Risk","Healthy"), cohorts, nomatch=F)]
        cohortTable <- matrix(nrow=length(clusterVals), ncol=length(cohorts))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(cohorts)){
            cohortTable[i, j] <- nrow(subset(ctClust, cohort==cohorts[j] & 
                                               kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(cohortTable) <- paste("cohort_", cohorts, sep="")
        rownames(cohortTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Cohort vs. Cluster", quote=F)
        print(cohortTable)
        # print(fisher.test(factor(ctClust$cohort), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(cohortTable, simulate.p.value=T))
        
        ## Bar Plots
        cohortTablePlots <- tableBarPlots(plotFactor = "cohort", plotFactorVals = cohorts, clusterVals = clusterVals, factorTable = cohortTable, colorKey = cohortColorKey)
        # print(cohortTablePlots[[4]])
      }
      
      ## all markers
      if("marker" %in% fisherTests){
        markers <- unique(ctClust$marker)
        markers <- markers[match(c("12-20", "13-21", "CP11","CP13", "CP18", "Alt1-9", "CLIP",
                                   "CXCR3+/PD1-", "ICOS+/PD1-", "PD1-/CXCR3-", 
                                   "PD1-/ICOS-", "PD1+/CXCR3+", "PD1+/ICOS-", 
                                   "PD1+/ICOS+", "Bulk_CD4+"), 
                                 markers, nomatch=F)]
        markerTable <- matrix(nrow=length(clusterVals), ncol=length(markers))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(markers)){
            markerTable[i, j] <- nrow(subset(ctClust, marker==markers[j] & 
                                               kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(markerTable) <- paste("marker_", markers, sep="")
        rownames(markerTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("marker vs. Cluster", quote=F)
        print(markerTable)
        # print(fisher.test(factor(ctClust$marker), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(markerTable, simulate.p.value=T))
        
        ## Bar Plots
        markerTablePlots <- tableBarPlots(plotFactor = "marker", plotFactorVals = markers, clusterVals = clusterVals, factorTable = markerTable, colorKey = markerColorKey)
        # print(markerTablePlots[[4]])
      }
      
      ## marker Types
      if("markerType" %in% fisherTests){
        markerTypes <- unique(ctClust$markerType)
        markerTypes <- markerTypes[match(c("Surface", "Tetramer", "Bulk_CD4+"), markerTypes, nomatch=F)]
        markerTypeTable <- matrix(nrow=length(clusterVals), ncol=length(markerTypes))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(markerTypes)){
            markerTypeTable[i, j] <- nrow(subset(ctClust, markerType==markerTypes[j] & 
                                               kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(markerTypeTable) <- paste("markerType_", markerTypes, sep="")
        rownames(markerTypeTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("markerType vs. Cluster", quote=F)
        print(markerTypeTable)
        # print(fisher.test(factor(ctClust$markerType), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(markerTypeTable, simulate.p.value=T))
        
        ## Bar Plots
        markerTypeTablePlots <- tableBarPlots(plotFactor = "markerType", plotFactorVals = markerTypes, clusterVals = clusterVals, factorTable = markerTypeTable, colorKey = markerTypeColorKey)
        # print(markerTypeTablePlots[[4]])
      }
      
      
      ng <- nullGrob()
      if("patients" %in% fisherTests){}
#       if("marker" %in% fisherTests){}
#       if("cohort" %in% fisherTests){}
      
      if("markerType" %in% fisherTests & "cohort" %in% fisherTests & "marker" %in% fisherTests){
        grid.arrange(grobs = list(markerTypeTablePlots[[3]], markerTypeTablePlots[[4]], 
                                  cohortTablePlots[[3]], cohortTablePlots[[4]], 
                                  markerTablePlots[[3]], markerTablePlots[[4]]), nrow=3)
        
        grid.arrange(grobs = list(markerTypeTablePlots[[1]], markerTypeTablePlots[[2]], 
                                  cohortTablePlots[[1]], cohortTablePlots[[2]], 
                                  markerTablePlots[[1]], markerTablePlots[[2]]), nrow=3)
      }
      
      else if("patients" %in% fisherTests & "marker" %in% fisherTests){
        grid.arrange(grobs = list(patientsTablePlots[[3]], patientsTablePlots[[4]], markerTablePlots[[3]], markerTablePlots[[4]], ng, ng), nrow=3)
        
        grid.arrange(grobs = list(patientsTablePlots[[1]], patientsTablePlots[[2]], markerTablePlots[[1]], markerTablePlots[[2]], ng, ng), nrow=3)
      }
      
      ## patients and marker
      if("patients.marker" %in% fisherTests){
        ctClust$patients.marker <- paste(ctClust$patients, ctClust$marker, sep = ".")
        patientss <- unique(ctClust$patients)
        patientss <- patientss[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY4","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), patientss, nomatch=F)]
        markers <- unique(ctClust$marker)
        markers <- markers[order(as.numeric(markers))]
        eg <- expand.grid(patientss, markers)
        patientss.markers <- sprintf('%s.%s', eg[,1], eg[,2])
        
        patients.markerTable <- matrix(nrow=length(clusterVals), ncol=length(patientss.markers))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(patientss.markers)){
            patients.markerTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                          patients.marker==patientss.markers[j]))
          }
        }
        colnames(patients.markerTable) <- patientss.markers
        rownames(patients.markerTable) <- paste("cluster_", clusterVals, sep="")
        patients.markerTable <- patients.markerTable[, which(colSums(patients.markerTable) > 0)]
        patientss.markers <- colnames(patients.markerTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient.Marker vs. Cluster", quote=F)
        # print(patientsTable)
        print(patients.markerTable, quote = F)
        # print(fisher.test(factor(ctClust$patients), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(patients.markerTable, simulate.p.value=T))
        
        ## Bar Plots
        patients.markerColors <- brewer.pal(12, "Paired")
        patients.markerColorKey <- cbind(patientss.markers, patients.markerColors[1:length(patientss.markers)])
        
        patients.markerTablePlots <- tableBarPlots(plotFactor = "patients.marker", plotFactorVals = patientss.markers, clusterVals = clusterVals, factorTable = patients.markerTable, colorKey = patients.markerColorKey)
        
        grid.arrange(grobs = list(patients.markerTablePlots[[1]], patients.markerTablePlots[[2]],  ng, ng, ng, ng), nrow=3)
        
        ctClust <- subset(ctClust, select = -patients.marker)
      }
      
      ## patients and cohort
      if("patients.cohort" %in% fisherTests){
        ctClust$patients.cohort <- paste(ctClust$patients, ctClust$cohort, sep = ".")
        patientss <- unique(ctClust$patients)
        patientss <- patientss[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY4","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), patientss, nomatch=F)]
        cohorts <- unique(ctClust$cohort)
        cohorts <- cohorts[match(c("Child","Adult","Risk","Healthy"), cohorts, nomatch=F)]
        eg <- expand.grid(patientss, cohorts)
        patientss.cohorts <- sprintf('%s.%s', eg[,1], eg[,2])
        
        patients.cohortTable <- matrix(nrow=length(clusterVals), ncol=length(patientss.cohorts))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(patientss.cohorts)){
            patients.cohortTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                     patients.cohort==patientss.cohorts[j]))
          }
        }
        colnames(patients.cohortTable) <- patientss.cohorts
        rownames(patients.cohortTable) <- paste("cluster_", clusterVals, sep="")
        patients.cohortTable <- patients.cohortTable[, which(colSums(patients.cohortTable) > 0)]
        patientss.cohorts <- colnames(patients.cohortTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient.Cohort vs. Cluster", quote=F)
        # print(patientsTable)
        print(patients.cohortTable, quote = F)
        # print(fisher.test(factor(ctClust$patients), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(patients.cohortTable, simulate.p.value=T))
        
        ## Bar Plots
        patients.cohortColors <- brewer.pal(12, "Paired")
        patients.cohortColorKey <- cbind(patientss.cohorts, patients.cohortColors[1:length(patientss.cohorts)])
        
        patients.cohortTablePlots <- tableBarPlots(plotFactor = "patients.cohort", plotFactorVals = patientss.cohorts, clusterVals = clusterVals, factorTable = patients.cohortTable, colorKey = patients.cohortColorKey)
        
        ctClust <- subset(ctClust, select = -patients.cohort)
      }
      
      ## marker & cohort
      if("marker.cohort" %in% fisherTests){
        ctClust$marker.cohort <- paste(ctClust$marker, ctClust$cohort, sep = ".")
        
        cohorts <- unique(ctClust$cohort)
        cohorts <- cohorts[match(c("Child","Adult","Risk","Healthy"), cohorts, nomatch=F)]
        markers <- unique(ctClust$marker)
        markers <- markers[order(as.numeric(markers))]
        
        eg <- expand.grid(markers, cohorts)
        markers.cohorts <- sprintf('%s.%s', eg[,1], eg[,2])
        
        marker.cohortTable <- matrix(nrow=length(clusterVals), ncol=length(markers.cohorts))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(markers.cohorts)){
            marker.cohortTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                         marker.cohort==markers.cohorts[j]))
          }
        }
        colnames(marker.cohortTable) <- markers.cohorts
        rownames(marker.cohortTable) <- paste("cluster_", clusterVals, sep="")
        marker.cohortTable <- marker.cohortTable[, which(colSums(marker.cohortTable) > 0)]
        markers.cohorts <- colnames(marker.cohortTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Marker.Cohort vs. Cluster", quote=F)
        # print(patientsTable)
        print(marker.cohortTable, quote = F)
        # print(fisher.test(factor(ctClust$patients), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(marker.cohortTable, simulate.p.value=T))
        
        ## Bar Plots
        marker.cohortColors <- brewer.pal(12, "Paired")
        marker.cohortColorKey <- cbind(markers.cohorts, marker.cohortColors[1:length(markers.cohorts)])
        
        marker.cohortTablePlots <- tableBarPlots(plotFactor = "marker.cohort", plotFactorVals = markers.cohorts, clusterVals = clusterVals, factorTable = marker.cohortTable, colorKey = marker.cohortColorKey)
        
        grid.arrange(grobs = list(patients.markerTablePlots[[1]], patients.markerTablePlots[[2]], patients.cohortTablePlots[[1]], patients.cohortTablePlots[[2]], marker.cohortTablePlots[[1]], marker.cohortTablePlots[[2]]), nrow=3)
        
        ctClust <- subset(ctClust, select = -marker.cohort)
      }
      
      
      ## patients & cohort & marker
      if("patients.cohort.marker" %in% fisherTests){
        ctClust$patients.cohort.marker <- paste(ctClust$patients, ctClust$cohort, ctClust$marker, sep = ".")
      
        cohorts <- unique(ctClust$cohort)
        cohorts <- cohorts[match(c("Child","Adult","Risk","Healthy"), cohorts, nomatch=F)]
        patientss <- unique(ctClust$patients)
        patientss <- patientss[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY4","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), patientss, nomatch=F)]
        markers <- unique(ctClust$marker)
        markers <- markers[order(as.numeric(markers))]
        
        eg <- expand.grid(patientss, cohorts, markers)
        patientss.cohorts.markers <- sprintf('%s.%s.%s', eg[,1], eg[,2], eg[,3])
        
        patients.cohort.markerTable <- matrix(nrow=length(clusterVals), ncol=length(patientss.cohorts.markers))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(patientss.cohorts.markers)){
            patients.cohort.markerTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                  patients.cohort.marker==patientss.cohorts.markers[j]))
          }
        }
        colnames(patients.cohort.markerTable) <- patientss.cohorts.markers
        rownames(patients.cohort.markerTable) <- paste("cluster_", clusterVals, sep="")
        patients.cohort.markerTable <- patients.cohort.markerTable[, which(colSums(patients.cohort.markerTable) > 0)]
        patientss.cohorts.markers <- colnames(patients.cohort.markerTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient.Cohort.Marker vs. Cluster", quote=F)
        # print(patientsTable)
        print(patients.cohort.markerTable, quote = F)
        # print(fisher.test(factor(ctClust$patients), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(patients.cohort.markerTable, simulate.p.value=T))
        
        ## Bar Plots
        patients.cohort.markerColors <- brewer.pal(12, "Set3")
        patients.cohort.markerColorKey <- cbind(patientss.cohorts.markers, patients.cohort.markerColors[1:length(patientss.cohorts.markers)])
        
        patients.cohort.markerTablePlots <- tableBarPlots(plotFactor = "patients.cohort.marker", plotFactorVals = patientss.cohorts.markers, clusterVals = clusterVals, factorTable = patients.cohort.markerTable, colorKey = patients.cohort.markerColorKey)
        
        ng <- nullGrob()
        grid.arrange(grobs = list(patients.cohort.markerTablePlots[[1]], patients.cohort.markerTablePlots[[2]], ng), nrow=3)
        
        ctClust <- subset(ctClust, select = -patients.cohort.marker)
      }
      
      ## mouse
      # if("mouse" %in% fisherTests){
      #   mouses <- unique(ctClust$mouse)
      #   mouseTable <- matrix(nrow=length(clusterVals), ncol=length(mouses))
      #   for(i in 1:length(clusterVals)){
      #     for(j in 1:length(mouses)){
      #       mouseTable[i, j] <- nrow(subset(ctClust, mouse==mouses[j] & 
      #                                           kmeans.cluster==clusterVals[i]))
      #     }
      #   }
      #   colnames(mouseTable) <- paste("mouse_", mouses, sep="")
      #   rownames(mouseTable) <- paste("cluster_", clusterVals, sep="")
      #   print("     ", quote=F)
      #   print("     ", quote=F)
      #   print("Mouse vs. Cluster", quote=F)
      #   print(mouseTable)
      #   # print(fisher.test(factor(ctClust$mouse), factor(ctClust$kmeans.cluster), workspace = 1000000000))
      #   print(chisq.test(mouseTable, simulate.p.value=T))
      #   
      # } 
    }
#  }
  
  if(cumulativeExpHist == T){
    #### print histogram of cumulative expression for each cell, colored by cluster ####
    #ctClust$ctSum <- rowSums(ctClust[,10:ncol(ctClust)])
    ctExp <- (ctClust[,10:ncol(ctClust)])^2
    ctClust$ctSum <- log2(rowSums(ctExp))
    ctSumMin <- min(ctClust$ctSum) - (min(ctClust$ctSum) %% 50)
    ctSumMax <- 50 + max(ctClust$ctSum) - (max(ctClust$ctSum) %% 50)
    ctClust <- ctClust[order(ctClust$kmeans.cluster),]
    
    hist <- ggplot(ctClust, aes(ctSum, fill=factor(kmeans.cluster))) + 
      geom_histogram(bins = 150) +
      scale_x_continuous(breaks = seq(ctSumMin, ctSumMax, 50), 
                         minor_breaks = seq(ctSumMin, ctSumMax, 10)) +
      scale_fill_brewer(palette = "Set3") +
      guides(fill=guide_legend(title="clusters")) +
      ggtitle("histogram of cumulative expression values per cell") +
      xlab("sum of normalized expression values") +
      theme(text=element_text(size=25),
            axis.title.x=element_text(vjust=-0.5),
            plot.margin=unit(c(8,2,8,2),"cm"))
    print(hist)
    
    ## for insulin paper
    # date <- gsub("-", "_", Sys.Date())
    # ggsave(paste("cumulative_exp_hist_", date, ".png", sep=""), hist, device = "png", width = 14, height = 13, path = "~/Documents/abe/biomark/qpcr/figures/insulin/results")
    
    ctClust$ctSum <- NULL
  }
    
    ## filter out low expression clusters
  if(filterClusters == T){
    totalCells <- nrow(ctClust)
    
    ctClust <- subset(ctClust, !(kmeans.cluster %in% clustersToRemove))
    
    remCells <- totalCells - nrow(ctClust)
    remClust <- paste(clustersToRemove, collapse = ", ")
    
    if(removalText == T){
      cat(paste(remCells, " / ", totalCells, " cells were removed due to low gene expression (cluster ID of removed clusters: ", remClust, ")", sep=""))
    }
  }
  
  return(ctClust)
  
}