#### cluster and filter to remove outliers

clusterFilter <- function(ctInput, testK = F, numCenters = 2, plotHeatmap=F, plotClustOnly=F, 
                          heatmapFactor, heatmapColorBy, heatmapTissueLabel,
                          # fisherTests = c("probe", "tissue", "age", "mouse"),
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
    ## label cells by age
    ## age -> markers
    if("age" %in% heatmapColorBy){
      ctClust$ageColor <- NA
      ctClust[which(ctClust$age=="12-20"), "ageColor"] <- "#9A8822"
      ctClust[which(ctClust$age=="13-21"), "ageColor"] <- "#F5CDB4"
      ctClust[which(ctClust$age=="CP11"), "ageColor"] <- "#F8AFA8"
      ctClust[which(ctClust$age=="CP13"), "ageColor"] <- "#FDDDA0"
      ctClust[which(ctClust$age=="CP18"), "ageColor"] <- "#74A089"
      ctClust[which(ctClust$age=="CXCR3+/PD1-"), "ageColor"] <- "#FF0000"
      ctClust[which(ctClust$age=="ICOS+/PD1-"), "ageColor"] <- "#556A5B"
      ctClust[which(ctClust$age=="PD1-/CXCR3-"), "ageColor"] <- "#50A45C"
      ctClust[which(ctClust$age=="PD1-/ICOS-"), "ageColor"] <- "#F2AD00"
      ctClust[which(ctClust$age=="PD1+/CXCR3+"), "ageColor"] <- "#F69100"
      ctClust[which(ctClust$age=="PD1+/ICOS-"), "ageColor"] <- "#C49647"
      ctClust[which(ctClust$age=="PD1+/ICOS+"), "ageColor"] <- "#5BBCD6"

      Ages <- ctClust$ageColor
      ctClust <- subset(ctClust, select=-c(ageColor))
      AnnoColors <- cbind(AnnoColors, Ages)
    }
    
    ## label cells by source: Child, Adult, Risk
    if("source" %in% heatmapColorBy){
      ctClust$sourceColor <- NA
      ctClust[which(ctClust$cellSource=="Child"), "sourceColor"] <- "#FDE725FF"
      ctClust[which(ctClust$cellSource=="Adult"), "sourceColor"] <- "#21908CFF"
      ctClust[which(ctClust$cellSource=="Risk"), "sourceColor"] <- "#440154FF"
      ctClust[which(ctClust$cellSource=="Healthy"), "sourceColor"] <- "gray48"
      Tissues <- ctClust$sourceColor
      ctClust <- subset(ctClust, select=-c(sourceColor))
      AnnoColors <- cbind(AnnoColors, Tissues)
    }
    
    ## label cells by probe
    ## probe -> patients
    if("probe" %in% heatmapColorBy){
      ctClust$probeColor <- NA
      ctClust[which(ctClust$probe=="Healthy"), "probeColor"] <- "gray19"
      ctClust[which(ctClust$probe=="RAD1"), "probeColor"] <- "#DD8D29"
      ctClust[which(ctClust$probe=="RAD2"), "probeColor"] <- "#E1C408"
      ctClust[which(ctClust$probe=="RAD3"), "probeColor"] <- "#84BB78"
      ctClust[which(ctClust$probe=="RAD4"), "probeColor"] <- "#859C78"
      ctClust[which(ctClust$probe=="RAD5"), "probeColor"] <- "#DB6E07"
      ctClust[which(ctClust$probe=="RAD6"), "probeColor"] <- "#B40F20"
      ctClust[which(ctClust$probe=="TEY14"), "probeColor"] <- "#798E87"
      ctClust[which(ctClust$probe=="TEY3"), "probeColor"] <- "#C27D38"
      ctClust[which(ctClust$probe=="TEY6"), "probeColor"] <- "#CCC591"
      ctClust[which(ctClust$probe=="TEY8"), "probeColor"] <- "#29211F"
      ctClust[which(ctClust$probe=="TNET1"), "probeColor"] <- "#F3DF6C"
      ctClust[which(ctClust$probe=="TNET2"), "probeColor"] <- "#CEAB07"
      ctClust[which(ctClust$probe=="TNET3"), "probeColor"] <- "#D5D5D3"
      ctClust[which(ctClust$probe=="TNET4"), "probeColor"] <- "#24281A"
      Probes <- ctClust$probeColor
      ctClust <- subset(ctClust, select=-c(probeColor))
      AnnoColors <- cbind(AnnoColors, Probes)
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
    
    ageColorKey <- rbind(c("12-20", "#9A8822"),
                         c("13-21", "#F5CDB4"),
                         c("CP11", "#F8AFA8"),
                         c("CP13", "#FDDDA0"),
                         c("CP18", "#74A089"),
                         c("CXCR3+/PD1-", "#FF0000"),
                         c("ICOS+/PD1-", "#556A5B"),
                         c("PD1-/CXCR3-", "#50A45C"),
                         c("PD1-/ICOS-", "#F2AD00"),
                         c("PD1+/CXCR3+", "#F69100"),
                         c("PD1+/ICOS-", "#C49647"),
                         c("PD1+/ICOS+", "#5BBCD6"))
    ageColorKey <- ageColorKey[which(ageColorKey[,1] %in% ctClust$age),]
    
    sourceColorKey <- rbind(c("Child", "#FDE725FF"),
                            c("Adult","#21908CFF"),
                            c("Risk","#440154FF"),
                            c("Healthy", "gray48"))
    sourceColorKey <- sourceColorKey[which(sourceColorKey[,1] %in% ctClust$cellSource),]
    
    probeColorKey <- rbind(c("Healthy", "gray19"),
                           c("RAD1", "#DD8D29"),
                           c("RAD2", "#E1C408"),
                           c("RAD3", "#84BB78"),
                           c("RAD4", "#859C78"),
                           c("RAD5", "#DB6E07"),
                           c("RAD6", "#B40F20"),
                           c("TEY14", "#798E87"),
                           c("TEY3", "#C27D38"),
                           c("TEY6", "#CCC591"),
                           c("TEY8", "#29211F"),
                           c("TNET1", "#F3DF6C"),
                           c("TNET2", "#CEAB07"),
                           c("TNET3", "#D5D5D3"),
                           c("TNET4", "#24281A"))
    probeColorKey <- probeColorKey[which(probeColorKey[,1] %in% ctClust$probe),]
    
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
#               main=paste("heatmap for", heatmapTissueLabel, sep=" "), 
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
                main=paste("heatmap for", heatmapTissueLabel, sep=" "), 
                xlab="log2(gene expression / Gapdh expression)", 
                ylab="Individual Cells\n\n\n")
      if("age" %in% heatmapColorBy){
        legend("topleft",
               legend=ageColorKey[,1],
               fill=ageColorKey[,2], 
               title="age",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
        insetX <- insetX + 0.12
      }
      if("source" %in% heatmapColorBy){
        legend("topleft",
               legend=sourceColorKey[,1],
               fill=sourceColorKey[,2], 
               title="tissue",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.05))
        insetX <- insetX + 0.12
      }
      if("probe" %in% heatmapColorBy){
        legend("topleft",
               legend=probeColorKey[,1],
               fill=probeColorKey[,2], 
               title="probe",
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
                  main=paste("heatmap for", heatmapTissueLabel, sep=" "), 
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
      
      ## probe
      if("probe" %in% fisherTests){
        probes <- unique(ctClust$probe)
        probes <- probes[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), probes, nomatch=F)]
        probeTable <- matrix(nrow=length(clusterVals), ncol=length(probes))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes)){
            probeTable[i, j] <- nrow(subset(ctClust, probe==probes[j] & 
                                                      kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(probeTable) <- paste("probe_", probes, sep="")
        rownames(probeTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient vs. Cluster", quote=F)
        # print(probeTable)
        print(probeTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probeTable, simulate.p.value=T))
        
        ## Bar Plots
        probeTablePlots <- tableBarPlots(plotFactor = "probe", plotFactorVals = probes, clusterVals = clusterVals, factorTable = probeTable, colorKey = probeColorKey)
      }
      
      
      ## patient cohorts
      if("tissue" %in% fisherTests){
        cellSources <- unique(ctClust$cellSource)
        cellSources <- cellSources[match(c("Child","Adult","Risk","Healthy"), cellSources, nomatch=F)]
        sourceTable <- matrix(nrow=length(clusterVals), ncol=length(cellSources))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(cellSources)){
            sourceTable[i, j] <- nrow(subset(ctClust, cellSource==cellSources[j] & 
                                               kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(sourceTable) <- paste("cellSource_", cellSources, sep="")
        rownames(sourceTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Cohort vs. Cluster", quote=F)
        print(sourceTable)
        # print(fisher.test(factor(ctClust$cellSource), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(sourceTable, simulate.p.value=T))
        
        ## Bar Plots
        tissueTablePlots <- tableBarPlots(plotFactor = "tissue", plotFactorVals = cellSources, clusterVals = clusterVals, factorTable = sourceTable, colorKey = sourceColorKey)
      }
      
      ## age
      if("age" %in% fisherTests){
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        ageTable <- matrix(nrow=length(clusterVals), ncol=length(ages))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(ages)){
            ageTable[i, j] <- nrow(subset(ctClust, age==ages[j] & 
                                                kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(ageTable) <- paste("age_", ages, sep="")
        rownames(ageTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Marker vs. Cluster", quote=F)
        print(ageTable)
        # print(fisher.test(factor(ctClust$age), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(ageTable, simulate.p.value=T))
        
        ## Bar Plots
        ageTablePlots <- tableBarPlots(plotFactor = "age", plotFactorVals = ages, clusterVals = clusterVals, factorTable = ageTable, colorKey = ageColorKey)
        
#         grid.arrange(grobs = list(probeTablePlots[[3]], probeTablePlots[[4]], tissueTablePlots[[3]], tissueTablePlots[[4]], ageTablePlots[[3]], ageTablePlots[[4]]), nrow=3)
#         
#         grid.arrange(grobs = list(probeTablePlots[[1]], probeTablePlots[[2]], tissueTablePlots[[1]], tissueTablePlots[[2]], ageTablePlots[[1]], ageTablePlots[[2]]), nrow=3)
      }
      ng <- nullGrob()
      if("probe" %in% fisherTests){}
#       if("age" %in% fisherTests){}
#       if("tissue" %in% fisherTests){}
      
      if("probe" %in% fisherTests & "tissue" %in% fisherTests & "age" %in% fisherTests){
        grid.arrange(grobs = list(probeTablePlots[[3]], probeTablePlots[[4]], tissueTablePlots[[3]], tissueTablePlots[[4]], ageTablePlots[[3]], ageTablePlots[[4]]), nrow=3)
        
        grid.arrange(grobs = list(probeTablePlots[[1]], probeTablePlots[[2]], tissueTablePlots[[1]], tissueTablePlots[[2]], ageTablePlots[[1]], ageTablePlots[[2]]), nrow=3)
      }
      
      if("probe" %in% fisherTests & "age" %in% fisherTests){
        grid.arrange(grobs = list(probeTablePlots[[3]], probeTablePlots[[4]], ageTablePlots[[3]], ageTablePlots[[4]], ng, ng), nrow=3)
        
        grid.arrange(grobs = list(probeTablePlots[[1]], probeTablePlots[[2]], ageTablePlots[[1]], ageTablePlots[[2]], ng, ng), nrow=3)
      }
      
      ## probe and age
      if("probe.age" %in% fisherTests){
        ctClust$probe.age <- paste(ctClust$probe, ctClust$age, sep = ".")
        probes <- unique(ctClust$probe)
        probes <- probes[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), probes, nomatch=F)]
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        eg <- expand.grid(probes, ages)
        probes.ages <- sprintf('%s.%s', eg[,1], eg[,2])
        
        probe.ageTable <- matrix(nrow=length(clusterVals), ncol=length(probes.ages))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes.ages)){
            probe.ageTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                          probe.age==probes.ages[j]))
          }
        }
        colnames(probe.ageTable) <- probes.ages
        rownames(probe.ageTable) <- paste("cluster_", clusterVals, sep="")
        probe.ageTable <- probe.ageTable[, which(colSums(probe.ageTable) > 0)]
        probes.ages <- colnames(probe.ageTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient.Marker vs. Cluster", quote=F)
        # print(probeTable)
        print(probe.ageTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probe.ageTable, simulate.p.value=T))
        
        ## Bar Plots
        probe.ageColors <- brewer.pal(12, "Paired")
        probe.ageColorKey <- cbind(probes.ages, probe.ageColors[1:length(probes.ages)])
        
        probe.ageTablePlots <- tableBarPlots(plotFactor = "probe.age", plotFactorVals = probes.ages, clusterVals = clusterVals, factorTable = probe.ageTable, colorKey = probe.ageColorKey)
        
        grid.arrange(grobs = list(probe.ageTablePlots[[1]], probe.ageTablePlots[[2]],  ng, ng, ng, ng), nrow=3)
        
        ctClust <- subset(ctClust, select = -probe.age)
      }
      
      ## probe and tissue
      if("probe.tissue" %in% fisherTests){
        ctClust$probe.tissue <- paste(ctClust$probe, ctClust$cellSource, sep = ".")
        probes <- unique(ctClust$probe)
        probes <- probes[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), probes, nomatch=F)]
        tissues <- unique(ctClust$cellSource)
        tissues <- tissues[match(c("Child","Adult","Risk","Healthy"), tissues, nomatch=F)]
        eg <- expand.grid(probes, tissues)
        probes.tissues <- sprintf('%s.%s', eg[,1], eg[,2])
        
        probe.tissueTable <- matrix(nrow=length(clusterVals), ncol=length(probes.tissues))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes.tissues)){
            probe.tissueTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                     probe.tissue==probes.tissues[j]))
          }
        }
        colnames(probe.tissueTable) <- probes.tissues
        rownames(probe.tissueTable) <- paste("cluster_", clusterVals, sep="")
        probe.tissueTable <- probe.tissueTable[, which(colSums(probe.tissueTable) > 0)]
        probes.tissues <- colnames(probe.tissueTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient.Cohort vs. Cluster", quote=F)
        # print(probeTable)
        print(probe.tissueTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probe.tissueTable, simulate.p.value=T))
        
        ## Bar Plots
        probe.tissueColors <- brewer.pal(12, "Paired")
        probe.tissueColorKey <- cbind(probes.tissues, probe.tissueColors[1:length(probes.tissues)])
        
        probe.tissueTablePlots <- tableBarPlots(plotFactor = "probe.tissue", plotFactorVals = probes.tissues, clusterVals = clusterVals, factorTable = probe.tissueTable, colorKey = probe.tissueColorKey)
        
        ctClust <- subset(ctClust, select = -probe.tissue)
      }
      
      ## age & tissue
      if("age.tissue" %in% fisherTests){
        ctClust$age.tissue <- paste(ctClust$age, ctClust$cellSource, sep = ".")
        
        cellSources <- unique(ctClust$cellSource)
        cellSources <- cellSources[match(c("Child","Adult","Risk","Healthy"), cellSources, nomatch=F)]
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        
        eg <- expand.grid(ages, cellSources)
        ages.tissues <- sprintf('%s.%s', eg[,1], eg[,2])
        
        age.tissueTable <- matrix(nrow=length(clusterVals), ncol=length(ages.tissues))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(ages.tissues)){
            age.tissueTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                         age.tissue==ages.tissues[j]))
          }
        }
        colnames(age.tissueTable) <- ages.tissues
        rownames(age.tissueTable) <- paste("cluster_", clusterVals, sep="")
        age.tissueTable <- age.tissueTable[, which(colSums(age.tissueTable) > 0)]
        ages.tissues <- colnames(age.tissueTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Marker.Cohort vs. Cluster", quote=F)
        # print(probeTable)
        print(age.tissueTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(age.tissueTable, simulate.p.value=T))
        
        ## Bar Plots
        age.tissueColors <- brewer.pal(12, "Paired")
        age.tissueColorKey <- cbind(ages.tissues, age.tissueColors[1:length(ages.tissues)])
        
        age.tissueTablePlots <- tableBarPlots(plotFactor = "age.tissue", plotFactorVals = ages.tissues, clusterVals = clusterVals, factorTable = age.tissueTable, colorKey = age.tissueColorKey)
        
        grid.arrange(grobs = list(probe.ageTablePlots[[1]], probe.ageTablePlots[[2]], probe.tissueTablePlots[[1]], probe.tissueTablePlots[[2]], age.tissueTablePlots[[1]], age.tissueTablePlots[[2]]), nrow=3)
        
        ctClust <- subset(ctClust, select = -age.tissue)
      }
      
      
      ## probe & tissue & age
      if("probe.tissue.age" %in% fisherTests){
        ctClust$probe.tissue.age <- paste(ctClust$probe, ctClust$cellSource, ctClust$age, sep = ".")
      
        cellSources <- unique(ctClust$cellSource)
        cellSources <- cellSources[match(c("Child","Adult","Risk","Healthy"), cellSources, nomatch=F)]
        probes <- unique(ctClust$probe)
        probes <- probes[match(c("Healthy","RAD1","RAD2","RAD3","RAD4","RAD5","RAD6","TEY3","TEY6","TEY8","TEY14","TNET1","TNET2","TNET3","TNET4"), probes, nomatch=F)]
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        
        eg <- expand.grid(probes, cellSources, ages)
        probes.tissues.ages <- sprintf('%s.%s.%s', eg[,1], eg[,2], eg[,3])
        
        probe.tissue.ageTable <- matrix(nrow=length(clusterVals), ncol=length(probes.tissues.ages))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes.tissues.ages)){
            probe.tissue.ageTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                  probe.tissue.age==probes.tissues.ages[j]))
          }
        }
        colnames(probe.tissue.ageTable) <- probes.tissues.ages
        rownames(probe.tissue.ageTable) <- paste("cluster_", clusterVals, sep="")
        probe.tissue.ageTable <- probe.tissue.ageTable[, which(colSums(probe.tissue.ageTable) > 0)]
        probes.tissues.ages <- colnames(probe.tissue.ageTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Patient.Cohort.Marker vs. Cluster", quote=F)
        # print(probeTable)
        print(probe.tissue.ageTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probe.tissue.ageTable, simulate.p.value=T))
        
        ## Bar Plots
        probe.tissue.ageColors <- brewer.pal(12, "Set3")
        probe.tissue.ageColorKey <- cbind(probes.tissues.ages, probe.tissue.ageColors[1:length(probes.tissues.ages)])
        
        probe.tissue.ageTablePlots <- tableBarPlots(plotFactor = "probe.tissue.age", plotFactorVals = probes.tissues.ages, clusterVals = clusterVals, factorTable = probe.tissue.ageTable, colorKey = probe.tissue.ageColorKey)
        
        ng <- nullGrob()
        grid.arrange(grobs = list(probe.tissue.ageTablePlots[[1]], probe.tissue.ageTablePlots[[2]], ng), nrow=3)
        
        ctClust <- subset(ctClust, select = -probe.tissue.age)
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