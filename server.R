#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# LOAD LIBRARIES
library(shiny)


# LOAD DATA
## META
  load(file = "Inputs/probeMapByUniqueGenes.Rdata")
  GeneMappingGroups=read.csv("Inputs/GeneMappingGroups.csv")
  
  
## CLINICAL DATA
  load("Inputs/clinData.Rdata")
  
## REDUCED HM450 & CPM DATA
  load(file = "Inputs/ReducedCPM.Rdata")
  load(file = "Inputs/ReducedHM450_Body.Rdata")
  load(file = "Inputs/ReducedHM450_TSS.Rdata")
  
  
myTheme=function(titleFont=15,generalFont=10,legendFont=8,axisFont=8,subtitleFont=12){
  theme(axis.text = element_text(face = "bold",color = "black",size = axisFont),
              plot.title = element_text(face="bold",colour = "darkblue",size = titleFont),
              title = element_text(face="bold",colour = "darkblue",size = generalFont), 
              legend.text = element_text(face="bold",color="darkblue",size=legendFont),
              plot.subtitle =element_text(color="red",size = subtitleFont))
}

  
# FUNCTIONS
metaDataPlot=function(type){
  library(ggplot2)
  if (type == "1"){
    p=ggplot(GeneMappingGroups,aes(reorder(Var1,Freq), Freq/sum(GeneMappingGroups$Freq), fill=Var1))+
      geom_bar(stat = "identity",show.legend = F)+xlab("Probe-Gene Mapping Categories")+ylab("Proportion of Probes")+geom_text(aes(label=Freq))
      theme_classic()+myTheme(generalFont = 20,axisFont = 15)
    return(p)
  }else{
    library(ComplexHeatmap)
    library(circlize)
    mat=as.data.frame(probeMappingTest)
    rownames(mat)=probeMappingTest$Gene
    mat=mat[-1]
    mat=as.matrix(mat)
    
    probeCounts=mat[,1]
    mat=mat[,-1]
    colnames(mat)=c(colnames(mat)[-c(5,6,7)],"5 UTR","3 UTR","1st Exon")
    
    column_ha = rowAnnotation(`Region Proportions` =anno_barplot( round(colSums(mat)/nrow(mat),2),gp = gpar(fill="brown")),width = unit(2,"cm"))
    row_ha = HeatmapAnnotation(`Probes per gene`= anno_barplot(probeCounts),height = unit(2,"cm"))
    marks = HeatmapAnnotation(foo = anno_mark(side = "left",at = which(probeCounts>300), 
                                          labels_gp = gpar(fontsize=10),labels = names(which(probeCounts>300))))
    
    ht=Heatmap(t(mat),col = colorRamp2(c(0,1), c("black", "yellow")),cluster_columns = F, cluster_rows = F,
            show_column_names = F,right_annotation = column_ha, top_annotation = row_ha,
            name = "Probe Proportion",column_title = "Genes",bottom_annotation = marks)
    return(ht)
  }
}


## GENE PLOTS
makeGenePlot=function(gene="CDH1", compareBy, DoColSplit=T,type="1"){
  print("here1")
  
  makeCompareByGroups=function(clinData, compareBy="ER Status"){
    if (compareBy == "ER Status"){
      return(subset(clinData,`ER IHC` %in% c("Positive","Negative","Adjacent Normal"))$Case.ID)
    }else if(compareBy == "Histology"){
      return(subset(clinData,`ER IHC` %in% c("Positive","Adjacent Normal") & `Final Pathology` %in% c("IDC","ILC","Adjacent Normal"))$Case.ID)
    }else if(compareBy == "PAM50"){
      return(subset(clinData,PAM50 %in% c("LumA","LumB","Her2","Basal","Normal","Adjacent Normal"))$Case.ID)
    }
  }
  
  makeColsplitGroups=function(compareBy="ER Status",sampleList){
    if (compareBy == "ER Status"){
      field="ER IHC"
    }else if(compareBy == "Histology"){
      field="Final Pathology"
    }else if(compareBy == "PAM50"){
      field="PAM50"
    }
    return(data.frame(factor(clinData[sampleList,field])))
  }
  
  makeHeatmapAnnotation=function(compareBy="ER Status",sampleList){
      if (compareBy == "ER Status"){
        annot=HeatmapAnnotation(ER = factor(clinData[sampleList,"ER IHC"]),na_col = "gray", height = unit(1, "cm"),
                                col = list(ER = c("Positive"="#FF6E40","Negative"="#1DE9B6","NA"="gray","Adjacent Normal"="#8E24AA")))
      }else if(compareBy == "Histology"){
        annot=HeatmapAnnotation(Histology = factor(clinData[sampleList,"Final Pathology"]),na_col = "gray", height = unit(1, "cm"),
                                col = list(Histology = c("IDC" = "#29B6F6", "ILC" = "orange", "IDC-L"="#4E342E", "Other"="#546E7A","Adjacent Normal"="#8E24AA")))
      }else if(compareBy == "PAM50"){
        annot=HeatmapAnnotation(PAM50 = factor(clinData[sampleList,"PAM50"]),na_col = "gray", height = unit(1, "cm"),
                                col = list(PAM50 = c("LumA"="#1A237E","LumB"="#29B6F6","Basal"="#FF6F00","Normal"="#64DD17","Her2"="#ff8a80","NA"="gray","Adjacent Normal"="#8E24AA")))
      }
      return(annot)
  }
  
  chooseClinAnnot=function(clinData,sampleList,compareBy){
    if (compareBy == "ER Status"){
      return(subset(clinData[sampleList,],`ER IHC` %in% c("Positive","Negative","Adjacent Normal"))[,c("Case.ID","ER IHC")])
    }else if(compareBy == "Histology"){
      return(subset(clinData[sampleList,],`ER IHC` %in% c("Positive","Adjacent Normal") & `Final Pathology` %in% c("IDC","ILC","Adjacent Normal"))[,c("Case.ID","Final Pathology")])
    }else if(compareBy == "PAM50"){
      return(subset(clinData[sampleList,],PAM50 %in% c("LumA","LumB","Her2","Basal","Normal","Adjacent Normal"))[,c("Case.ID","PAM50")])
    }
  }
  
  ggplotRegression <- function (fit) {
    labs(title=gene,subtitle = paste0("Adj R2 = ",round(summary(fit)$adj.r.squared, 2),
                                      ", Slope =",round(fit$coef[[2]], 2),
                                      ", P =",round(summary(fit)$coef[2,4], 3)))
  }
  
  
  
  # SAMPLE LIST
  sampleList=Reduce(intersect, list(makeCompareByGroups(clinData,compareBy),colnames(TCGA_BRCA_CPM_Reduced),colnames(HM450_TSS_Avg_Reduced)))
  
  print("here1")
  
  matHM450=cbind(HM450_TSS_Avg_Reduced[gene,sampleList],HM450_Body_Avg_Reduced[gene,sampleList]); 
  matCPM=as.matrix(TCGA_BRCA_CPM_Reduced[gene,sampleList]);
  
  if (type == "1"){
    library(ComplexHeatmap)
    library(circlize)
    print("here2")
    
    mat=apply(cbind(matCPM,matHM450),2, FUN = function(x) scale(x,center = T, scale = T))
    colnames(mat)=paste0(c("Gene Expression","TSS Methylation","Gene Body Methylation"))
    
    column_ha=makeHeatmapAnnotation(compareBy,sampleList)
    
    colSplit=NULL
    if (DoColSplit){ colSplit=makeColsplitGroups(compareBy,sampleList) }
    
    ht=Heatmap(t(mat),row_names_side = "left",show_row_names = T, name = "Signal Intensity",column_title  = "Samples",
            top_annotation = column_ha,show_heatmap_legend = T,column_split = colSplit)
    
    print("here3")
    return(ht)
    }
  
  
  library(RColorBrewer)
  myColors <- brewer.pal(6,"Set1")
  
  genePlot=function(df,cols=NULL){
    p=ggplot(df,aes(df[,1],df[,2]))+geom_point(aes(color=df[,3]),alpha=0.5,size=4)+theme_classic()+
      guides(size=10, color=guide_legend(title = compareBy))+ geom_smooth(method = "lm", color="black",se = T,size=1)+
      ggplotRegression(lm(df[,1] ~ df[,2], data = df))
    if(length(cols)>0){ 
      p=p+scale_colour_manual(values = cols)
    }else(
      p=p+scale_colour_manual(values =  myColors[unique(df[,3])])
    )
    return(p)
  }
  
  
  
  if (type == "2"){
    annot=chooseClinAnnot(clinData,sampleList,compareBy)
    p2DF=data.frame(exp=matCPM[annot[,1],1], TSS=matHM450[annot[,1],1], Body=matHM450[annot[,1],2],annot=annot[,2])
    
    
    if (!DoColSplit){
      # TSS MAIN
      TSSplot=genePlot(p2DF[,c("TSS","exp","annot")])+xlab("TSS Methylation (B-Value)")+ylab("Expression (Counts Per Million)")+myTheme(generalFont = 15,axisFont = 10)
      return(TSSplot)
    }else{
      # TSS SUB
      TSSsubplots=lapply(unique(p2DF$annot), function(x){
          col=myColors[match(x,unique(p2DF$annot))]
          genePlot(p2DF[p2DF$annot %in% x,][,c("TSS","exp","annot")],cols = col)+xlab("TSS Methylation (B-Value)")+ylab("Expression (Counts Per Million)")+
            myTheme(generalFont = 8,subtitleFont = 8)
      })
      names(TSSsubplots)=unique(p2DF$annot)
      return(do.call(ggpubr::ggarrange,TSSsubplots))
    }
  }
  
  
  if (type == "3"){
    annot=chooseClinAnnot(clinData,sampleList,compareBy)
    p2DF=data.frame(exp=matCPM[annot[,1],1], TSS=matHM450[annot[,1],1], Body=matHM450[annot[,1],2],annot=annot[,2])
    
    if (!DoColSplit){
      # BODY MAIN
      Bodyplot=genePlot(p2DF[,c("Body","exp","annot")])+xlab("Body Methylation (B-Value)")+ylab("Expression (Counts Per Million)")+myTheme(generalFont = 15,axisFont = 10)
      return(Bodyplot)
    }else{
      # BODY SUB
    Bodysubplots=lapply(unique(p2DF$annot), function(x){
      col=myColors[match(x,unique(p2DF$annot))]
      genePlot(p2DF[p2DF$annot %in% x,][,c("Body","exp","annot")],cols = col)+xlab("Body Methylation (B-Value)")+ylab("Expression (Counts Per Million)")+
        myTheme(generalFont = 8,subtitleFont = 8)
    })
    names(Bodysubplots)=unique(p2DF$annot)  
    return(do.call(ggpubr::ggarrange,Bodysubplots))
    }
  }

}


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  # META DATA PLOT FUNCTION
  mdpwidth=reactive({
    return(ifelse(input$metaDataPlotType == "1",700,1300))
  })
  mdpheight=reactive({
    return(ifelse(input$metaDataPlotType == "1",700,700))
  })
  
  output$metaDataPlot <- renderPlot(res = 120,width =mdpwidth, height = mdpheight,{
    # metaDataPlot function
    metaDataPlot(type=input$metaDataPlotType)
  })
  
  
  # GENE PLOT FUNCTION
  gpwidth=reactive({
    return(ifelse(input$genePlotType == "1",1300,ifelse(input$DoColSplit, 1100, 800)))
  })
  gpheight=reactive({
    return(ifelse(input$genePlotType == "1",600,ifelse(input$DoColSplit, 700,600)))
  })

  output$genePlot <- renderPlot(res = 120,width =gpwidth, height = gpheight,{
    # genePlot function
    makeGenePlot(gene=input$geneID, compareBy=input$comparyBy, DoColSplit=input$DoColSplit, type=input$genePlotType)
  })

})
