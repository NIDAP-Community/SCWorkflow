#' @title Filter & QC Samples 
#' @description Filters cells and Genes for each sample and generates QC Plots 
#' to evaluate data before and after filtering. 
#' @details This is Step 2 in the basic Single-Cell RNA-seq workflow. Multiple 
#' cell and gene filters can be selected to remove poor quality data and noise. 
#' Workflows can use this downstream of any Seurat Object. This tool is 
#' typically the second step in the Single Cell Workflow.
#' 
#' @param object a list of seurat objects for each sample.
#' @param min.cells Filter out genes found in less than this number of cells.
#'  E.g. Setting to 20 will remove genes found in fewer than 3 cells of 
#'  a sample. (Default: 20)
#' @param filter.vdj.genes If FALSE to remove VDJ genes from the scRNA
#'  transcriptome assay. This is to prevent clustering bias in T-cells of the 
#'  same clonotype. Only recommended if you are also doing TCR-seq. 
#'  (Default: FALSE)
#' @param nfeature.limits Filter out cells where the number of genes found in each 
#'  cell exceed the selected lower or upper limits. 
#'  Usage c(lower limit, Upper Limit). E.g. setting to c(200,1000) will remove 
#'  cells that have fewer than 200 genes or more than 1000 genes 
#'  for each sample. (Default: c(NA, NA))
#' @param mad.nfeature.limits Set filter limits based on how many Median 
#' Absolute Deviations an outlier cell will have. Calculated from the median 
#' gene number for all cells in your sample. Usage c(lower limit, Upper Limit)
#'  E.g. setting to c(3,5) will remove all cells with more than 3 absolute 
#'  deviations less than the median or 5 absolute deviations greater than the 
#'  median. (Default: c(5,5))
#' @param ncounts.limits Filter out cells where the total number of molecules 
#'  (umi) detected within a cell exceed the selected limits.
#'  Usage c(lower limit, Upper Limit). E.g. setting to c(200,100000) will remove 
#'  cells that have fewer than 200 or greater than 100000 molecules. 
#'  (Default: c(NA, NA))
#' @param mad.ncounts.limits Set filter limits based on how many Median Absolute 
#'  Deviations an outlier cell will have. Calculated from the median number of 
#'  molecules for all cells in your sample. Usage c(lower limit, Upper Limit)
#'  E.g. setting to c(3,5) will remove all cells with more than 3 absolute 
#'  deviations less than the median or with more than 5 absolute deviations 
#'  greater than the median. (Default: c(5,5))
#' @param mitoch.limits Filter out cells whose proportion of mitochondrial genes
#'  exceed the selected lower or upper limits. 
#'  Usage c(lower limit, Upper Limit). E.g. setting to c(0,8) will not set the 
#'  lower limit and removes cells with more than 8% mitochondrial RNA. 
#'  (Default: c(NA,25))
#' @param mad.mitoch.limits Set filter limits based on how many Median Absolute 
#'  Deviations an outlier cell will have. Calculated from the Median percentage 
#'  of mitochondrial RNA for all cells in your sample. 
#'  Usage c(lower limit, Upper Limit). E.g. setting to c(NA,3) will not set a 
#'  lower limit and remove all cells with more than 3 absolute deviations 
#'  greater than the median. (Default: c(NA,3))
#' @param complexity.limits Complexity represents Number of genes detected per 
#' UMI. The more genes detected per UMI, the more complex the data. 
#' Filter out cells whose Complexity exceed the selected lower or upper limits.
#' Cells that have a high number of UMIs but only a low number of genes could 
#' be dying cells, but also could represent a population of a low complexity 
#' cell type (i.e red blood cells). We suggest that you set the lower limit to 
#' 0.8 if samples have suspected RBC contamination. 
#' Usage c(lower limit, Upper Limit). E.g. setting to c(0.8,0) will not set 
#' an upper limit and removes cells with complexity less than 0.8.
#' (Default: c(NA,NA))
#' @param mad.complexity.limits Set filter limits based on how many Median 
#' Absolute Deviations an outlier cell will have. Calculated from the Median 
#' complexity for all cells in your sample. 
#' Usage c(lower limit, Upper Limit). E.g. setting to c(5,NA) will not set an 
#' upper limit and remove all cells with  more than 5 absolute deviations 
#' less than the median. (Default: c(5,NA))
#' @param topNgenes.limits Filter Cells based on the percentage of total counts 
#' in top N most highly expressed genes. Outlier cells would have a high 
#' percentage of counts in just a few genes and should be removed. 
#' The same considerations outlined in "complexity.limits" should be taken for 
#' this filter. 
#' Usage c(lower limit, Upper Limit). E.g. setting to c(NA,50) will not set a 
#' lower limit and remove cells with greater than 50% of reads in the top N 
#' genes. (Default: c(NA,NA))
#' @param mad.topNgenes.limitsSet Filter limits based on how many Median 
#' Absolute Deviations an outlier cell will have. Calculated from the Median 
#' percentage of counts in the top N Genes.
#' Usage c(lower limit, Upper Limit). E.g. setting to c(5,5) will remove all 
#' cells with more than 5 absolute deviations greater than or 5 absolute 
#' deviations less than the median percentage. (Default: c(5,5))
#' @param n.topgnes Select the number of top highly expressed genes used to 
#' calculate the percentage of reads found in these genes. 
#' E.g. a value of 20 calculates the percentage of reads found in the top 20 
#' most highly expressed Genes.
#' (Default: 20)
#' @param do.doublets.fitler Use scDblFinder to identify and remove doublet 
#' cells. Doublets are defined as two cells that are sequenced under the same 
#' cellular barcode, for example, if they were captured in the same droplet.
#' (Default: TRUE)
#' 
#' 
#' @importFrom Seurat CreateAssayObject Idents as.SingleCellExperiment 
#' @importFrom Seurat AddMetaData SCTransform FindVariableFeatures
#' @importFrom Seurat RunPCA RunUMAP
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange rename
#' @importFrom ggplot2 ggplot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr annotate_figure get_legend ggarrange
#' @importFrom gridExtra arrangeGrob tableGrob
#' @importFrom scDblFinder scDblFinder
#' @importFrom stringr str_split_fixed
#' @importFrom stats mad median
#' @importFrom grid grobHeight textGrob grid.newpage gTree grid.draw

#' 
#' @export
#' 
#' @return Seurat Object and QC plots

filterQC <- function(object,
                     
                     ## Filter Samples
                     min.cells = 20,
                     filter.vdj.genes=F,
                     nfeature.limits=c(NA,NA),
                     mad.nfeature.limits=c(5,5),
                     ncounts.limits=c(NA,NA),
                     mad.ncounts.limits=c(5,5),
                     mitoch.limits = c(NA,25),
                     mad.mitoch.limits = c(NA,3),
                     complexity.limits = c(NA,NA),
                     mad.complexity.limits = c(5,NA),
                     topNgenes.limits = c(NA,NA),
                     mad.topNgenes.limits = c(5,5),
                     n.topgnes=20,
                     do.doublets.fitler=T,
                     
                     ## dim Reduction settings
                     plot.outliers="None", #options(None,UMAP,tSNE) 
                     group.column = NA,
                     nfeatures = 2000,
                     low.cut = 0.1,
                     high.cut = 8,
                     low.cut.disp = 1,
                     high.cut.disp = 100000,
                     selection.method = "vst",
                     npcs = 30,
                     vars_to_regress=NULL,
                     seed.for.PCA = 42,
                     seed.for.TSNE = 1,
                     seed.for.UMAP = 42
                     
                     
                     
){
  
  
  
  ## --------- ##
  ## Functions ####
  ## --------- ##
  
  ### Helper Functions #####
  
  .topNGenes <- function(so,n.topgnes) { 
    ##Extract counts table
    counts_matrix = GetAssayData(so, slot="counts")
    
    ## calculate Counts in Top n genes 
    tbl=  apply(counts_matrix,2,function(i){
      cnts=i[order(i,decreasing=T)]
      
      t20=sum(cnts[1:n.topgnes])
      total=sum(cnts)
      
      pertop20=(t20/total)*100
      return(pertop20)
    })
    
    ### add to metadata  
    tbl=as.data.frame(tbl)
    so_out=AddMetaData(so, tbl, col.name = 'pct_counts_in_top_N_genes')
    
    return(so_out)  
  }
  
  .rowMaxs=function(mtx){
    apply(mtx, 1, max,na.rm=T)
  }
  .rowMins=function(mtx){
    apply(mtx, 1, min,na.rm=T)
  }
  
  .madCalc <- function(so,column,limits){
    stdev <- mad(so@meta.data[,column])
    med <- median(so@meta.data[,column])
    
    minlim <- med-(limits[1]*stdev)
    maxlim <- med+(limits[2]*stdev)
    # gl <- format(round(maxlim,0),nsmall=0)## Remove used in testing
    
    return(c(minlim,maxlim))
    
  }
  
  .checkLimits <- function(limits){
    minlim=limits[1]
    maxlim=limits[2]
    if(is.numeric(minlim)==F |is.na(minlim) |is.null(minlim)){minlim=-Inf}
    if(is.numeric(maxlim)==F|is.na(maxlim) |is.null(maxlim)){maxlim=Inf}
    return(c(minlim,maxlim))
  }  
  
  
  
  ### Plotting Functions ####
  
  
  #### Filter QC plots  ####
  
  .plotViolin2=function(count.df,value){
    axis.lab = unique(count.df$filt)
    ylabs=gsub(" \\(", "\n\\(",value)
    ylabs=gsub(paste0(" Top",n.topgnes), paste0("\nTop",n.topgnes),ylabs)
    
    ### Set up table fore cut off lines
    ## clean up cutoff values
    for (v in c(value)) {
      
      count.df[,paste0(v,'_Filters')]=gsub("Low:|High:","",
                                           count.df[,paste0(v,'_Filters')])
      count.df[,paste0('MAD ',v)]=gsub("Low:|High:","",
                                       count.df[,paste0('MAD ',v)])
      
      count.df[,paste0(v,'_Filters')]=gsub("\n",",",
                                           count.df[,paste0(v,'_Filters')])
      count.df[,paste0('MAD ',v)]=gsub("\n",",",
                                       count.df[,paste0('MAD ',v)])
      
      count.df[,paste0(v,'_Filters')]=gsub("-Inf|Inf","NA",
                                           count.df[,paste0(v,'_Filters')])
      count.df[,paste0('MAD ',v)]=gsub("-Inf|Inf","NA",
                                       count.df[,paste0('MAD ',v)])
    }
    
    ####### Y axis
    ## get sample cut off values
    count.df.lim=count.df[,c('Sample',paste0(value,'_Filters'),
                             paste0('MAD ',value))]%>%unique
    
    ## convert to numeric values
    count.df.lim_MAD=str_split_fixed(count.df.lim[,paste0('MAD ',value)],
                                     pattern=',',2)%>%as.data.frame
    rownames(count.df.lim_MAD)=count.df.lim$Sample
    colnames(count.df.lim_MAD)=c('Low','High')
    count.df.lim_MAD$Low=as.numeric(count.df.lim_MAD$Low)
    count.df.lim_MAD$High=as.numeric(count.df.lim_MAD$High)
    count.df.lim_cut=str_split_fixed(count.df.lim[,paste0(value,'_Filters')],
                                     pattern=',',2)%>%as.data.frame
    rownames(count.df.lim_cut)=count.df.lim$Sample
    colnames(count.df.lim_cut)=c('Low','High') 
    count.df.lim_cut$Low=as.numeric(count.df.lim_cut$Low)%>%suppressWarnings()
    count.df.lim_cut$High=as.numeric(count.df.lim_cut$High)%>%suppressWarnings()
    count.df.lim$Low=.rowMaxs(cbind(count.df.lim_cut[,'Low'],
                                   count.df.lim_MAD[,'Low']))#,na.rm=T)
    count.df.lim[count.df.lim$Low<0,'Low']=0
    count.df.lim$High=.rowMins(cbind(count.df.lim_cut[,'High'],
                                    count.df.lim_MAD[,'High']))#,na.rm=T)
    
    
    
    g <- ggplot(count.df, aes(x=filt, y=.data[[value]])) +
      # ggtitle(paste(name,count.df$variable[1])) +
      theme_classic()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=rel(1.5)),
            legend.title=element_blank(), 
            axis.text=element_text(size=10),
            axis.title.x=element_blank(),
            # axis.title.y=element_blank(),
            # axis.text.x=element_text(angle=45,hjust=1),
            axis.text.x=element_blank(),
            plot.title = element_text(size = 12, face = "bold")) +
      geom_violin(aes(fill=filt)) +  
      # scale_fill_manual(values = c("#00AFBB", "#FC4E07")) + 
      geom_boxplot(width=0) +
      # scale_x_discrete(limits = as.vector(axis.lab))+
      labs( y = ylabs)+
      geom_hline(data = count.df.lim, aes(yintercept = Low,
                                          linetype='Filter\nLimits'), 
                 col = 'grey')+
      geom_hline(data = count.df.lim, aes(yintercept = High),
                 linetype='dashed', col = 'grey')+
      scale_linetype_manual(name = "Filter\nLimits", values = 'dashed',
                            guide = guide_legend(override.aes = 
                                                   list(color = c("grey"),
                                                        size=2)))+
      facet_wrap(~Sample,nrow=1)%>%suppressWarnings()
  
    return(g)
  }
  
  
  
  .plotScatter2=function(count.df,value){
    # count.df$filt=factor(count.df$filt,levels = c('filt','raw'))
    count.df$filt=factor(count.df$filt,levels = c('raw','filt'))
    ylabs=gsub(" \\(", "\n\\(",value)
    ylabs=gsub(paste0(" Top",n.topgnes), paste0("\nTop",n.topgnes),ylabs)
    
    ### Set up table fore cut off lines
    ## clean up cutoff values
    for (v in c('UMI Count (nCount_RNA)',value)) {
      
    count.df[,paste0(v,'_Filters')]=gsub("Low:|High:","",
                                         count.df[,paste0(v,'_Filters')])
    count.df[,paste0('MAD ',v)]=gsub("Low:|High:","",
                                     count.df[,paste0('MAD ',v)])
    
    count.df[,paste0(v,'_Filters')]=gsub("\n",",",
                                         count.df[,paste0(v,'_Filters')])
    count.df[,paste0('MAD ',v)]=gsub("\n",",",count.df[,paste0('MAD ',v)])
    
    count.df[,paste0(v,'_Filters')]=gsub("-Inf|Inf","NA",
                                         count.df[,paste0(v,'_Filters')])
    count.df[,paste0('MAD ',v)]=gsub("-Inf|Inf","NA",
                                     count.df[,paste0('MAD ',v)])
    }
    
####### Y axis
    ## get sample cut off values
    count.df.lim=count.df[,c('Sample',paste0(value,'_Filters'),
                             paste0('MAD ',value))]%>%unique
    
    ## convert to numeric values
    count.df.lim_MAD=str_split_fixed(count.df.lim[,paste0('MAD ',value)],
                                     pattern=',',2)%>%as.data.frame
      rownames(count.df.lim_MAD)=count.df.lim$Sample
      colnames(count.df.lim_MAD)=c('Low','High')
      count.df.lim_MAD$Low=as.numeric(count.df.lim_MAD$Low)
      count.df.lim_MAD$High=as.numeric(count.df.lim_MAD$High)
    count.df.lim_cut=str_split_fixed(count.df.lim[,paste0(value,'_Filters')],
                                     pattern=',',2)%>%as.data.frame
      rownames(count.df.lim_cut)=count.df.lim$Sample
      colnames(count.df.lim_cut)=c('Low','High') 
      count.df.lim_cut$Low=as.numeric(count.df.lim_cut$Low)
      count.df.lim_cut$High=as.numeric(count.df.lim_cut$High)
    count.df.lim$Low=.rowMaxs(cbind(count.df.lim_cut[,'Low'],
                                   count.df.lim_MAD[,'Low']))#,na.rm=T)
    count.df.lim[count.df.lim$Low<0,'Low']=0
    count.df.lim$High=.rowMins(cbind(count.df.lim_cut[,'High'],
                                    count.df.lim_MAD[,'High']))#,na.rm=T)
    
####### X axis
    valueX='UMI Count (nCount_RNA)'
    ## get sample cut off values
    count.df.limX=count.df[,c('Sample',paste0(valueX,'_Filters'),
                              paste0('MAD ',valueX))]%>%unique
    
    ## convert to numeric values
    count.df.limX_MAD=str_split_fixed(count.df.limX[,paste0('MAD ',valueX)],
                                      pattern=',',2)%>%as.data.frame
      rownames(count.df.limX_MAD)=count.df.limX$Sample
      colnames(count.df.limX_MAD)=c('Low','High')
      count.df.limX_MAD$Low=as.numeric(count.df.limX_MAD$Low)
      count.df.limX_MAD$High=as.numeric(count.df.limX_MAD$High)
    count.df.limX_cut=str_split_fixed(count.df.limX[,paste0(valueX,'_Filters')],
                                      pattern=',',2)%>%as.data.frame
      rownames(count.df.limX_cut)=count.df.limX$Sample
      colnames(count.df.limX_cut)=c('Low','High') 
      count.df.limX_cut$Low=as.numeric(count.df.limX_cut$Low)
      count.df.limX_cut$High=as.numeric(count.df.limX_cut$High)
    count.df.limX$Low=.rowMaxs(cbind(count.df.limX_cut[,'Low'],
                                    count.df.limX_MAD[,'Low']))#,na.rm=T)
    count.df.limX[count.df.limX$Low<0,'Low']=0
    count.df.limX$High=.rowMins(cbind(count.df.limX_cut[,'High'],
                                     count.df.limX_MAD[,'High']))#,na.rm=T)
        
    
    
    g <- count.df%>%arrange(filt)%>%
      ggplot(aes(x=`UMI Count (nCount_RNA)`,y=.data[[value]],color=filt)) + 
      geom_point(size = 0.5) +
      
      theme_classic() +
      theme(
        # strip.background =element_rect(fill="grey"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=rel(1.5)),
            legend.title=element_blank(), 
            # axis.title.x=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            axis.text=element_text(size=10)
      ) +       
      guides(colour = guide_legend(override.aes = list(size=5)))+
      labs( y = ylabs)+
      geom_hline(data = count.df.lim, aes(yintercept = Low,
                                          linetype='Filter\nLimits'), 
                 col = 'grey')+
      geom_hline(data = count.df.lim, aes(yintercept = High),
                 linetype='dashed', col = 'grey')+
      geom_vline(data = count.df.limX, aes(xintercept = Low),linetype='dashed', 
                 col = 'grey')+
      geom_vline(data = count.df.limX, aes(xintercept = High),
                 linetype='dashed', col = 'grey')+
      scale_linetype_manual(name = "Filter Limits", values = 'dashed',
                            guide = guide_legend(override.aes = 
                                                   list(color = c("grey"),
                                                        size=2)))+
      facet_wrap(~Sample,nrow=1,)%>%suppressWarnings()
    
    return(g)
  }
  
  
  
  .combinePlots=function(plot.list){
    plot.list.mod=plot.list
    for (x in c(2:length(plot.list))) {
      plot.list.mod[[x]]=plot.list.mod[[x]]+theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
    }
      for (x in c(1:(length(plot.list)-1))) {
        plot.list.mod[[x]]=plot.list.mod[[x]]+theme(
          axis.title.x = element_blank()
        )
      }
    return(plot.list.mod)
  }
  
  
  .runTsnepPlot= function(filterCat,filterM,so,reduction){
    if (reduction=="umap") {
      qcFiltr.df.plot <- as.data.frame(so@reductions$umap@cell.embeddings)
    }else{
      qcFiltr.df.plot <- as.data.frame(so@reductions$tsne@cell.embeddings)
    }
    
    filterM=filterM[,filterCat,drop=F]
    
    qcFiltr.df.plot[,filterCat]=NA
    qcFiltr.df.plot[rownames(qcFiltr.df.plot)%in%
                   rownames(filterM[filterM[,1]==T,1,drop=F]),
                 colnames(filterM)]="Normal"
    qcFiltr.df.plot[rownames(qcFiltr.df.plot)%in%
                   rownames(filterM[filterM[,1]==F,1,drop=F]==F),
                 colnames(filterM)]="Removed"
    qcFiltr.df.plot[,colnames(filterM)]=
      factor(qcFiltr.df.plot[,colnames(filterM)],levels=c("Normal",'Removed'))
    
    g <- .plotTsne(qcFiltr.df.plot,so@project.name,filterCat)
    return(g) 
  }
  
  .plotTsne <- function(qcFiltr.df.plot,name,var){
    g <- ggplot(qcFiltr.df.plot[order(qcFiltr.df.plot[,var]),]) +
      geom_point(mapping = aes_string(x=colnames(qcFiltr.df.plot[,1,drop=F]),
                                      y=colnames(qcFiltr.df.plot[,2,drop=F]),
                                      color=var),
                 size = 1) + 
      theme_classic() + 
      ggtitle(paste(name,"\n",var)) +
      scale_color_manual(values = c("grey","red")) +
  #scale_colour_gradient2(midpoint = mid[i],low="blue",mid="grey",high="red") + 
      theme(legend.title=element_blank())%>%suppressWarnings()
  }
  
  #### Post Filter plots ####
  .plotScatterPost2=function(count.df,xaxis,yaxis){	
    xlab = as.character(xaxis)	
    ylab = as.character(yaxis)	
    ylab=gsub(" \\(", "\n\\(",ylab)
    ylab=gsub(paste0(" Top",n.topgnes), paste0("\nTop",n.topgnes),ylab)
    
    name = paste(ylab,"vs.",xlab)          
    g =ggplot(count.df, aes(x=.data[[xaxis]], y=.data[[yaxis]],color = Sample))+
      geom_point(size = 0.5) + 
      theme_classic() +
      theme(legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      guides(colour = guide_legend(override.aes = list(size=5))) +
      scale_color_manual(values = col2) +
      labs( x = xlab, y = ylab)%>%suppressWarnings()
    
    return(g)
  }
  
  .plotHistPost2 <- function(count.df,xaxis){	
    g=ggplot(count.df) + 
      theme_bw() +
      geom_density(aes(x = .data[[xaxis]], colour = Sample)) +
      # labs(x = NULL) +
      theme(legend.position='right',
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=12),
            legend.title=element_blank(), 
            # axis.title.x=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            axis.text=element_text(size=10)
      )  + 
      # ggtitle(xaxis) +
      scale_x_continuous(trans='log10') + 
      scale_color_manual(values = col2) %>% 
      suppressMessages()%>%suppressWarnings()
    return(g)
  }

  
  .plotViolinPost2=function(count.df,yaxis){
    axis.lab = unique(count.df$Sample)
    ylabs=gsub(" \\(", "\n\\(",yaxis)
    ylabs=gsub(paste0(" Top",n.topgnes), paste0("\nTop",n.topgnes),ylabs)
    
    
    g=ggplot(count.df, aes(x=Sample, y=(.data[[yaxis]]))) +
      ggtitle(yaxis) +
      theme_classic()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=rel(1)),
            legend.title=element_blank(), 
            axis.text=element_text(size=10),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            plot.title = element_text(size = 15, face = "bold")) +
      labs( y = ylabs)+
      geom_violin(aes(fill=as.factor(Sample))) +  
      scale_fill_manual(values = col2) +
      geom_boxplot(width=0) +
      scale_x_discrete(limits = as.vector(axis.lab))%>%suppressWarnings() 
    return(g)
    
  }
  
  
  ### Process SO object #####
  
  .seuratObject <- function(so.list,i) {
    ## i should be list of SO objects that have gone through pre-processing
    ## and samples should be Log normalized 
    # so.nf=so.nf.list[[i]]
    so.nf <- so.list[[i]]
    
    ## Optional: TSNE/UMAP for figures ####
    
    if(plot.outliers!="none"){
      so.nf.qcFiltr <- SCTransform(so.nf,do.correct.umi = TRUE,
                                vars.to.regress=vars_to_regress, 
                                return.only.var.genes = FALSE)
      so.nf.qcFiltr = FindVariableFeatures(object = so.nf.qcFiltr, 
                             nfeatures = nfeatures, 
                             mean.cutoff = c(low.cut, high.cut), 
                             dispersion.cutoff=c(low.cut.disp,high.cut.disp), 
                             selection.method=selection.method, verbose = FALSE)
      so.nf.qcFiltr <- RunPCA(object = so.nf.qcFiltr, 
                           npcs = npcs, verbose = FALSE,
                           seed.use = seed.for.PCA)
      if (plot.outliers=="umap") {
        
        so.nf.qcFiltr <- RunUMAP(object = so.nf.qcFiltr, 
                              reduction = "pca", 
                              dims = 1:npcs, 
                              seed.use=seed.for.UMAP)
        qcFiltr.df.plot = as.data.frame(
          so.nf.qcFiltr@reductions$umap@cell.embeddings)
        so.nf=AddMetaData(so.nf,qcFiltr.df.plot) 
        
      }else{
        
        so.nf.qcFiltr <- RunTSNE(object = so.nf.qcFiltr, 
                              reduction = "pca", 
                              dim.embed = 2, 
                              dims = 1:npcs, 
                              seed.use = seed.for.TSNE)
        qcFiltr.df.plot = as.data.frame(
          so.nf.qcFiltr@reductions$tsne@cell.embeddings)
        so.nf=AddMetaData(so.nf,qcFiltr.df.plot) 
      }
    }
    
    ## Filtering 
    ## set SO for filtering
    so <- so.nf
    
    
    ### Gene Filters ####    
    ### Remove VDJ genes 
    VDJgenesOut=c()
    if (filter.vdj.genes==TRUE) {
      allGenes = rownames(so)
      VDJgenes = c("TRBV","TRAV","TRBD","TRAJ","TRBJ")
      print("Removing VDJ genes. Genes removed...")
      print(length(allGenes))
      for (j in VDJgenes) {
        j=toupper(j)
        print(allGenes[grepl(j, toupper(allGenes))])
        VDJgenesOut=c(VDJgenesOut,allGenes[grepl(j, toupper(allGenes))])
        
        allGenes = allGenes[!grepl(j, toupper(allGenes))]  
      }
      so <- subset(so,features = allGenes)
    }    
    
    ### min # of Cells per gene
    ##CreateSeuratObject can apply min.cells but we do it in this module not in 
    # processing module.
    # CreateSeuratObject also calcs nfeature and nCounts. Testing showed that 
    # nfeature and nCounts calculations are not effected by min.cells when using 
    # CreateSeuratObject
    
    gene.cell.count=apply(so@assays$RNA@counts,1,function(x){sum(x>0)})
    so=subset(so, features=names(gene.cell.count)[(gene.cell.count>=min.cells)])
    
    
    ### Cell filters  ####
    ## Caluclate filter Metrics
    
    ## calculate Counts in Top 20 Genes
    so=.topNGenes(so,n.topgnes)
    
    ## Counts(umi) Filter 
    mad.ncounts.limits=.madCalc(so,'nCount_RNA',mad.ncounts.limits)
    mad.ncounts.limits=.checkLimits(mad.ncounts.limits)
    ncounts.limits=.checkLimits(ncounts.limits)
    mad.ncounts.limits=mad.ncounts.limits%>%round(digits=2)
    ncounts.limits=ncounts.limits%>%round(digits=2)
    
    ncounts.filter=((so@meta.data$nCount_RNA >= 
                      max(ncounts.limits[1],mad.ncounts.limits[1])) &
                     (so@meta.data$nCount_RNA <= 
                        min(ncounts.limits[2], mad.ncounts.limits[2]))
    ) 
    
    ## Gene Filter (nFeatrue)
    mad.nfeature.limits=.madCalc(so,'nFeature_RNA',mad.nfeature.limits)
    mad.nfeature.limits=.checkLimits(mad.nfeature.limits)
    nfeature.limits=.checkLimits(nfeature.limits)
    mad.nfeature.limits=mad.nfeature.limits%>%round(digits=2)
    nfeature.limits=nfeature.limits%>%round(digits=2)
    
    nfeature.filter= ((so@meta.data$nFeature_RNA >= 
                     max(nfeature.limits[1],mad.nfeature.limits[1])) &
                    (so@meta.data$nFeature_RNA <= 
                       min(nfeature.limits[2], mad.nfeature.limits[2]))
    )  
    
    ## Mitoc Fiter
    mad.mitoch.limits =.madCalc(so,'percent.mt',mad.mitoch.limits)  
    mad.mitoch.limits=.checkLimits(mad.mitoch.limits)
    mitoch.limits=.checkLimits(mitoch.limits)
    mad.mitoch.limits=mad.mitoch.limits%>%round(digits=2)
    mitoch.limits=mitoch.limits%>%round(digits=2)
    
    mitochPer.filter= ((so@meta.data$percent.mt >= 
                          max(mitoch.limits[1],mad.mitoch.limits[1])) &
                         (so@meta.data$percent.mt <= 
                            min(mitoch.limits[2], mad.mitoch.limits[2]))
    )    
    
    
    ## Complexity Filter
    mad.complexity.limits=.madCalc(so,'log10GenesPerUMI',mad.complexity.limits)
    mad.complexity.limits=.checkLimits(mad.complexity.limits)
    complexity.limits=.checkLimits(complexity.limits)
    mad.complexity.limits=mad.complexity.limits%>%round(digits=2)
    complexity.limits=complexity.limits%>%round(digits=2)
    
    complexity.filter= ((so@meta.data$log10GenesPerUMI >= 
                           max(complexity.limits[1],mad.complexity.limits[1])) &
                          (so@meta.data$log10GenesPerUMI <= 
                           min(complexity.limits[2], mad.complexity.limits[2]))
    )
    
    
    ## Top N Filter
    #  mad.topNgenes.limits = c(NA,5)
    # topNgenes.limits=c(NA,20)
    mad.topNgenes.limits=
      .madCalc(so,'pct_counts_in_top_N_genes',mad.topNgenes.limits)
    mad.topNgenes.limits=
      .checkLimits(mad.topNgenes.limits)
    topNgenes.limits=
      .checkLimits(topNgenes.limits)
    mad.topNgenes.limits=mad.topNgenes.limits%>%round(digits=2)
    topNgenes.limits=topNgenes.limits%>%round(digits=2)
    
    topN.filter= ((so@meta.data$pct_counts_in_top_N_genes >= 
                      max(topNgenes.limits[1],mad.topNgenes.limits[1])) &
                     (so@meta.data$pct_counts_in_top_N_genes <= 
                        min(topNgenes.limits[2], mad.topNgenes.limits[2]))
    )    
    
    ## doublets Filter
    if(do.doublets.fitler==T){
      doublets.fitler <- so@meta.data$Doublet%in%"singlet"
    }else{
      doublets.fitler=rep(TRUE,nrow(so@meta.data))
    }    
    


        ### Combine filters ####
    filterIndex <-ncounts.filter & 
      nfeature.filter & 
      mitochPer.filter & 
      complexity.filter & 
      topN.filter & 
      doublets.fitler
    filterIndex=as.data.frame(filterIndex)
    rownames(filterIndex) = rownames(so@meta.data)
    
    filter_matrix <- cbind(ncounts.filter,
                           nfeature.filter,
                           mitochPer.filter,
                           complexity.filter,
                           topN.filter,
                           doublets.fitler)
    rownames(filter_matrix)=rownames(so@meta.data)
    
    so.nf=AddMetaData(so.nf,filter_matrix,colnames(filter_matrix))
    
    
    ## print Filter resutls
    cat("\n\n")
    cat(i,":\n")
    cat(paste0("# of Cells before Filter: ",nrow(filter_matrix),"\n"))
    cat(paste0("# of Cells after all Filters: ",sum(filterIndex),"\n"))
    perc.remain <- (sum(filterIndex)/nrow(filter_matrix))*100
    perc.remain <- formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Percent Remaining: " ,perc.remain,"% \n\n")) 
    cat(paste0("# of Cells removed by individual Filters: \n"))
    print(colSums(filter_matrix==F))
    cat("\n\n")
    
    ## create Filter resutls table
    filtSum=matrix(nrow=1,ncol=3)
    colnames(filtSum)=c("Cells before Filtering",
                        "Cells after all Filters",
                        "Percent Remaining")
      filtSum[,"Cells before Filtering"]=nrow(filter_matrix)
      filtSum[,"Cells after all Filters"]=sum(filterIndex)
      filtSum[,"Percent Remaining"]=perc.remain
    topN.filterRename=paste0('% Counts in Top',n.topgnes,' Genes filter')
    filtTbl=colSums(filter_matrix==F)%>%t()%>%as.data.frame()
    filtTbl=rename(filtTbl,
                   'UMI Count (nCount_RNA)' = 'ncounts.filter',
                   'Gene Count (nFeature_RNA)' ='nfeature.filter',
                   '% Mitochondrial Genes (percent.mt)'='mitochPer.filter',
                   'Complexity (log10GenesPerUMI)'='complexity.filter',
                   'DoubletFinder (scDblFinder)'='doublets.fitler',
                   !!topN.filterRename :='topN.filter')
        colnames(filtTbl)=paste('Cells removed by ' ,colnames(filtTbl))
        
        filtSum=cbind(filtSum,filtTbl)
        filtSum[,"VDJ Genes Removed"]=VDJgenesOut%>%length
        rownames(filtSum)=i
        # tableGrob(filtSum)%>%plot
        
    ########################################################## #    
    ## create Filter Limits table
        topN.filterRename=paste0('% Counts in Top',n.topgnes,' Genes')
  cat('VDJ Genes Removed: ',length(VDJgenesOut), '\n')      
  cat('Minimum Cells per Gene: ',min.cells,'\n')
  cat('UMI Count (nCount_RNA) Limits: ',ncounts.limits,'\n')
  cat('MAD UMI Count (nCount_RNA) Limits: ',mad.ncounts.limits,'\n')
  cat('Gene Count (nFeature_RNA) Limits: ',nfeature.limits,'\n')
  cat('MAD Gene Count (nFeature_RNA) Limits: ',mad.nfeature.limits,'\n')
  cat('% Mitochondrial Genes (percent.mt) Limits: ',mitoch.limits,'\n')
  cat('MAD % Mitochondrial Genes (percent.mt) Limits: ',mad.mitoch.limits,'\n')
  cat('Complexity (log10GenesPerUMI) Limits: ',complexity.limits,'\n')
  cat('MAD Complexity (log10GenesPerUMI) Limits: ',mad.complexity.limits,'\n')
  cat(topN.filterRename,' Limits: ',topNgenes.limits,'\n')
  cat('MAD ',topN.filterRename,' Limits: ',mad.topNgenes.limits,'\n')
  cat('Doublets Filter: ',do.doublets.fitler,'\n')
  
    
    
    ## create Filter Limits table
    FiltLmts=matrix(nrow=1,ncol=12)
    colnames(FiltLmts)=c("Minimum Cells per Gene",
                        "UMI Count (nCount_RNA)",
                        "MAD UMI Count (nCount_RNA)",
                        "Gene Count (nFeature_RNA)",
                        "MAD Gene Count (nFeature_RNA)",
                        "% Mitochondrial Genes (percent.mt)",
                        "MAD % Mitochondrial Genes (percent.mt)",
                        "Complexity (log10GenesPerUMI)",
                        "MAD Complexity (log10GenesPerUMI)",
                        paste0(topN.filterRename,"") ,
                        paste0('MAD ',topN.filterRename,'') ,
                        "DoubletFinder (scDblFinder)"
                        )

    FiltLmts[,'Minimum Cells per Gene']=min.cells
    FiltLmts[,'UMI Count (nCount_RNA)']=
      paste0(c("Low:","High:"),ncounts.limits)%>%paste(collapse = "\n")
    FiltLmts[,'MAD UMI Count (nCount_RNA)']=
      paste0(c("Low:","High:"),mad.ncounts.limits)%>%paste(collapse = "\n")
    FiltLmts[,'Gene Count (nFeature_RNA)']=
      paste0(c("Low:","High:"),nfeature.limits)%>%paste(collapse = "\n")
    FiltLmts[,'MAD Gene Count (nFeature_RNA)']=
      paste0(c("Low:","High:"),mad.nfeature.limits)%>%paste(collapse = "\n")
    FiltLmts[,'% Mitochondrial Genes (percent.mt)']=
      paste0(c("Low:","High:"),mitoch.limits)%>%paste(collapse = "\n")
    FiltLmts[,'MAD % Mitochondrial Genes (percent.mt)']=
      paste0(c("Low:","High:"),mad.mitoch.limits)%>%paste(collapse = "\n")
    FiltLmts[,'Complexity (log10GenesPerUMI)']=
      paste0(c("Low:","High:"),complexity.limits)%>%paste(collapse = "\n")
    FiltLmts[,'MAD Complexity (log10GenesPerUMI)']=
      paste0(c("Low:","High:"),mad.complexity.limits)%>%paste(collapse = "\n")
    FiltLmts[,paste0(topN.filterRename,'')]=
      paste0(c("Low:","High:"),topNgenes.limits)%>%paste(collapse = "\n")
    FiltLmts[,paste0('MAD ',topN.filterRename,'')]=
      paste0(c("Low:","High:"),mad.topNgenes.limits)%>%paste(collapse = "\n")
    FiltLmts[,'DoubletFinder (scDblFinder)']=do.doublets.fitler
    rownames(FiltLmts)=i
    
    ### Apply Filters ####
    
    so <- subset(so, cells = 
                   rownames(filterIndex[which(filterIndex[,1]==T), ,drop=F])
    )
    
    ## Select filters that remove cells
    filter_matrix=
      filter_matrix[,is.na(apply(
        filter_matrix,2,function(i){match(FALSE,i)}))==F
      ]
    
    
    ## Optional: create TSNE plots ####
    if(plot.outliers!="none"){   filter.plot.cols=3
    gtsne.list=lapply(colnames(filter_matrix),
                      function(x){.runTsnepPlot(x,filter_matrix,
                                                so.nf.qcFiltr,plot.outliers)})
    gtsne.all <- arrangeGrob(grobs = gtsne.list, ncol = filter.plot.cols)
    
    so2.list <- list(Filtered=so,
                     FilteringMeta=so.nf@meta.data,
                     FiltLmts=FiltLmts,
                     filtSum=filtSum,
                     TSNEfilter=gtsne.all
    )
    }else{
      so2.list <- list(Filtered=so,
                       FilteringMeta=so.nf@meta.data,
                       FiltLmts=FiltLmts,
                       filtSum=filtSum) 
      
    }
    
    return(so2.list)
    
           
  }
  
  
  
  ## --------------- ##
  ## Main Code Block ####
  ## --------------- ##
  ## make sure that plot.outliers is character not boolean and case insensitive
  plot.outliers =tolower(as.character(plot.outliers))

  
  ## Caluclate filter Metrics
  so.nf.list=lapply(names(object), function(i){
    so=object[[i]]
    
    ## calculate Counts in Top 20 Genes
    ##calculated after min.cell filter as well
    so=.topNGenes(so,n.topgnes) 
    
    ## Annotate Doublets: ####
    ## Gene filter does not effect doublet ident and so not recalculated
    
    if( do.doublets.fitler==T){
      sce <- as.SingleCellExperiment(so)
      set.seed(123)
      sce.dbl <- scDblFinder(sce)%>%suppressWarnings()
      sce.class <- sce.dbl$scDblFinder.class
      so <- AddMetaData(so,sce.class,"Doublet")
    }else{ print('doublets Identification not Run')}
    return(so)
  })
  names(so.nf.list)=names(object)
  
  
  ### Run Filtering Function ####

  so.f.out <- lapply(names(so.nf.list), 
                     function(i){.seuratObject(so.nf.list,i)})
  names(so.f.out)=names(so.nf.list)

  
  #### Get Filtered SO
  so.f.list <- lapply(so.f.out,function(x) {x[['Filtered']]})
  
  #### Get Cell filtering Metadata
  so.nf.list.meta=lapply(names(so.f.out),
                         function(x){so.f.out[[x]][['FilteringMeta']]})
  names(so.nf.list.meta)=names(so.f.out)
  
  #### Get filtering Summary Tables
  
  filtSum=lapply(so.f.out,`[[`,'filtSum')%>%do.call(rbind,.)%>%as.data.frame()
  FiltLmts=lapply(so.f.out,`[[`,'FiltLmts')%>%do.call(rbind,.)%>%as.data.frame()
  
  
  
  
  ### Collect QC figures ####

  ##Polychorme will create a palette of random colores based on a few seed 
  ## colors. Here i used RcolorBrewer set1 color palette as the seed. Can not 
  ## Reproduce exact pallet each time so creating 1 pallet that will be 
  ## hardcoded in the Filter_QC and Combine_Normalize functions
  
  # library(Polychrome)
  # P100 = createPalette(100,  brewer.pal(5, "Set1"), target = 'normal')
  # swatch(P100)
  
  col2=c("#E6161C","#2E7AB6","#53B251","#9A53A6","#FC7800","#FFDA16","#FC00B2",
         "#7622FE","#7D532E","#00FBF4","#F126FE","#DCDEE6","#22FD0D","#FD839B",
         "#005D51","#E2E19E","#C6F022","#FBAD70","#008AFE","#F8AED8","#A5001C",
         "#B42685","#5835B4","#CEC0FF","#654B63","#FF7DF1","#F91C7F","#5C6D16",
         "#16FABB","#7DCCAF","#8CD7F4","#C86EFC","#FFCDBD","#882632","#B0921C",
         "#B65500","#FD0DDB","#65FE84","#9F9CFF","#948C7E","#8FAA00","#758EA1",
         "#C5166C","#93456A","#0DC4FB","#AF0DB9","#F5AAFD","#554984","#C57A7A",
         "#0D9D9E","#B7EEAD","#FD7168","#B700FB","#85EC00","#3B8655","#AE8EB5",
         "#E7E965","#FE87D2","#FC9B16","#AFF27A","#8876FF","#C3D2B4","#FBD092",
         "#1626DB","#FE5ACD","#E4AA0D","#A96A0D","#B3966A","#00E5FD","#004F6E",
         "#7D00A3","#FE947F","#BC87E9","#4F4726","#ED3D5A","#C79BA5","#5CA8FD",
         "#005FBF","#1CBB84","#B2EDED","#FD5C32","#FCD6FF","#F668AB","#97C1FB",
         "#7F1CD7","#94AF66","#BE16A7","#DE78A1","#22B80D","#7568BB","#82F4A3",
         "#709686","#FC6086","#BB60A3","#9F1658","#EE6CFF","#2E8BA8","#47E9C8",
         "#8387BA","#BCCDEC" )
  # swatch(col2)
  
  ##Alternative is to combine RcolorBrewer palettes  to create larger palette
  ## More limited on number of colors (61)
  
  # display.brewer.all()
  # col2 <- c(brewer.pal(8, "Set1"),
  #           brewer.pal(12, "Paired"),
  #           brewer.pal(7, "Dark2"),
  #           brewer.pal(6, "Accent"),
  #           brewer.pal(12, "Set3"),
  #           brewer.pal(8,"Set2"),            
  #           brewer.pal(11, "Spectral")
  #           )%>%unique
  
  
  
  features=c("orig.ident",
             "nCount_RNA","nFeature_RNA",
             "percent.mt","log10GenesPerUMI",
             "pct_counts_in_top_N_genes","DoubletFinder (scDblFinder)")
  v=features[features%in%c('nCount_RNA','nFeature_RNA',
                           'percent.mt','log10GenesPerUMI')]
  
  
  
  #### Combine meta.data tables ####
  filtTable.all=data.frame()
  rawTable.all=data.frame()
  for(s in names(so.f.list)) {
    filtTable=so.f.list[[s]]@meta.data
    filtTable=filtTable[,colnames(filtTable)%in%features,drop=F]
    filtTable$Sample=s
    filtTable$filt='filt'
    
    rawTable=so.nf.list[[s]]@meta.data
    rawTable=rawTable[,colnames(rawTable)%in%features,drop=F]
    rawTable$Sample=s
    rawTable$filt='raw'
    
    filtTable.all=rbind(filtTable.all,filtTable)
    rawTable.all=rbind(rawTable.all,rawTable)
  }
  
  table.meta=rbind(filtTable.all,rawTable.all)
  table.meta$nFeature_RNA=as.numeric(table.meta$nFeature_RNA)
  table.meta$filt=factor(table.meta$filt,levels = c('raw','filt'))
  
  topN.filterRename=paste0('% Counts in Top',n.topgnes,' Genes')
  
  table.meta=rename(table.meta,
                    'UMI Count (nCount_RNA)' = 'nCount_RNA',
                    'Gene Count (nFeature_RNA)' ='nFeature_RNA',
                    '% Mitochondrial Genes (percent.mt)'='percent.mt',
                    'Complexity (log10GenesPerUMI)'='log10GenesPerUMI',
                    !!topN.filterRename :='pct_counts_in_top_N_genes'
  )
  


  #add Sample Cutoffs
  
  table.meta= merge(table.meta,FiltLmts,
        by.x='Sample',by.y='row.names',
        suffix=c("","_Filters"),all.x=T)
        
  
  v=c('UMI Count (nCount_RNA)' ,
      'Gene Count (nFeature_RNA)',
      '% Mitochondrial Genes (percent.mt)',
      'Complexity (log10GenesPerUMI)',
      topN.filterRename)
  
  #### Create Filter QC Plots ####
  
  ### Violin Plots
  
  ### Violin Plots for each meta.data metric  
  violin.list=lapply(v,function(x){.plotViolin2(table.meta,x)})
  names(violin.list)=v
  colnames(table.meta)
  
  ### Combine Violin Plots
  violin.list.mod=.combinePlots(violin.list)
  
  violin.grob=ggarrange(plotlist=violin.list.mod,
                        nrow=length(violin.list),
                        common.legend = T,
                        legend = 'right')  
  
  
  ### Scatter Plots

  ### scatter Plots for each meta.data metric  
  scatter.list=lapply(v[!v%in%'UMI Count (nCount_RNA)'],
                      function(x){.plotScatter2(table.meta,x)})
  names(scatter.list)=v[!v%in%'UMI Count (nCount_RNA)']
  
  ### Combine scatter Plots
  scatter.list.mod=.combinePlots(scatter.list)
  
  scatter.grob=ggarrange(plotlist=scatter.list.mod,
                         nrow=length(scatter.list),
                         common.legend = T,
                         legend = 'right')  
  
  
  #### Create Post Filter QC figure ####
  
  qc.df.post=table.meta[table.meta$filt=='filt',]
  
  
  ### Post Filter Summary - Scatter
  
  scatter.allsamples=lapply(v,
                      function(y){.plotScatterPost2(
                        qc.df.post,'UMI Count (nCount_RNA)',y)})
  names(scatter.allsamples)=v
  
  scatter.allsamples.grob=ggarrange(plotlist=scatter.allsamples,
                                    ncol=1,
                                    common.legend = T,
                                    legend = 'right')
  
  
  ### Post Filter Summary - histogram
  
  hist.allsamples=lapply(v,function(x){.plotHistPost2(qc.df.post,x)})
  names(hist.allsamples)=v
  
  hist.allsamples.grob=ggarrange(plotlist=hist.allsamples,
                                 ncol=1,
                                 common.legend = T,
                                 legend = 'right') %>% 
    suppressMessages()%>%suppressWarnings()
  
  
  ### Post Filter Summary - Violin
  
  violin.allsamples=lapply(v,function(x){.plotViolinPost2(qc.df.post,x)})
  names(violin.allsamples)=v
  violin.allsamples.grob=ggarrange(plotlist=violin.allsamples,
                                   ncol=1,
                                   common.legend = T,
                                   legend = 'right')
  
  violin.allsamples.grob=annotate_figure(violin.allsamples.grob, 
                                       top = text_grob("Filtered QC Summary", 
                                                   face = "bold", size = 14))
  
  
  ### Post Filter Summary - combined Scatter + Histogram
  
  postFilter.grobs=ggarrange(
    ggarrange(plotlist=scatter.allsamples,
              ncol=1,legend = 'none'),
    ggarrange(plotlist=hist.allsamples,
              ncol=1,legend = 'none'),
    # ggarrange(plotlist=violin.allsamples,ncol=1,legend = 'none'),
    legend.grob=get_legend(scatter.allsamples[[1]]),
    ncol=2,
    legend='right') %>% 
    suppressMessages()%>%suppressWarnings()
  
  postFilter.grobs=annotate_figure(postFilter.grobs, 
                                   top = text_grob("Filtered QC Summary", 
                                                   face = "bold", size = 14))
  
  
  ### TSNE Plots
  if (plot.outliers!="none") {
    qcFiltr.grobs <- lapply(so.f.out,function(x) x[['TSNEfilter']])
  }
  
  
  
  
  # cellcount.nf <- lapply(so.nf.list, function(x) dim(x)[2])
  # cellcount.f <- lapply(so.f.list, function(x) dim(x)[2])
  # sum.before <- sum(unlist(cellcount.nf))
  # sum.after <- sum(unlist(cellcount.f))
  # cat("\n\nTotal number of cells before filtering:", sum.before, "\n")
  # cat("Total number of cells after filtering: ", sum.after,"\n\n")
  # cat("Percentage cells remaining after filtering:", 
  #     (sum.after/sum.before)*100,"\n")
  
  gc(full = TRUE) 
  

  
  ### Output ####
  out=list(object=so.f.list,
       FilteringTables=list(
         FilteringMeta=so.nf.list.meta,
         FilteringCounts=filtSum,
         FilteringLimits=FiltLmts
       ),
       plots=list(
         ViolinPlotCombine=violin.grob,
         ViolinPlot=violin.list,
         ScatterPlotCombine=scatter.grob,
         Scatter=scatter.list,
         
         PostFilterCombined=postFilter.grobs,
         ViolinPostFilter=violin.allsamples.grob,
         ScatterPostFilter=scatter.allsamples.grob,
         HistogramPostFilter=hist.allsamples.grob
       )
  )
  
  if(plot.outliers!='none'){
    out$plots=c(out$plots,TSNEFeature=qcFiltr.grobs)
    
  } 
    return(out)

  
}


