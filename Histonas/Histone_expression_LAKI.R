library(ggplot2)
library(ggsignif)

#cargamos el resultado de DESeq y los id y los nombres de las histonas, los genes y los pseudogenes
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_all.obj")
load("/home/antotartier/data/histones/gene_sets/mouse/gene_set_histonas_raton.obj")
load("/home/antotartier/data/histones/gene_sets/mouse/histone_genes.obj")
load("/home/antotartier/data/histones/gene_sets/mouse/histone_pseudogenes.obj")


#lo convertimos en un dataframe
res_all_df<-lapply(res_all,as.data.frame)
for(x in seq_along(res_all_df)){
  rownames(res_all_df[[x]])<-gsub("\\.[0-9]+$","",rownames(res_all_df[[x]]))#eliminamos el punto de los ensemble IDs
}

#nos quedamos con las histonas
res_hist<-list()
for(name in names(res_all_df)){
  x<-res_all_df[[name]]
  x<-x[rownames(x)%in%gene_set$Ensembl_Accession_ID,]
  res_hist[[name]]<-x
} 

#sustituimos los id por los nombres
id_to_name<-setNames(gene_set$Marker_Symbol,gene_set$Ensembl_Accession_ID)
for(x in seq_along(res_hist)){
  rownames(res_hist[[x]])<-id_to_name[rownames(res_hist[[x]])]
}


#hacemos una lista de los nombres de las histonas cluster y no cluster
cluster<-gene_set$Marker_Symbol[grepl("cluster",gene_set$Marker_Symbol) & !grepl("lncRNA",gene_set$BioTypes)]
cluster<-cluster[cluster %in% histone_genes$Marker_Symbol]
no_cluster<-gene_set$Marker_Symbol[!grepl("cluster",gene_set$Marker_Symbol) & !grepl("lncRNA",gene_set$BioTypes)]
no_cluster<-no_cluster[!(no_cluster %in% histone_pseudogenes$Marker_Symbol)]

#hacemos una lista de los nombres de las histonas por familia
H1<-histone_genes$Marker_Symbol[grepl("H1",histone_genes$Marker_Symbol)]
H2A<-histone_genes$Marker_Symbol[grepl("(?i)H2a",histone_genes$Marker_Symbol)]
H2B<-histone_genes$Marker_Symbol[grepl("H2b",histone_genes$Marker_Symbol)]
H3<-histone_genes$Marker_Symbol[grepl("H3",histone_genes$Marker_Symbol)]
H4<-histone_genes$Marker_Symbol[grepl("H4",histone_genes$Marker_Symbol)]

#extraemos los genes up-regulados
up_reg<-list()
for(name in names(res_hist)){
  x<-res_hist[[name]]
  up<-x[x$log2FoldChange>0 & x$padj<0.05,]
  up_reg[[name]]<-up
}
#definimos un patrón para seleccionar los genes dependientes de la replicación
dpts<-c("cluster","H2Ax","H4c16")


#clasificamos los genes
for(x in seq_along(res_hist)){
y<-res_hist[[x]]
res_hist[[x]]<-y %>%
  mutate(type=ifelse(rownames(y)%in%histone_pseudogenes$Marker_Symbol, "pseudogen",
                     ifelse(rownames(y)%in%histone_genes$Marker_Symbol,"gen", "Lnc")),
         clustering=ifelse(rownames(y)%in%cluster, "cluster",
                           ifelse(rownames(y)%in%no_cluster,"No cluster", "unknown")),
         familia=ifelse(rownames(y)%in%H1, "H1",
                        ifelse(rownames(y)%in%H2A,"H2a",
                               ifelse(rownames(y)%in%H2B,"H2b",
                                      ifelse(rownames(y)%in%H3,"H3",
                                             ifelse(rownames(y)%in%H4,"H4", "Unknown"))))),
         reg=ifelse(grepl(paste(dpts, collapse = "|"),rownames(y)) & !grepl("H3f4",rownames(y)),"HDR",
                    ifelse(!grepl(paste(dpts, collapse = "|"),rownames(y)) | grepl("H3f4",rownames(y)),"HIR"
                           ,"unknown")))
}

#nos quedamos con los resultados estadísticamente significativos, y eliminamos los posibles NA (normalmente outliers)

res_hist_sig<-list()
for(name in names(res_hist)){
  x<-res_hist[[name]]
  x<-x[x$padj<0.05,]
  res_hist_sig[[name]]<-x
}
res_hist_sig<-lapply(res_hist_sig,na.omit)

#dejamos también resultados sin el filtro de la significación pero sin NA
res_hist<-lapply(res_hist,na.omit)

#eliminamos los pseudogenes
res_hist_sig_nopseudo<-list()
for(name in names(res_hist_sig)){
  x<-res_hist_sig[[name]]
  x<-x[x$type=="gen",]
  res_hist_sig_nopseudo[[name]]<-x
}

res_hist_nopseudo<-list()
for(name in names(res_hist)){
  x<-res_hist[[name]]
  x<-x[x$type=="gen",]
  res_hist_nopseudo[[name]]<-x
}

#creamos una función que haga los boxplots
myBoxplot<-function(data,group,sizesig,pointsize,filter=F,comp_vec){
  group <- enquo(group)
  group_col <- ensym(group)  # Capture the group column symbolically
  
  generate_comparisons <- function(levels) {
    combn(levels, 2, simplify = FALSE)
  }#genera todas las posibles comparaciones por pares a partir de los niveles
  
  levels_list <- unique(data[[as.character(group_col)]])  # Extract unique levels
  comparisons <- generate_comparisons(levels_list)
  
  plot<-ggplot(data, aes(x=!!group, y=log2FoldChange, fill=!!group)) +
    geom_boxplot(alpha = 0.80, coef=0,outlier.shape = NA) +
    stat_summary(fun.data = mean_sd, geom = "errorbar", width=0.25) +
    geom_point(aes(fill = !!group), size = pointsize, shape = 21, position = position_jitterdodge()) +
    theme(axis.title = element_text(size = 14),text = element_text(size = 18),
          axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, hjust = 0.5))
  if(!filter){
    plot<-plot + geom_signif(comparisons = comparisons, test = "t.test",y_position = max(data$log2FoldChange)+0.1,step_increase = 0.1,
                             textsize = sizesig)}
  else{
    plot<-plot +geom_signif(comparisons = list(comp_vec), test = "t.test",y_position = max(data$log2FoldChange)+0.1,step_increase = 0.1,
                            textsize = sizesig)}#este condicional permite quedarse solo con ciertas comparaciones
  
  return(plot)
}


#hacemos los boxplots, menos el de pseudogenes porque no hay significativos
for(x in seq_along(res_hist_sig)){
  y<-res_hist_sig[[x]]
  z<-sum(y$type=="pseudogen")
  print(z)
}

#juntamos los datos
res_hist_sig_nopseudo_df<-lapply(res_hist_sig_nopseudo,as.data.frame)
res_hist_sig_nopseudo_df<-bind_rows(res_hist_sig_nopseudo_df,.id = "tejido")

res_hist_nopseudo_df<-lapply(res_hist_nopseudo,as.data.frame)
res_hist_nopseudo_df<-bind_rows(res_hist_nopseudo_df,.id = "tejido")

#hacemos los gráficos
clustering<-myBoxplot(res_hist_sig_nopseudo_df,group = clustering,sizesig = 4)+facet_wrap(~tejido,scales = "free")
family_sig<-myBoxplot(res_hist_sig_nopseudo_df,group = familia,sizesig = 4,filter=T,comp_vec=NULL)+facet_wrap(~tejido,scales = "free_y")
reg_sig<-myBoxplot(res_hist_sig_nopseudo_df,group = reg,sizesig = 4)+facet_wrap(~tejido,scales = "free_y")

family<-myBoxplot(res_hist_nopseudo_df,group = familia,sizesig = 4,pointsize = 1,filter=T,comp_vec=NULL)+facet_wrap(~tejido,scales = "free_y")
reg<-myBoxplot(res_hist_nopseudo_df,group = reg,sizesig = 4,pointsize = 1)+facet_wrap(~tejido,scales = "free_y")

#guardamos los plots que nos interesan para el TFM
ggsave("/home/antotartier/data/histones/LAKI/boxplots_family.pdf",family)
ggsave("/home/antotartier/data/histones/LAKI/boxplots_reg.pdf",reg)



#el principal regulador de la expresion de las histonas dependientes de la replicación es NPAT, un activador de la transcripción de las mismas que se activa durante la fase S del ciclo celular gracias a la fosforilación del mismo por parte de CDK2/CycE. Dados los cambios en la expresión de las histonas observados procedo a analizar la expresión de estos genes en ambas líneas celulares en busca de una posible causa para la down-regulación observada
# load("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/DESeq_conj.obj")
# rownames(DESeq)<-gsub("\\.[0-9]+$","",rownames(DESeq))
# NPAT<-plotCounts(DESeq, gene="ENSG00000149308", intgroup=c("Passage","Cells"),returnData = T)
# NPAT_plot<-ggplot(NPAT, aes(x=Passage,y=count,colour=Cells))+geom_point()+ggtitle("Count plot NPAT")
# CDK2<-plotCounts(DESeq, gene="ENSG00000123374", intgroup=c("Passage","Cells"),returnData = T)
# CDK2_plot<-ggplot(CDK2, aes(x=Passage,y=count,colour=Cells))+geom_point()+ggtitle("Count plot CDK2")
# CycE<-plotCounts(DESeq, gene="ENSG00000105173", intgroup=c("Passage","Cells"),returnData = T)
# CycE_plot<-ggplot(CycE, aes(x=Passage,y=count,colour=Cells))+geom_point()+ggtitle("Count plot CycE")
# 
# #guardamos todos los plots
# png(filename = "~/histones/senescencia/NPAT_CDK2_CycE_cts.png", height =5000 , width =18000 ,res =1100 )
# grid.arrange(NPAT_plot,CDK2_plot,CycE_plot,ncol=3)
# dev.off()

#los resultados muestran cambios consistentes en ambas líneas celulares en la CDK2, concretamente la downregulación de la misma con un LFC de -1,4. Dado que si está disminye también disminuye la fosforilación de NPAT, este cambio podría ser uno de los responsables de la disminución de la expresión de las histonas. Además, teniendo en cuenta que la disminución es significativamente mayor en las histonas en cluster, que en su mayoría son dependientes de la replicación (aquellas reguladas por NPAT, CDK2 y CycE), esta causa resulta verosimil. No obstante, pese a que la expresión de las CDKs teóricamente permanece constante a lo largo del ciclo celular, para confirmar que este cambio está asociado a la senescencia es necesario descartar que este sea consecuencia del arresto de la proliferación que se produce durante la senescencia. Para ello realizaremos este análisis en datos de arresto por falta de suero en el que se produce un arresto de la proliferación independiente a la senscencia.

#----------------------------------------------------------------------------------------------------------------------------------------
#calculamos el LFC medio de las histonas para cada tejido
res_hist<-na.omit(res_hist)
LFC_hist<-list()
for(name in names(res_hist)){
x<-res_hist[[name]]
y<-mean(x$log2FoldChange)
LFC_hist[[name]]<-data.frame(tejido=name,LFC=y)
}
LFC_hist<-do.call(rbind,LFC_hist)

#machos
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_M.obj")
#lo convertimos en un dataframe
res_M_df<-lapply(res_M,as.data.frame)
res_M_df<-lapply(res_M_df,na.omit)#eliminamos los NA
for(x in seq_along(res_M_df)){
  rownames(res_M_df[[x]])<-gsub("\\.[0-9]+$","",rownames(res_M_df[[x]]))#eliminamos el punto de los ensemble IDs
}
#nos quedamos con las histonas
res_hist_M<-list()
for(name in names(res_M_df)){
  x<-res_M_df[[name]]
  x<-x[rownames(x)%in%gene_set$Ensembl_Accession_ID,]
  res_hist_M[[name]]<-x
} 
#calculamos el LFC medio de las histonas para cada tejido
LFC_hist_M<-list()
for(name in names(res_hist_M)){
  x<-res_hist_M[[name]]
  y<-mean(x$log2FoldChange)
  LFC_hist_M[[name]]<-data.frame(tejido=name,LFC=y)
}
LFC_hist_M<-do.call(rbind,LFC_hist_M)

#hembras
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_F.obj")
#lo convertimos en un dataframe
res_H_df<-lapply(res_F,as.data.frame)
res_H_df<-lapply(res_H_df,na.omit)#eliminamos los NA
for(x in seq_along(res_H_df)){
  rownames(res_H_df[[x]])<-gsub("\\.[0-9]+$","",rownames(res_H_df[[x]]))#eliminamos el punto de los ensemble IDs
}
#nos quedamos con las histonas
res_hist_H<-list()
for(name in names(res_H_df)){
  x<-res_H_df[[name]]
  x<-x[rownames(x)%in%gene_set$Ensembl_Accession_ID,]
  res_hist_H[[name]]<-x
} 
#calculamos el LFC medio de las histonas para cada tejido
LFC_hist_H<-list()
for(name in names(res_hist_H)){
  x<-res_hist_H[[name]]
  y<-mean(x$log2FoldChange)
  LFC_hist_H[[name]]<-data.frame(tejido=name,LFC=y)
}
LFC_hist_H<-do.call(rbind,LFC_hist_H)

#guardamos los resultados
save(LFC_hist,file="/home/antotartier/data/histones/LAKI/LFC_hist.obj")
save(LFC_hist_M,file="/home/antotartier/data/histones/LAKI/LFC_hist_M.obj")
save(LFC_hist_H,file="/home/antotartier/data/histones/LAKI/LFC_hist_F.obj")

#----------------------------------------------------------------------------------------------------------------------------------------
#Repetimos lo mismo pero con los resultados estadísiticamente significativos

LFC_hist_sig<-list()
for(name in names(res_hist_sig)){
  x<-res_hist_sig[[name]]
  y<-mean(x$log2FoldChange)
  LFC_hist_sig[[name]]<-data.frame(tejido=name,LFC=y)
}
LFC_hist_sig<-do.call(rbind,LFC_hist_sig)

#machos
#nos quedamos con las histonas con cambio significativos
res_hist_M_sig<-list()
for(name in names(res_M_df)){
  x<-res_M_df[[name]]
  x<-x[rownames(x)%in%gene_set$Ensembl_Accession_ID,]
  x<-x[x$padj<0.05,]
  res_hist_M_sig[[name]]<-x
}
#calculamos el LFC medio de las histonas para cada tejido
LFC_hist_M_sig<-list()
for(name in names(res_hist_M_sig)){
  x<-res_hist_M_sig[[name]]
  y<-mean(x$log2FoldChange)
  LFC_hist_M_sig[[name]]<-data.frame(tejido=name,LFC=y)
}
LFC_hist_M_sig<-do.call(rbind,LFC_hist_M_sig)

#hembras
#nos quedamos con las histonas con cambio significativos
res_hist_H_sig<-list()
for(name in names(res_H_df)){
  x<-res_H_df[[name]]
  x<-x[rownames(x)%in%gene_set$Ensembl_Accession_ID,]
  x<-x[x$padj<0.05,]
  res_hist_H_sig[[name]]<-x
} 
#calculamos el LFC medio de las histonas para cada tejido
LFC_hist_H_sig<-list()
for(name in names(res_hist_H_sig)){
  x<-res_hist_H_sig[[name]]
  y<-mean(x$log2FoldChange)
  LFC_hist_H_sig[[name]]<-data.frame(tejido=name,LFC=y)
}
LFC_hist_H_sig<-do.call(rbind,LFC_hist_H_sig)

#guardamos los resultados
save(LFC_hist_sig,file="/home/antotartier/data/histones/LAKI/LFC_hist_sig.obj")
save(LFC_hist_M_sig,file="/home/antotartier/data/histones/LAKI/LFC_hist_M_sig.obj")
save(LFC_hist_H_sig,file="/home/antotartier/data/histones/LAKI/LFC_hist_F_sig.obj")