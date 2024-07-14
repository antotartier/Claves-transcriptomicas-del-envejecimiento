#cargamos las libreríass necesarias
library(openxlsx)
library(dplyr)
library(fgsea)
library(data.table)
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(ggsignif)

#funciones

#heatmap
myHeatmap<-function(data,title){
  heatmap(as.matrix(data),
          scale = "row",  
          margins = c(8, 8),      
          cexRow = 0.5,             # Row text size
          cexCol = 1
  )
}

#genera el rank, a partir del df del resultado DESeq
get_rank<-function(x){
  res_gsea<-as.data.frame(x)
  res_gsea<-res_gsea[res_gsea$baseMean>10,]#eliminamos genes con pocas lecturas
  rownames(res_gsea)<-gsub("\\.[0-9]+$","",rownames(res_gsea))#eliminamos el punto de los ensemble IDs
  res_gsea<-mutate(res_gsea,rank= -log10(pvalue+min(res_gsea$pvalue[res_gsea$pvalue!=0], na.rm= T))*sign(log2FoldChange))
  res_gsea<-na.omit(res_gsea)#eliminamos na y duplicados
  res_gsea<-distinct(res_gsea)
  res_gsea_vector<-res_gsea$rank
  names(res_gsea_vector)<-row.names(res_gsea)
  res_gsea_vector<-sort(res_gsea_vector, decreasing = T)
  return(res_gsea_vector)
}

#boxplots con prueba T
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

#-------------------------------------------------------------------------------

#HEATMAP DE LAS HISTONAS


#cargamos el sample_info (nos quedamos con las muestras de hígado), el gene_set de histonas y los resultados de DESEq
gene_set<-read.xlsx("/home/antotartier/data/histones/gene_sets/mouse/gene_set_histonas_raton.xlsx")
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/DDS_exp_all.obj")


#identificamos los genes cluster, los pseudogenes y los lncRNA
gene_set <- gene_set %>%
  mutate(Marker_Symbol = ifelse(grepl("pseudogene", BioTypes), paste0(Marker_Symbol, "_pseudogene"),Marker_Symbol))%>%
  mutate(Marker_Symbol = ifelse(grepl("cluster", Marker_Name),paste0(Marker_Symbol,"_cluster"),Marker_Symbol))%>%
  mutate(Marker_Symbol = ifelse(grepl("lncRNA", BioTypes),paste0(Marker_Symbol,"_lnc"),Marker_Symbol))
save(gene_set,file = "/home/antotartier/data/histones/gene_sets/mouse/gene_set_histonas_raton.obj")


#normalizamos los counts y los extraemos
vsd<-lapply(DDS_exp,vst,blind=F)
counts<-lapply(vsd,assay)


#eliminamos el punto de los ensemble IDs y nos quedamos con los de las histonas
for(x in 1:length(counts)){
  rownames(counts[[x]])<-gsub("\\.[0-9]+$","",rownames(counts[[x]]))
}
cts_hist<-list()
for(name in names(counts)){
  y<-counts[[name]]
  y<-y[rownames(y) %in% gene_set$Ensembl_Accession_ID,]
  cts_hist[[name]]<-y
}


#sustituimos los id por los nombres
id_to_name<-setNames(gene_set$Marker_Symbol,gene_set$Ensembl_Accession_ID)
for(x in 1:length(cts_hist)){
  rownames(cts_hist[[x]])<-id_to_name[rownames(cts_hist[[x]])]
}


#sustituimos el número del tubo por la condición
sample_info<-sample_info%>%
  mutate(Genotipo=substr(SampleID,nchar(SampleID)-2,nchar(SampleID)-1)) %>%#los ultimos 2 cracteres antes del ultimo son el genotipo
  mutate(Sexo=substr(SampleID,nchar(SampleID),nchar(SampleID)))%>% # el ultimo caracter es el sexo
  mutate(Condición=paste(Genotipo,Sexo,sep="_"))%>%
  mutate(Condición= stri_replace_all_fixed(Condición,pattern = c("KO","F"),replacement = c("LAKI","H"),vectorize_all = F))
ausent <- c("Tube_39","Tube_56","Tube_57","Tube_66","Tube_87") #error while transfering data
FACE <- c("Tube_115","Tube_116","Tube_117","Tube_118","Tube_119","Tube_120")
sample_info<-sample_info[!grepl(paste(ausent,collapse = "|"),sample_info$SampleName)]
sample_info<-sample_info[!grepl(paste(FACE,collapse = "|"),sample_info$SampleName)]
sample_info$SampleName<-gsub("Tube_","",sample_info$SampleName)

tube_to_condition<-setNames(sample_info$Condición,sample_info$SampleName)
for(x in 1:length(cts_hist)){
  colnames(cts_hist[[x]])<-tube_to_condition[colnames(cts_hist[[x]])]
}


#hacemos el heatmap y guardamos los resultados
png(filename = "/home/antotartier/data/histones/LAKI/heatmap_LAKI_Riñón.png", width = 9600, height = 10000, res = 1100)
myHeatmap(cts_hist$Riñón,"Riñón")
dev.off()

png(filename = "/home/antotartier/data/histones/LAKI/heatmap_LAKI_Ileon.png", width = 9600, height = 10000, res = 1100)
myHeatmap(cts_hist$Íleon,"Íleon")
dev.off()

png(filename = "/home/antotartier/data/histones/LAKI/heatmap_LAKI_Músculo.png", width = 9600, height = 10000, res = 1100)
myHeatmap(cts_hist$Músculo,"Riñón")
dev.off()

png(filename = "/home/antotartier/data/histones/LAKI/heatmap_LAKI_Hígado.png", width = 9600, height = 10000, res = 1100)
myHeatmap(cts_hist$Hígado,"Hígado")
dev.off()

png(filename = "/home/antotartier/data/histones/LAKI/heatmap_LAKI_Colon.png", width = 9600, height = 10000, res = 1100)
myHeatmap(cts_hist$Colon,"Colon")
dev.off()

png(filename = "/home/antotartier/data/histones/LAKI/heatmap_LAKI_Corazón.png", width = 9600, height = 10000, res = 1100)
myHeatmap(cts_hist$Corazón,"Corazón")
dev.off()


pdf("/home/antotartier/data/histones/LAKI/heatmap_LAKI_Riñón.pdf", width = 8,height = 8)
myHeatmap(cts_hist$Riñón,"Riñón")
dev.off()

pdf("/home/antotartier/data/histones/LAKI/heatmap_LAKI_Ileon.pdf", width = 8,height = 8)
myHeatmap(cts_hist$Íleon,"Íleon")
dev.off()

pdf("/home/antotartier/data/histones/LAKI/heatmap_LAKI_Músculo.pdf", width = 8,height = 8)
myHeatmap(cts_hist$Músculo,"Riñón")
dev.off()

pdf("/home/antotartier/data/histones/LAKI/heatmap_LAKI_Hígado.pdf", width = 8,height = 8)
myHeatmap(cts_hist$Hígado,"Hígado")
dev.off()

pdf("/home/antotartier/data/histones/LAKI/heatmap_LAKI_Colon.pdf", width = 8,height = 8)
myHeatmap(cts_hist$Colon,"Colon")
dev.off()

pdf("/home/antotartier/data/histones/LAKI/heatmap_LAKI_Corazón.pdf", width = 8,height = 8)
myHeatmap(cts_hist$Corazón,"Corazón")
dev.off()


#-------------------------------------------------------------------------------

#GSEA HISTONAS


#cargamos el gene_set y el resultado DESEq
load("/home/antotartier/data/histones/gene_sets/mouse/gene_set_histonas_raton.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_all.obj")
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")


#generamos el gene set de todos los genes de histonas
histone_gene_set<-list(histonas = gene_set$Ensembl_Accession_ID)


#separamos genes y pseudogenes
histone_genes<-gene_set[!grepl("pseudogene",gene_set$BioTypes) & !grepl("lncRNA",gene_set$BioTypes),]#filtarmos genes
histone_pseudogenes<-gene_set[grepl("pseudogene",gene_set$Marker_Symbo) &  !grepl("lncRNA",gene_set$BioTypes),]#filtramos pseudogenes


#obtenemos los id y los guardamos
gene_id<-histone_genes$Ensembl_Accession_ID
pseudogene_id<-histone_pseudogenes$Ensembl_Accession_ID
save(histone_genes,file = "/home/antotartier/data/histones/gene_sets/mouse/histone_genes.obj")
save(histone_pseudogenes,file = "/home/antotartier/data/histones/gene_sets/mouse/histone_pseudogenes.obj")


#generamos el gene set de genes y pseudogenes
gene_pseudo_gene_set<-list(genes_histonas = gene_id, pseudogenes_histonas = pseudogene_id)


#generamos el rank
rank<-lapply(res_all,get_rank)


#ejecutamos fgsea

#Primero para todas las histonas
fgseaRes <- lapply(rank,fgsea,pathways = histone_gene_set)
fgseaRes_sig<-list()
for(name in names(fgseaRes)){
  x<-fgseaRes[[name]]
  y<-x[x$padj<0.05,]
  fgseaRes_sig[[name]]<-y
}

#Dividiendo en genes y pseudogenes
fgseaRes_groups <- lapply(rank,fgsea,pathways = gene_pseudo_gene_set)
fgseaRes_groups_sig<-list()
for(name in names(fgseaRes_groups)){
  x<-fgseaRes_groups[[name]]
  y<-x[x$padj<0.05,]
  fgseaRes_groups_sig[[name]]<-y
}


#hacemos los plots del enriquecimiento en genes pues en pseudogenes no fue significativo en ningún caso
fgseaRes_groups_sig_df<-lapply(fgseaRes_groups_sig,as.data.frame)
fgseaRes_groups_sig_df<-bind_rows(fgseaRes_groups_sig_df,.id = "Tejido")
GSEA_all_plot<-ggplot(fgseaRes_groups_sig_df,aes(x=NES,y=Tejido))+
  geom_col(aes(fill=padj))+
  ggtitle("GSEA genes de histonas LAKI")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado")



#Análisis en machos y hembras por separado
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_M.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_F.obj")


#generamos los rank
rank_M<-lapply(res_M,get_rank)
rank_F<-lapply(res_F,get_rank)


#ejecutamos fgsea

#Primero para todas las histonas
fgseaRes_M <- lapply(rank_M,fgsea,pathways = histone_gene_set)
fgseaRes_sig_M<-list()
for(name in names(fgseaRes_M)){
  x<-fgseaRes_M[[name]]
  y<-x[x$padj<0.05,]
  fgseaRes_sig_M[[name]]<-y
}

fgseaRes_F <- lapply(rank_F,fgsea,pathways = histone_gene_set)
fgseaRes_sig_F<-list()
for(name in names(fgseaRes_F)){
  x<-fgseaRes_F[[name]]
  y<-x[x$padj<0.05,]
  fgseaRes_sig_F[[name]]<-y
}

#Dividiendo en genes y pseudogenes
fgseaRes_groups_M <- lapply(rank_M,fgsea,pathways = gene_pseudo_gene_set)
fgseaRes_groups_sig_M<-list()
for(name in names(fgseaRes_groups_M)){
  x<-fgseaRes_groups_M[[name]]
  y<-x[x$padj<0.05,]
  fgseaRes_groups_sig_M[[name]]<-y
}

fgseaRes_groups_F <- lapply(rank_F,fgsea,pathways = gene_pseudo_gene_set)
fgseaRes_groups_sig_F<-list()
for(name in names(fgseaRes_groups_F)){
  x<-fgseaRes_groups_F[[name]]
  y<-x[x$padj<0.05,]
  fgseaRes_groups_sig_F[[name]]<-y
}

#hacemos los plots del enriquecimiento en genes pues en pseudogenes no fue significativo en ningún caso
fgseaRes_groups_sig_M_df<-lapply(fgseaRes_groups_sig_M,as.data.frame)
fgseaRes_groups_sig_M_df<-bind_rows(fgseaRes_groups_sig_M_df,.id = "Tejido")
GSEA_M_plot<-ggplot(fgseaRes_groups_sig_M_df,aes(x=NES,y=Tejido))+
  geom_col(aes(fill=padj))+
  ggtitle("GSEA genes de histonas LAKI machos")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado")

fgseaRes_groups_sig_F_df<-lapply(fgseaRes_groups_sig_F,as.data.frame)
fgseaRes_groups_sig_F_df<-bind_rows(fgseaRes_groups_sig_F_df,.id = "Tejido")
GSEA_F_plot<-ggplot(fgseaRes_groups_sig_F_df,aes(x=NES,y=Tejido))+
  geom_col(aes(fill=padj))+
  ggtitle("GSEA genes de histonas LAKI hembras")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado")


#guardamos todos los plots
png(filename = "/home/antotartier/data/histones/LAKI/GSEA_plots_LAKI.png", height =3200 , width =12000 ,res =1100 )
grid.arrange(GSEA_all_plot,GSEA_M_plot,GSEA_F_plot,ncol=3)
dev.off()

ggsave("/home/antotartier/data/histones/LAKI/GSEA_plots_LAKI.pdf",GSEA_all_plot)
ggsave("/home/antotartier/data/histones/LAKI/GSEA_plots_LAKI_M.pdf",GSEA_M_plot)
ggsave("/home/antotartier/data/histones/LAKI/GSEA_plots_LAKI_F.pdf",GSEA_F_plot)


#-------------------------------------------------------------------------------

#ANÁLISIS DE EXPRESIÓN DIFERENCIAL


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


#eliminamos los posibles NA
res_hist<-lapply(res_hist,na.omit)


#eliminamos los pseudogenes
res_hist_nopseudo<-list()
for(name in names(res_hist)){
  x<-res_hist[[name]]
  x<-x[x$type=="gen",]
  res_hist_nopseudo[[name]]<-x
}


#juntamos los datos
res_hist_nopseudo_df<-lapply(res_hist_nopseudo,as.data.frame)
res_hist_nopseudo_df<-bind_rows(res_hist_nopseudo_df,.id = "tejido")


#hacemos los gráficos
family<-myBoxplot(res_hist_nopseudo_df,group = familia,sizesig = 4,pointsize = 1,filter=T,comp_vec=NULL)+facet_wrap(~tejido,scales = "free_y")
reg<-myBoxplot(res_hist_nopseudo_df,group = reg,sizesig = 4,pointsize = 1)+facet_wrap(~tejido,scales = "free_y")

#guardamos los plots 
ggsave("/home/antotartier/data/histones/LAKI/boxplots_family.pdf",family)
ggsave("/home/antotartier/data/histones/LAKI/boxplots_reg.pdf",reg)


#-------------------------------------------------------------------------------

#CALCULAMOS EL LFC MEDIO DE LAS HISTONAS EN LOS DIFERENTES TEJIDOS PARA UTILIZARLO EN LOS ANÁLISIS DE CORRELACIÓN

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
