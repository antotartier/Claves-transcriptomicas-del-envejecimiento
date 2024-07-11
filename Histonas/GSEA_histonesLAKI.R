library(openxlsx)
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(DESeq2)

#cargamos el gene_set y el resultado DESEq
load("/home/antotartier/data/histones/gene_sets/mouse/gene_set_histonas_raton.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_all.obj")
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")


#generamos el gene set
histone_gene_set<-list(histonas = gene_set$Ensembl_Accession_ID)

#hacemos un gene set para genes y otro para pseudogenes
histone_genes<-gene_set[!grepl("pseudogene",gene_set$BioTypes) & !grepl("lncRNA",gene_set$BioTypes),]#filtarmos genes
histone_pseudogenes<-gene_set[grepl("pseudogene",gene_set$Marker_Symbo) &  !grepl("lncRNA",gene_set$BioTypes),]#filtramos pseudogenes
#obtenemos los id y los guardamos
gene_id<-histone_genes$Ensembl_Accession_ID
pseudogene_id<-histone_pseudogenes$Ensembl_Accession_ID
save(histone_genes,file = "/home/antotartier/data/histones/gene_sets/mouse/histone_genes.obj")
save(histone_pseudogenes,file = "/home/antotartier/data/histones/gene_sets/mouse/histone_pseudogenes.obj")
#generamos el gene set
gene_pseudo_gene_set<-list(genes_histonas = gene_id, pseudogenes_histonas = pseudogene_id)

#generamos el rank, a partir del df del resultado DESeq
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


#hacemos los plots
fgseaRes_groups_sig_df<-lapply(fgseaRes_groups_sig,as.data.frame)
fgseaRes_groups_sig_df<-bind_rows(fgseaRes_groups_sig_df,.id = "Tejido")
GSEA_all_plot<-ggplot(fgseaRes_groups_sig_df,aes(x=NES,y=Tejido))+
  geom_col(aes(fill=padj))+
  ggtitle("GSEA genes de histonas LAKI")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado")

#---------------------------------------------------------------------------------------------------------------------
#Análisis en machos y hembras por separado
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_M.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_F.obj")
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

#hacemos los plots
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
