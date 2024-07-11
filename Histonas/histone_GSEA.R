library(openxlsx)
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)


#cargamos el resultado de DESEq
load("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/DESeq_res_conj.obj")
#cargamos el input para hacer el geneset de las histonas
histones<-read.xlsx("/home/antotartier/data/histones/gene_sets/human/gene_set_histonas.xlsx")

#incluimos en le nombre si es un psudogen y si pertenece a un cluster
histones <- histones %>%
  mutate(symbol = ifelse(grepl("pseudogene", name), paste0(symbol, "_pseudogene"),symbol))
histones <- histones %>%
  mutate(symbol = ifelse(grepl("cluster", name),paste0(symbol,"_cluster"),symbol))
#generamos los id y los nombres
histone_name_id<-histones[,c(2,20)]
histone_id<-noquote(histones$ensembl_gene_id)
#guardamos la lista de id de histonas y nombres de histonas
save(histone_id,file = "/home/antotartier/data/histones/gene_sets/human/histone_id.obj")
save(histone_name_id,file = "/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")
#generamos el gene set
histone_gene_set<-list(histonas = histone_id)

#hacemos un gene set para genes y otro para pseudogenes
histone_genes<-histones[!grepl("pseudogene",histones$name) & is.na(histones$pseudogene.org),]#filtarmos genes
histone_pseudogenes<-histones[grepl("pseudogene",histones$name) | !is.na(histones$pseudogene.org),]#filtramos pseudogenes
#obtenemos los id y los guardamos
gene_id<-noquote(histone_genes$ensembl_gene_id)
pseudogene_id<-noquote(histone_pseudogenes$ensembl_gene_id)
save(histone_genes,file = "/home/antotartier/data/histones/gene_sets/human/histone_genes.obj")
save(histone_pseudogenes,file = "/home/antotartier/data/histones/gene_sets/human/histone_pseudogenes.obj")
#generamos el gene set
gene_pseudo_gene_set<-list(genes_histonas = gene_id, pseudogenes_histonas = pseudogene_id)

#generamos el rank, a partir del df del resultado DESeq
get_rank<-function(x){
  res_gsea<-x[x$baseMean>10,]#eliminamos genes con pocas lecturas
  rownames(res_gsea)<-gsub("\\.[0-9]+$","",rownames(res_gsea))#eliminamos el punto de los ensemble IDs
  res_gsea<-mutate(res_gsea,rank= -log10(pvalue+min(res_gsea$pvalue[res_gsea$pvalue!=0], na.rm= T))*sign(log2FoldChange))
  res_gsea<-na.omit(res_gsea)#eliminamos na y duplicados
  res_gsea<-distinct(res_gsea)
  res_gsea_vector<-res_gsea$rank
  names(res_gsea_vector)<-row.names(res_gsea)
  res_gsea_vector<-sort(res_gsea_vector, decreasing = T)
  return(res_gsea_vector)
}
res_gsea<-as.data.frame(res)
rank<-get_rank(res_gsea)


#ejecutamos fgsea
#Primero para todas las histonas
fgseaRes <- fgsea(pathways = histone_gene_set, 
                  stats    = rank)
fgseaRes_sig<-fgseaRes[fgseaRes$padj<0.05,]#la intersección con el rank es de 80 en vez de 82 (que es el nº de histonas detectadas en el experimento), pues 2 de ellas tienen un basemean<10 y por tanto son eliminadas al generar el rank, concretamente la H3P31 y la H2BC20P
#Dividiendo en genes y pseudogenes
fgseaRes_groups<-fgsea(pathways = gene_pseudo_gene_set, 
                       stats    = rank)
fgseaRes_groups_sig<-fgseaRes_groups[fgseaRes_groups$padj<0.05,]

#hacemos los plots
all<-plotEnrichment(histone_gene_set[["histonas"]], rank) + labs(title="Histonas")#Primero para todas las histonas
#Dividiendo en genes y pseudogenes
genes<-plotEnrichment(gene_pseudo_gene_set[["genes_histonas"]], rank) + labs(title="Genes Histonas")
pseudogenes<-plotEnrichment(gene_pseudo_gene_set[["pseudogenes_histonas"]], rank) + labs(title="Pseudogenes Histonas")

#guardamos todos los plots
png(filename = "/home/antotartier/data/histones/senescencia/plotEnrichment_senescencia.png", height =3200 , width =12000 ,res =1100 )
grid.arrange(all,genes,pseudogenes,ncol=3)
dev.off()

#plot TFM
plot_TFM<-plotEnrichment(gene_pseudo_gene_set[["genes_histonas"]], rank)+xlab("Ranking")+ylab("Enrichment Score")+theme(axis.title = element_text(size = 15),plot.title = element_text(hjust = 0.5,size = 20))+ ggtitle("Genes Histonas")

png(filename = "/home/antotartier/data/histones/senescencia/plotEnrichment_senescencia_hit_gen.png", height =4000 , width =6000 ,res =1100 )
plot_TFM
dev.off()

ggsave("/home/antotartier/data/histones/senescencia/plotEnrichment_senescencia_hit_gen.pdf",plot_TFM)
