#hacemos el heatmap

library(openxlsx)
library(dplyr)

#definimos el out path
out_path<-"/home/antotartier/data/histones/fib_HGPS/"

#cargamos el vsd de DESeq y los id de las histonas y los nombres
load("/home/antotartier/data/transposones/fib_HGPS/results/TE_count/vsd_HGPS.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")
sample_info<-fread("/home/antotartier/data/fib_HGPS/FASTQ/Sra_meta.tsv")
sample_info$meta<-paste(sample_info$Genotipo,sample_info$Sexo,sample_info$Edad,sep = "_")

#sacamos la matriz de counts normalizados y nos quedamos con los de las histonas
counts<-assay(vsd)
rownames(counts)<-gsub("\\.[0-9]+$","",rownames(counts))#eliminamos el punto de los ensemble IDs
counts_hist<-counts[rownames(counts) %in% histone_name_id$ensembl_gene_id,]

#sustituimos los id por los nombres
id_to_name<-setNames(histone_name_id$symbol,histone_name_id$ensembl_gene_id)
rownames(counts_hist)<-id_to_name[rownames(counts_hist)]

#sustituimos en las columnas el nombre de la muestra por la condición
tube_to_condition<-setNames(sample_info$meta,sample_info$SRA)
colnames(counts_hist)<-tube_to_condition[colnames(counts_hist)]

#hacemos el heatmap
png(filename = paste0(out_path,"heatmap_fib_HGPS.png"), width = 9600, height = 12800, res = 1100)
heatmap(counts_hist,
        scale = "row",  
        margins = c(8, 8),      
        cexRow = 0.5,             # Row text size
        cexCol = 1)           # Column text size
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------
#hacemos el GSEA

#cargamos el resultado de DESEq
load("/home/antotartier/data/transposones/fib_HGPS/results/TE_count/res_fib_HGPS.obj")

#generamos el gene set
histone_gene_set<-list(histonas = histone_name_id$ensembl_gene_id)

#cargamos la lista de genes y  pseudogenes
load("/home/antotartier/data/histones/gene_sets/human/histone_genes.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_pseudogenes.obj")
#generamos el gene set
gene_pseudo_gene_set<-list(genes_histonas = histone_genes$ensembl_gene_id, pseudogenes_histonas = histone_pseudogenes$ensembl_gene_id)

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
fgseaRes_sig<-fgseaRes[fgseaRes$padj<0.05,]

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
png(filename = paste0(out_path,"GSEA_fib_HGPS.png"), height =3200 , width =12000 ,res =1100 )
grid.arrange(all,genes,pseudogenes,ncol=3)
dev.off()

ggsave(paste0(out_path,"GSEA_fib_HGPS.pdf"),genes)
