#cargamos las librerías necesarias
library(openxlsx)
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggsignif)

#funciones

#genera el rank, a partir del df del resultado DESeq
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

#boxplots con prueba T
myBoxplot<-function(data,group,titulo,sizesig=0,heightsig=0,axsize,axtextsize,filter=F,comp_vec){
  group <- enquo(group)
  group_col <- ensym(group)  # Capture the group column symbolically
  
  generate_comparisons <- function(levels) {
    combn(levels, 2, simplify = FALSE)
  }#genera todas las posibles comparaciones por pares a partir de los niveles
  
  levels_list <- unique(data[[as.character(group_col)]])  # Extract unique levels
  comparisons <- generate_comparisons(levels_list)
  
  plot<-ggplot(data, aes(x=!!group, y=log2FoldChange, fill=!!group)) +
    geom_boxplot(alpha = 0.80, coef=0,outlier.shape = NA) +
    stat_summary(fun.data = mean_sd, geom = "errorbar", width=0.25)+
    geom_point(aes(fill = !!group), size = 2, shape = 21, position = position_jitterdodge()) +
    ggtitle(titulo)+
    ylab(expression("log"[2]*" Fold Change"))+
    theme(axis.title = element_text(size = axsize),text = element_text(size = axtextsize,face = "bold"),
          axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, hjust = 0.5))
  if(!filter){
    plot<-plot + geom_signif(comparisons = comparisons, test = "t.test",y_position = heightsig,step_increase = 0.1,
                             textsize = sizesig)}
  else{
    plot<-plot +geom_signif(comparisons = list(comp_vec), test = "t.test",y_position = heightsig,step_increase = 0.1,
                            textsize = sizesig)}#este condicional permite quedarse solo con ciertas comparaciones
  
  return(plot)
}

#-------------------------------------------------------------------------------

#GSEA DE LAS HISTONAS

#cargamos el resultado de DESEq
load("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/DESeq_res_conj.obj")


#cargamos el input para hacer el geneset de las histonas, generado previamente a partir de la lista de genes de HGNC
histones<-read.xlsx("/home/antotartier/data/histones/gene_sets/human/gene_set_histonas.xlsx")


#incluimos en le nombre de los genes de histonas si es un psudogen y si pertenece a un cluster
histones <- histones %>%
  mutate(symbol = ifelse(grepl("pseudogene", name), paste0(symbol, "_pseudogene"),symbol))
histones <- histones %>%
  mutate(symbol = ifelse(grepl("cluster", name),paste0(symbol,"_cluster"),symbol))


#generamos los id y los nombres y los guardamos
histone_name_id<-histones[,c(2,20)]
histone_id<-noquote(histones$ensembl_gene_id)

save(histone_id,file = "/home/antotartier/data/histones/gene_sets/human/histone_id.obj")
save(histone_name_id,file = "/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")


#generamos el gene set
histone_gene_set<-list(histonas = histone_id)

#hacemos un gene set para genes y otro para pseudogenes
histone_genes<-histones[!grepl("pseudogene",histones$name) & is.na(histones$pseudogene.org),]#filtarmos genes
histone_pseudogenes<-histones[grepl("pseudogene",histones$name) | !is.na(histones$pseudogene.org),]#filtramos pseudogenes


#obtenemos los id de genes y pseudogenes y los guardamos
gene_id<-noquote(histone_genes$ensembl_gene_id)
pseudogene_id<-noquote(histone_pseudogenes$ensembl_gene_id)
save(histone_genes,file = "/home/antotartier/data/histones/gene_sets/human/histone_genes.obj")
save(histone_pseudogenes,file = "/home/antotartier/data/histones/gene_sets/human/histone_pseudogenes.obj")


#generamos el gene set de genes y pseudogenes
gene_pseudo_gene_set<-list(genes_histonas = gene_id, pseudogenes_histonas = pseudogene_id)

#generamos el rank
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

#Primero para todas las histonas
all<-plotEnrichment(histone_gene_set[["histonas"]], rank) + labs(title="Histonas")

#Dividiendo en genes y pseudogenes
genes<-plotEnrichment(gene_pseudo_gene_set[["genes_histonas"]], rank)+xlab("Ranking")+ylab("Enrichment Score")+theme(axis.title = element_text(size = 15),plot.title = element_text(hjust = 0.5,size = 20))+ ggtitle("Genes Histonas")
pseudogenes<-plotEnrichment(gene_pseudo_gene_set[["pseudogenes_histonas"]], rank) + labs(title="Pseudogenes Histonas")

#guardamos todos los plots
png(filename = "/home/antotartier/data/histones/senescencia/plotEnrichment_senescencia.png", height =3200 , width =12000 ,res =1100 )
grid.arrange(all,genes,pseudogenes,ncol=3)
dev.off()

png(filename = "/home/antotartier/data/histones/senescencia/plotEnrichment_senescencia_hit_gen.png", height =4000 , width =6000 ,res =1100 )
genes
dev.off()
ggsave("/home/antotartier/data/histones/senescencia/plotEnrichment_senescencia_hit_gen.pdf",genes)


#-------------------------------------------------------------------------------

#HEATMAP DE LAS HISTONAS


#cargamos el vsd de DESeq y los id de las histonas y los nombres
load("~/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/vsd_conj.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")
sample_info<-read.xlsx("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/sample_info.xlsx")


#sacamos la matriz de counts normalizados y nos quedamos con los de las histonas
counts<-assay(vsd)
rownames(counts)<-gsub("\\.[0-9]+$","",rownames(counts))#eliminamos el punto de los ensemble IDs
counts_hist<-counts[rownames(counts) %in% histone_name_id$ensembl_gene_id,c(1:3,7:12,16:18)]


#sustituimos los id por los nombres
id_to_name<-setNames(histone_name_id$symbol,histone_name_id$ensembl_gene_id)
rownames(counts_hist)<-id_to_name[rownames(counts_hist)]


#sustituimos en las columnas los tubos por la condición
sample_info<-sample_info[c(1:3,7:9,19:21,25:27),]
sample_info <- sample_info %>%
  mutate(
    condition = case_when(
      Passage == "p14" ~ paste0(Cells, "_Temprano"),
      Passage == "p20" ~ paste0(Cells, "_Medio"),
      Passage == "p24" ~ paste0(Cells, "_Tardío"),
      Passage == "p21" ~ paste0(Cells, "_Temprano"),
      Passage == "p27" ~ paste0(Cells, "_Medio"),
      Passage == "p34" ~ paste0(Cells, "_Tardío"),
      TRUE ~ NA_character_
    ),
    ID = paste0("Tube_",ID,".cntTable")
  )#creamos la columna condition
tube_to_condition<-setNames(sample_info$condition,sample_info$ID)
colnames(counts_hist)<-tube_to_condition[colnames(counts_hist)]


#hacemos el heatmap
heatplot<-heatmap(counts_hist,
                  scale = "row",  
                  margins = c(8, 8),      
                  cexRow = 0.5,             # Row text size
                  cexCol = 1.5)

png(filename = "/home/antotartier/data/histones/senescencia/heatmap_senescencia.png", width = 9600, height = 12800, res = 1100)
heatplot
dev.off()

pdf(file =  "/home/antotartier/data/histones/senescencia/heatmap_senescencia.pdf",height = 10,width = 10)
heatmap(counts_hist,
        scale = "row",  
        margins = c(8, 8),      
        cexRow = 0.5,             # Row text size
        cexCol = 1.5)
dev.off()


#-------------------------------------------------------------------------------

#ANÁLISIS DE EXPRESIÓN DIFERENCIAL

#cargamos el resultado de DESeq y los id y los nombres de las histonas, los genes y los pseudogenes
load("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/DESeq_res_conj.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_genes.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_pseudogenes.obj")


#lo convertimos en un dataframe
res_df<-as.data.frame(res)
rownames(res_df)<-gsub("\\.[0-9]+$","",rownames(res_df))#eliminamos el punto de los ensemble IDs


#nos quedamos con las histonas
res_hist<-res_df[rownames(res_df)%in%histone_name_id$ensembl_gene_id,]


#sustituimos los id por los nombres
id_to_name<-setNames(histone_name_id$symbol,histone_name_id$ensembl_gene_id)
rownames(res_hist)<-id_to_name[rownames(res_hist)]


#hacemos una lista de los nombres de las histonas cluster y no cluster
cluster<-histone_name_id$symbol[grepl("cluster",histone_name_id$symbol)]
cluster<-cluster[cluster %in% histone_genes$symbol]
no_cluster<-histone_name_id$symbol[!grepl("cluster",histone_name_id$symbol)]
no_cluster<-no_cluster[!(no_cluster %in% histone_pseudogenes$symbol)]


#hacemos una lista de los nombres de las histonas por familia
H1<-histone_genes$symbol[grepl("H1",histone_genes$symbol)]
H2A<-histone_genes$symbol[grepl("H2A",histone_genes$symbol)]
H2B<-histone_genes$symbol[grepl("H2B",histone_genes$symbol)]
H3<-histone_genes$symbol[grepl("H3",histone_genes$symbol)]
H4<-histone_genes$symbol[grepl("H4",histone_genes$symbol)]


#definimos un patrón para seleccionar los genes dependientes de la replicación
dpts<-c("cluster","H2AX","H4C16")


#clasificamos los genes
res_hist<-res_hist %>%
  mutate(type=ifelse(rownames(res_hist)%in%histone_pseudogenes$symbol, "pseudogen",
                     ifelse(rownames(res_hist)%in%histone_genes$symbol,"gen", "unknown")),
         clustering=ifelse(rownames(res_hist)%in%cluster, "cluster",
                           ifelse(rownames(res_hist)%in%no_cluster,"No cluster", "unknown")),
         familia=ifelse(rownames(res_hist)%in%H1, "H1",
                        ifelse(rownames(res_hist)%in%H2A,"H2A",
                               ifelse(rownames(res_hist)%in%H2B,"H2B",
                                      ifelse(rownames(res_hist)%in%H3,"H3",
                                             ifelse(rownames(res_hist)%in%H4,"H4", "Unknown"))))),
         reg=ifelse(grepl(paste(dpts, collapse = "|"),rownames(res_hist)) & !grepl("H3-4",rownames(res_hist)),"HDR",
                    ifelse(!grepl(paste(dpts, collapse = "|"),rownames(res_hist)) | grepl("H3-4",rownames(res_hist)),"HIR"
                           ,"unknown")))


#eliminamos los pseudogenes
res_hist_nopseudo<-res_hist[res_hist$type=="gen",]


#hacemos los plots
reg<-myBoxplot(data = res_hist_nopseudo,group = reg,titulo = NULL,sizesig = 4,heightsig = 0.4,axsize = 18,axtextsize = 18)

png(filename = "/home/antotartier/data/histones/senescencia/boxplot_reg.png", height =5000 , width =5000 ,res =1100 )
reg
dev.off()

fam<-myBoxplot(data = res_hist_nopseudo,group = familia,titulo = NULL,axsize = 18,axtextsize = 18,filter = T,comp_vec = NULL)

png(filename = "/home/antotartier/data/histones/senescencia/boxplot_fam.png", height =5000 , width =8000 ,res =1100 )
fam
dev.off()

ggsave("/home/antotartier/data/histones/senescencia/boxplot_reg.pdf",reg)
ggsave("/home/antotartier/data/histones/senescencia/boxplot_fam.pdf",fam)