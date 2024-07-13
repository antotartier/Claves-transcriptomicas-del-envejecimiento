library(dplyr)
library(DESeq2)
library(ggplot2)
library(fgsea)

#definimos el out_path
out_path<-"/home/antotartier/data/velocidad_transcrip/LAKI/regulators/"

#hacemos el GSEA en los datos de LAKI
#cargamos los gene sets de ratón
reg_pos<-fread("/home/antotartier/data/velocidad_transcrip/LAKI/regulators/GOBP_POSITIVE_REGULATION_OF_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION.v2023.2.Mm.tsv")
reg_pos_pol2<-fread("/home/antotartier/data/velocidad_transcrip/LAKI/regulators/GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_ELONGATION_BY_RNA_POLYMERASE_II.v2023.2.Mm.tsv")
reg_neg<-fread("/home/antotartier/data/velocidad_transcrip/LAKI/regulators/GOBP_NEGATIVE_REGULATION_OF_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION.v2023.2.Mm.tsv")
reg_neg_pol2<-fread("/home/antotartier/data/velocidad_transcrip/LAKI/regulators/GOBP_NEGATIVE_REGULATION_OF_TRANSCRIPTION_ELONGATION_BY_RNA_POLYMERASE_II.v2023.2.Mm.tsv")

#cargamos el gtf
gtf_mm<-import("/data/genomes/GRCm39_gencode/gencode.vM32.primary_assembly.annotation.gtf")
gtf_mm<-as.data.frame(gtf_mm)
gtf_mm<-gtf_mm[gtf_mm$type=="gene",]
colnames(gtf_mm)[colnames(gtf_mm)=="gene_name"]<-"GENE_SYMBOLS"

#cargamos los resultados conjuntos de machos y hembras
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_all.obj")

#merge de los gene_set y el gtf
reg_pos<-merge(reg_pos,gtf_mm,by="GENE_SYMBOLS")
reg_pos_pol2<-merge(reg_pos_pol2,gtf_mm,by="GENE_SYMBOLS")
reg_neg<-merge(reg_neg,gtf_mm,by="GENE_SYMBOLS")
reg_neg_pol2<-merge(reg_neg_pol2,gtf_mm,by="GENE_SYMBOLS")

#hacemos el gene term to gene
term2gene<-list(Positivos=reg_pos$gene_id,Negativos=reg_neg$gene_id,Positivos_Pol2=reg_pos_pol2$gene_id,Negativos_Pol2=reg_neg_pol2$gene_id)

#hacmeos los rank
get_rank<-function(x){
  res_gsea<-as.data.frame(x)
  res_gsea<-res_gsea[res_gsea$baseMean>10,]#eliminamos genes con pocas lecturas
  res_gsea<-mutate(res_gsea,rank= -log10(pvalue+min(res_gsea$pvalue[res_gsea$pvalue!=0], na.rm= T))*sign(log2FoldChange))
  res_gsea<-na.omit(res_gsea)#eliminamos na y duplicados
  res_gsea<-distinct(res_gsea)
  res_gsea_vector<-res_gsea$rank
  names(res_gsea_vector)<-row.names(res_gsea)
  res_gsea_vector<-sort(res_gsea_vector, decreasing = T)
  return(res_gsea_vector)
}
rank_list<-lapply(res_all,get_rank)

#hacemos el GSEA, no sale nada significativo
GSEA_list<-lapply(rank_list,function(x){fgsea(pathways = term2gene,stats = x)})
GSEA_list_sig<-list()
for(name in names(GSEA_list)){
  x<-GSEA_list[[name]]
  y<-x[x$padj<0.05]
  GSEA_list_sig[[name]]<-y
}

#--------------------------------------------------------------------------------------------------------------------
#probamos en machos y hembras por separado

#cargamos los resultados de machos y hembras
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_M.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_F.obj")

#hacemos los ranking
rank_M<-lapply(res_M,get_rank)
rank_F<-lapply(res_F,get_rank)

#hacemos el GSEA, no sale nada significativo
GSEA_M<-lapply(rank_M,function(x){fgsea(pathways = term2gene,stats = x)})
GSEA_F<-lapply(rank_F,function(x){fgsea(pathways = term2gene,stats = x)})
GSEA_M_sig<-list()
for(name in names(GSEA_M)){
  x<-GSEA_M[[name]]
  y<-x[x$padj<0.05]
  GSEA_M_sig[[name]]<-y
}
GSEA_F_sig<-list()
for(name in names(GSEA_F)){
  x<-GSEA_F[[name]]
  y<-x[x$padj<0.05]
  GSEA_F_sig[[name]]<-y
}

#----------------------------------------------------------------------------------------------------------------------
#Miramos que reguladores están alterados de forma estadísticamente significativa

#cargamos el resultado conjunto de machos y hembras
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_all.obj")

#filtramos los resultados y extraemos los de los reguladores
for(name in names(res_all)){
x<-res_all[[name]]
y<-na.omit(x)
sig<-y[y$padj<0.05,]
res_all[[name]]<-sig
}

reg_pos_res<-list()
for(name in names(res_all)){
  x<-res_all[[name]]
  y<-x[rownames(x) %in% reg_pos$gene_id,]
  reg_pos_res[[name]]<-y
}

reg_pos_pol2_res<-list()
for(name in names(res_all)){
  x<-res_all[[name]]
  y<-x[rownames(x) %in% reg_pos_pol2$gene_id,]
  reg_pos_pol2_res[[name]]<-y
}

reg_neg_res<-list()
for(name in names(res_all)){
  x<-res_all[[name]]
  y<-x[rownames(x) %in% reg_neg$gene_id,]
  reg_neg_res[[name]]<-y
}

reg_neg_pol2_res<-list()
for(name in names(res_all)){
  x<-res_all[[name]]
  y<-x[rownames(x) %in% reg_neg_pol2$gene_id,]
  reg_neg_pol2_res[[name]]<-y
}

#unimos los resultados y sustituimos el emsemble ID por el gene symbol
reg_pos_res_df<-lapply(reg_pos_res,as.data.frame)
id_to_symbol_pos<-setNames(reg_pos$GENE_SYMBOLS,reg_pos$gene_id)
for(name in names(reg_pos_res_df)){
x<-reg_pos_res_df[[name]]
rownames(x)<-id_to_symbol_pos[rownames(x)]
reg_pos_res_df[[name]]<-x
}
reg_pos_res_df<-bind_rows(reg_pos_res_df,.id = "órgano")

reg_pos_pol2_res_df<-lapply(reg_pos_pol2_res,as.data.frame)
id_to_symbol_pos_pol2<-setNames(reg_pos_pol2$GENE_SYMBOLS,reg_pos_pol2$gene_id)
for(name in names(reg_pos_pol2_res_df)){
  x<-reg_pos_pol2_res_df[[name]]
  rownames(x)<-id_to_symbol_pos_pol2[rownames(x)]
  reg_pos_pol2_res_df[[name]]<-x
}
reg_pos_pol2_res_df<-bind_rows(reg_pos_pol2_res_df,.id = "órgano")

reg_neg_res_df<-lapply(reg_neg_res,as.data.frame)
id_to_symbol_neg<-setNames(reg_neg$GENE_SYMBOLS,reg_neg$gene_id)
for(name in names(reg_neg_res_df)){
  x<-reg_neg_res_df[[name]]
  rownames(x)<-id_to_symbol_neg[rownames(x)]
  reg_neg_res_df[[name]]<-x
}
reg_neg_res_df<-bind_rows(reg_neg_res_df,.id = "órgano")

reg_neg_pol2_res_df<-lapply(reg_neg_pol2_res,as.data.frame)
id_to_symbol_neg_pol2<-setNames(reg_neg_pol2$GENE_SYMBOLS,reg_neg_pol2$gene_id)
for(name in names(reg_neg_pol2_res_df)){
  x<-reg_neg_pol2_res_df[[name]]
  rownames(x)<-id_to_symbol_neg_pol2[rownames(x)]
  reg_neg_pol2_res_df[[name]]<-x
}
reg_neg_pol2_res_df<-bind_rows(reg_neg_pol2_res_df,.id = "órgano")

#guardamos los datos en excel
res<-list(Positivos=reg_pos_res_df,Positivos_pol2=reg_pos_pol2_res_df,Negativos=reg_neg_res_df,Negativos_pol2=reg_neg_pol2_res_df)
write.xlsx(x = res,file = paste0(out_path,"Results.xlsx"),row.names=TRUE)

#exploramos la expresión de leo1
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/DDS_exp_all.obj")
Leo1<-lapply(DDS_exp,plotCounts,gene="ENSMUSG00000042487.7",intgroup=(c("sexo","genotipo")),returnData=T)
Leo1<-bind_rows(Leo1,.id = "tejido")

Leo1_plot<-ggplot(Leo1,aes(y=count,x=genotipo,fill=genotipo))+
geom_boxplot(alpha = 0.80, coef=0,outlier.shape = NA)+
stat_summary(fun.data = mean_sd, geom = "errorbar", width=0.25)+
geom_point(aes(color=genotipo),size = 2, shape = 21, ,color="black",position = position_jitterdodge())+
geom_signif(comparisons = list(c("KO","WT")),test = "t.test")+
ylab("Lecturas normalizadas")+
facet_wrap(~tejido,scales = "free_y")

ggsave(paste0(out_path,"Leo1_plot.pdf"),Leo1_plot)
