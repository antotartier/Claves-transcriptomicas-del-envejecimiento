library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(ggpubr)
library(openxlsx)
library(readr)
library(tidyr)
library(textclean)
library(rtracklayer)
library(forcats)
library(fgsea)
library(patchwork)

#definimos las variables que vamos a utilizar
out_path<-"/home/antotartier/data/transposones/LAKI/resultados_anlaisis/"
counts_path<- Sys.glob("/home/antotartier/data/transposones/LAKI/TEtranscripts/TE_count/*.cntTable")

#Fusionamos los archivos de count de todas las muestras
countData<-read_delim(counts_path, id = "sample", col_names = F, skip = 1)
countData<-pivot_wider(countData,names_from = sample, values_from = X2)
colnames(countData)<-c("gene_id", gsub("Tube_|\\.cntTable","",basename(colnames(countData[2:116]))))
countData<-data.frame(countData, row.names = 1,check.names = F)
countData<-countData[,order(as.numeric(colnames(countData)))]
countData<-countData[,!colnames(countData)%in%c("115","116","117","118","119","120")]#eliminamos las muestras de face


#generamos el coldata 
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")
ausent<-c("Tube_39","Tube_56","Tube_57","Tube_66","Tube_87") #no conseguimos transferirlas
colData<-sample_info[!grepl(paste(ausent,collapse = "|"),sample_info$SampleName)]
colData<-colData %>% mutate(genotipo=str_extract(SampleID,"KO|WT"),
                            sexo=ifelse(grepl("M",SampleID),"M","H"))
colData<-colData[1:109,c(1:2,4:5)]#eliminamos las 6 ultimas muestras que son de FACE
colData$SampleName<-gsub("Tube_","",colData$SampleName)
colData<-data.frame(colData,row.names = 1)
colData$genotipo<-as.factor(colData$genotipo)
colData$genotipo<-relevel(colData$genotipo,ref = "WT")
colData$sexo<-as.factor(colData$sexo)


#hacemos un checkeo de la cadidad de las muestras contando el número de lecturas de los genes codificantes (CDG), pues son los que mejor nos indican la calidad, observamos que la calidad del tubo 86 es claramente peor de lo que la del resto
cdg<-countData[!grepl(":",rownames(countData)),]
quality<-colSums(cdg)
quality<-data.frame(Lecturas_totales=quality,Muestras=rownames(colData))
quality$Muestras <- factor(quality$Muestras, levels = sort(unique(as.numeric(quality$Muestras))))
Q<-ggplot(quality,aes(x=Muestras,y=Lecturas_totales))+
  geom_bar(stat = "identity",width = 1,color="black")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

png(file=paste0(out_path,"total_lectures.png"),width = 6000, height = 4000, res = 1100)
Q
dev.off()
ggsave(paste0(out_path,"total_lectures.pdf"),Q)


#generamos el dataset con todas las muestras juntas 
DDS_dataset<-DESeqDataSetFromMatrix(countData,colData,design = ~sexo + genotipo)

#separamos el data set en función del órgano y lo guardamos en una lista
organs<-unique(DDS_dataset$Tissue)
DDS_dataset_organs<-list()
for(x in organs){
  subset_data<-DDS_dataset[,DDS_dataset$Tissue==x]
  DDS_dataset_organs[[x]]<-subset_data
}
names(DDS_dataset_organs)<-c("Íleon", "Riñón", "Músculo", "Hígado",  "Colon",  "Corazón" )

#filtramos los genes con al menos 10 lecturas en el grupo discreto (KO/WT) de menor tamaño
for(x in seq_along(DDS_dataset_organs)){
  y <-DDS_dataset_organs[[x]]
  smallestGroupSize <- min(c(sum(y$genotipo=="KO"),sum(y$genotipo=="WT")))
  y <- y[rowSums(counts(y) >= 10) >= smallestGroupSize,]
  DDS_dataset_organs[[x]]<-y
}

#realizamos el experimento DESeq
DDS_exp<-lapply(DDS_dataset_organs,DESeq)
DDS_exp_TE<-list()
for(name in names(DDS_exp)){
  x<-DDS_exp[[name]]
  TE<-x[grepl(":",rownames(x)),]
  DDS_exp_TE[[name]]<-TE
}


#determinamos si el sexo influye en la expresión de los transposones
vsd_TE<-lapply(DDS_exp_TE,vst,blind=F)
plotPCA_TE <- function(vsd, intgroup,title,colors) {
  pcaData <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes_string(x = "PC1", y = "PC2", color = intgroup)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors)+
    xlab(paste0("PC1: ", percentVar[1], "% Variación")) +
    ylab(paste0("PC2: ", percentVar[2], "% Variación")) +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))+
    ylim(min(pmin(pcaData$PC1,pcaData$PC2)),max(pmax(pcaData$PC1,pcaData$PC2)))
}


PCA_TE_sex<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "sexo", title = tissue, colors=c("H" = "#AD88C6", "M" = "#F6995C"))
}, vsd_TE, names(vsd_TE), SIMPLIFY = FALSE)
PCA_TE_sex$Colon<-PCA_TE_sex$Colon+ geom_label_repel(data=PCA_TE_sex$Colon$data[PCA_TE_sex$Colon$data$name==86,],aes( label=86))#anotamos la muestra de baja calidad
PCA_TE_plot_sex<-wrap_plots(PCA_TE_sex) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

png(file=paste0(out_path,"PCA_TE_sex.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_sex
dev.off()
ggsave(paste0(out_path,"PCA_TE_sex.pdf"),PCA_TE_plot_sex,height = 8,width = 8)


PCA_TE_genotype<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "genotipo", title = tissue, colors=c("KO" = "#FF90BC", "WT" = "#A5DD9B"))
}, vsd_TE, names(vsd_TE), SIMPLIFY = FALSE)
PCA_TE_genotype$Colon<-PCA_TE_genotype$Colon+ geom_label_repel(data=PCA_TE_genotype$Colon$data[PCA_TE_genotype$Colon$data$name==86,],aes( label=86))#anotamos la muestra de baja calidad
PCA_TE_plot_genotype<-wrap_plots(PCA_TE_genotype) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

png(file=paste0(out_path,"PCA_TE_genotype.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_genotype
dev.off()
ggsave(paste0(out_path,"PCA_TE_genotype.pdf"),PCA_TE_plot_genotype,height = 8,width = 8)

#NOTA: observamos como la muestra 86 cuya calidad es claramente inferior es un outlier

#salvo en riñón e hígado no parece que influya mucho el sexo, sin embargo por si acaso llevaremos a cabo el análisis por separado
DDS_dataset_M<-list()
for(name in names(DDS_dataset_organs)){
  x<-DDS_dataset_organs[[name]]
  M<-x[,x$sexo=="M"]
  design(M)<-formula(~genotipo)#ya no tenemos el sexo como fuente de variación
  DDS_dataset_M[[name]]<-M
}

DDS_dataset_F<-list()
for(name in names(DDS_dataset_organs)){
  x<-DDS_dataset_organs[[name]]
  Fem<-x[,x$sexo=="H"]
  design(Fem)<-formula(~genotipo)
  DDS_dataset_F[[name]]<-Fem
}

#llevamos a cabo el análisis DESeq
DDS_exp_M<-lapply(DDS_dataset_M,DESeq)
DDS_exp_F<-lapply(DDS_dataset_F,DESeq)

#extraigo los TE
DDS_exp_TE_M<-list()
for(name in names(DDS_exp_M)){
  x<-DDS_exp_M[[name]]
  TE<-x[grepl(":",rownames(x)),]
  DDS_exp_TE_M[[name]]<-TE
}

DDS_exp_TE_F<-list()
for(name in names(DDS_exp_F)){
  x<-DDS_exp_F[[name]]
  TE<-x[grepl(":",rownames(x)),]
  DDS_exp_TE_F[[name]]<-TE
}

#normalizamos los datos y hacemos los PCA de cada género
vsd_TE_M<-lapply(DDS_exp_TE_M,vst,blind=F)
vsd_TE_F<-lapply(DDS_exp_TE_F,vst,blind=F)

#PCA machos
PCA_TE_M<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "genotipo", title = tissue, colors=c("KO" = "#FF90BC", "WT" = "#A5DD9B"))
}, vsd_TE_M, names(vsd_TE_M), SIMPLIFY = FALSE)
PCA_TE_M$Colon<-PCA_TE_M$Colon+ geom_label_repel(data=PCA_TE_M$Colon$data[PCA_TE_M$Colon$data$name==86,],aes( label=86))#anotamos la muestra de baja calidad
PCA_TE_plot_M<-wrap_plots(PCA_TE_M) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
png(file=paste0(out_path,"PCA_TE_male.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_M
dev.off()

#PCA hembras
PCA_TE_F<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "genotipo", title = tissue, colors=c("KO" = "#FF90BC", "WT" = "#A5DD9B"))
}, vsd_TE_F, names(vsd_TE_F), SIMPLIFY = FALSE)
PCA_TE_plot_F<-wrap_plots(PCA_TE_F) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
png(file=paste0(out_path,"PCA_TE_female.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_F
dev.off()


#sevobserva claramente que la muestra 86 es un outlier y por tanto se elimina del experimento, eliminandolo del dataset, guardándolo y generando los nuevos experimentos
DDS_dataset<-DDS_dataset[,!colnames(DDS_dataset)==86]
save(DDS_dataset,file = paste0(out_path,"DDS_dataset_all.obj"))
organs<-unique(DDS_dataset$Tissue)
DDS_dataset_organs<-list()
for(x in organs){
  subset_data<-DDS_dataset[,DDS_dataset$Tissue==x]
  DDS_dataset_organs[[x]]<-subset_data
}
names(DDS_dataset_organs)<-c("Íleon", "Riñón", "Músculo", "Hígado",  "Colon",  "Corazón" )
#filtramos los genes con al menos 10 lecturas en el grupo discreto (KO/WT) de menor tamaño
for(x in seq_along(DDS_dataset_organs)){
  y <-DDS_dataset_organs[[x]]
  smallestGroupSize <- min(c(sum(y$genotipo=="KO"),sum(y$genotipo=="WT")))
  y <- y[rowSums(counts(y) >= 10) >= smallestGroupSize,]
  DDS_dataset_organs[[x]]<-y
}
DDS_dataset_M<-list()
for(name in names(DDS_dataset_organs)){
  x<-DDS_dataset_organs[[name]]
  M<-x[,x$sexo=="M"]
  design(M)<-formula(~genotipo)#ya no tenemos el sexo como fuente de variación
  DDS_dataset_M[[name]]<-M
}
DDS_dataset_F<-list()
for(name in names(DDS_dataset_organs)){
  x<-DDS_dataset_organs[[name]]
  Fem<-x[,x$sexo=="H"]
  design(Fem)<-formula(~genotipo)
  DDS_dataset_F[[name]]<-Fem
}
DDS_exp<-lapply(DDS_dataset_organs,DESeq)
DDS_exp_M<-lapply(DDS_dataset_M,DESeq)
DDS_exp_F<-lapply(DDS_dataset_F,DESeq)

#hacemos otra vez el PCA de todos y de los machos pero sin el outlier
DDS_exp_TE<-list()
for(name in names(DDS_exp)){
  x<-DDS_exp[[name]]
  TE<-x[grepl(":",rownames(x)),]
  DDS_exp_TE[[name]]<-TE
}

vsd_TE<-lapply(DDS_exp_TE,vst,blind=F)

PCA_TE_sex<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "sexo", title = tissue, colors=c("H" = "#AD88C6", "M" = "#F6995C"))
}, vsd_TE, names(vsd_TE), SIMPLIFY = FALSE)
PCA_TE_plot_sex<-wrap_plots(PCA_TE_sex) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
png(file=paste0(out_path,"PCA_TE_sex_sinout.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_sex
dev.off()

PCA_TE_genotype<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "genotipo", title = tissue, colors=c("KO" = "#FF90BC", "WT" = "#A5DD9B"))
}, vsd_TE, names(vsd_TE), SIMPLIFY = FALSE)
PCA_TE_plot_genotype<-wrap_plots(PCA_TE_genotype) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
png(file=paste0(out_path,"PCA_TE_genotype_sinout.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_genotype
dev.off()

DDS_exp_TE_M<-list()
for(name in names(DDS_exp_M)){
  x<-DDS_exp_M[[name]]
  TE<-x[grepl(":",rownames(x)),]
  DDS_exp_TE_M[[name]]<-TE
}

DDS_exp_TE_F<-list()
for(name in names(DDS_exp_F)){
  x<-DDS_exp_F[[name]]
  TE<-x[grepl(":",rownames(x)),]
  DDS_exp_TE_F[[name]]<-TE
}

vsd_TE_M<-lapply(DDS_exp_TE_M,vst,blind=F)
vsd_TE_F<-lapply(DDS_exp_TE_F,vst,blind=F)

PCA_TE_M<-mapply(function(vsd, tissue) {
  plotPCA_TE(vsd, intgroup = "genotipo", title = tissue, colors=c("KO" = "#FF90BC", "WT" = "#A5DD9B"))
}, vsd_TE_M, names(vsd_TE_M), SIMPLIFY = FALSE)
PCA_TE_plot_M<-wrap_plots(PCA_TE_M) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
png(file=paste0(out_path,"PCA_TE_male_sinout.png"),width = 8000, height = 8000, res = 1100)
PCA_TE_plot_M
dev.off()

#proseguimos con el análisis DESeq con la función de resultados
res_all<-lapply(DDS_exp,results,contrast=c("genotipo","KO","WT"))
res_M<-lapply(DDS_exp_M,results)
res_F<-lapply(DDS_exp_F,results)

#guardamos los resultados
#guardamos el resultado
save(DDS_exp,file = paste0(out_path,"DDS_exp_all.obj"))
save(res_all,file = paste0(out_path,"res_all.obj"))
save(res_M,file = paste0(out_path,"res_M.obj"))
save(res_F,file = paste0(out_path,"res_F.obj"))

res_all_df<-lapply(res_all,as.data.frame)
write.xlsx(res_all_df,file = paste0(out_path,"res_all.xlsx"),rowNames=TRUE)

#nos quedamos con los transpsosones
get_TE<-function(x){
  TE<-x[grepl(":",rownames(x)),]
  return(TE)
}
TE_all<-lapply(res_all,get_TE)
TE_M<-lapply(res_M,get_TE)
TE_F<-lapply(res_F,get_TE)

#separamos los genes codificantes
get_codif<-function(x){
  cod<-x[!grepl(":",rownames(x)),]
  return(cod)
}
codif_all<-lapply(res_all,get_codif)
codif_M<-lapply(res_M,get_codif)
codif_F<-lapply(res_F,get_codif)

#imprimimos un pequeño resumen de los resultados 
for (i in 1:6){
  #both
  print(paste(names(codif_all[i]), "Coding genes both",sep = " "))
  print(DESeq2:::summary.DESeqResults(codif_all[[i]],alpha = 0.05))
  print(paste(names(TE_all[i]), "Transposons both",sep = " "))
  print(DESeq2:::summary.DESeqResults(TE_all[[i]],alpha = 0.05))
  #males
  print(paste(names(codif_M[i]), "Coding genes males",sep = " "))
  print(DESeq2:::summary.DESeqResults(codif_M[[i]],alpha = 0.05))
  print(paste(names(TE_M[i]), "Transposons males",sep = " "))
  print(DESeq2:::summary.DESeqResults(TE_M[[i]],alpha = 0.05))
  #females
  print(paste(names(codif_F[i]), "Coding genes females",sep = " "))
  print(DESeq2:::summary.DESeqResults(codif_F[[i]],alpha = 0.05))
  print(paste(names(TE_F[i]), "Transposons females",sep = " "))
  print(DESeq2:::summary.DESeqResults(TE_F[[i]],alpha = 0.05))
}

#hacemos los volcano plots de los diferentes tejidos
myVolcano<-function(data,biosig,sig,labsize,pointsize,topn,titlesize,legtextsize){
  data<-na.omit(data)
  data <- as.data.frame(data) %>%
    mutate(Expresión = case_when(log2FoldChange >= biosig & padj <= sig ~ "Sobreexpresado",
                                 log2FoldChange <= -biosig & padj <= sig ~ "Infraexpresado",
                                 TRUE ~ "ns")) 
  # Ensure legend consistency
  data$Expresión <- factor(data$Expresión, levels = c("Sobreexpresado", "Infraexpresado", "ns"))
  data <- data %>%
    mutate(Expresión = fct_relevel(Expresión, "Sobreexpresado", "Infraexpresado"))
  #set axis limits
  maxLFC<-max(abs(data$log2FoldChange))
  #top genes
  top_up<-slice_max(data,n=topn,order_by = log2FoldChange)
  top_up<-top_up[top_up$log2FoldChange>1 & top_up$padj<0.05,]
  top_down<-slice_min(data,n=topn,order_by = log2FoldChange)
  top_down<-top_down[top_down$log2FoldChange<(-1) & top_down$padj<0.05,]
  top_genes<-rbind(top_up,top_down)
  #plot
  p<-ggplot(data = data,
         aes(x = log2FoldChange,
             y = -log10(padj))) + 
    geom_point(aes(colour = Expresión), 
               alpha = 0.2, 
               shape = 16,
               size = pointsize) + 
    geom_hline(yintercept = -log10(sig),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-biosig, biosig),
               linetype = "dashed") +
    geom_label_repel(data = top_genes,     
                     aes(label = rownames(top_genes)),
                     force = 2,
                     nudge_y = 1,
                     max.overlaps = 100000,
                     size=labsize) +
    scale_colour_manual(values = c("Sobreexpresado"="blue","Infraexpresado"="red","ns"="#3b3b3b")) + 
    scale_x_continuous(breaks = seq(round(-maxLFC),round(maxLFC),1),    
                       limits = c(-maxLFC,maxLFC))+
    theme(legend.title = element_blank(),legend.text = element_text(size = legtextsize),legend.position = "top",strip.text = element_text(size = titlesize))
  return(p)
}

#Todos
TE_all_df<-lapply(TE_all, as.data.frame)
TE_all_df<-bind_rows(TE_all_df,.id = "tejido")
Volc_all_plot<-myVolcano(TE_all_df,biosig=1,sig=0.05,labsize=1.5,pointsize=3,topn=0,titlesize = 15,legtextsize = 15)+facet_wrap(~tejido,scales = "free")
ggsave(filename = paste0(out_path,"Volc_TE_all.pdf"),plot = Volc_all_plot,height = 9,width = 8.5)

#machos
TE_M_df<-lapply(TE_M, as.data.frame)
TE_M_df<-bind_rows(TE_M_df,.id = "tejido")
Volc_M_plot<-myVolcano(TE_M_df,biosig=1,sig=0.05,labsize=1.5,pointsize=3,topn=0,titlesize = 15,legtextsize = 15)+facet_wrap(~tejido,scales = "free")

png(file=paste0(out_path,"Volc_TE_male.png"),width = 12000, height = 8000, res = 1100)
Volc_M_plot
dev.off()
ggsave(paste0(out_path,"Volc_TE_male.pdf"),Volc_M_plot,height = 7,width = 8)


#hembras
TE_F_df<-lapply(TE_F, as.data.frame)
TE_F_df<-bind_rows(TE_F_df,.id = "tejido")
Volc_F_plot<-myVolcano(TE_F_df,biosig=1,sig=0.05,labsize=1.5,pointsize=3,topn=0,titlesize = 15,legtextsize = 15)+facet_wrap(~tejido,scales = "free")

png(file=paste0(out_path,"Volc_TE_female.png"),width = 12000, height = 8000, res = 1100)
Volc_F_plot
dev.off()
ggsave(paste0(out_path,"Volc_TE_female.pdf"),Volc_F_plot,height = 7,width = 8)

#--------------------------------------------------------------------------------------------------------------------
  
#para explorar la expresión a nivel de familia llevamos a cabo un gsea
#función que crea el rank a partir del df del resultado DESeq, se usa el pvalor en vez dle ajustado para evitar los empates y se suma el pvalor mínimo para evitar los infinitos
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

#creamos el term to gene
rmsk <- "/data/genomes/GRCm39_gencode/GRCm39_GENCODE_rmsk_TE.gtf"
rmsk <- import(rmsk, format = "gtf")
term_to_gene<- data.frame(
  term = rmsk$family_id,
  gene_name = paste(rmsk$gene_id, rmsk$family_id, rmsk$class_id, sep = ":")
) %>%
  unique()#nos quedamos con los nombres únicos de los transposones

#genramos los rank 
rank_all <- lapply(res_all,get_rank)
rank_M <- lapply(res_M,get_rank)
rank_F <- lapply(res_F,get_rank)


set.seed(10)#creamos una semilla para obtener siempre los mismos resultados
#creamos una función que intersecte el term to gene y el rank y después ejecute la función GSEA
GSEA_transp<-function(x){
  term_to_gene<-term_to_gene[term_to_gene$gene_name %in% names(x),]#solo los trsndposones con lecturas en el experimento
  GSEA<-GSEA(x,TERM2GENE = term_to_gene,seed=T,pvalueCutoff=1)
  return(GSEA)
}
#ejecutamos la función
GSEA_all<-lapply(rank_all,GSEA_transp)
GSEA_M<-lapply(rank_M,GSEA_transp)
GSEA_F<-lapply(rank_F,GSEA_transp)


#filtramos los resultados significativos
GSEA_all_sig<-lapply(GSEA_all,function(df) {
  filter(df, p.adjust < 0.05)
})
GSEA_M_sig<-lapply(GSEA_M,function(df) {
  filter(df, p.adjust < 0.05)
})
GSEA_F_sig<-lapply(GSEA_F,function(df) {
  filter(df, p.adjust < 0.05)
})

#hacemos los plots
myGSEAplot<-function(data,titlesize,legsize){
  ggplot(as.data.frame(data), aes(x = NES, y = fct_reorder(ID,NES), fill = p.adjust)) + geom_bar(stat = "identity", color="black") + labs(x = "NES", y = "FAMILIAS",size=20) + scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado" )+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.title = element_text(size = legsize),strip.text = element_text(size = titlesize))+facet_wrap(~tejido)
}


GSEA_all_sig<-lapply(GSEA_all_sig,as.data.frame)
GSEA_all_sig<-bind_rows(GSEA_all_sig,.id="tejido")
GSEA_all_plot<-myGSEAplot(GSEA_all_sig,titlesize = 15,legsize = 15)
png(file=paste0(out_path,"GSEA_TE_all.png"),width = 10000, height = 8000, res = 1100)
GSEA_all_plot
dev.off()
ggsave(paste0(out_path,"GSEA_TE_all.pdf"),GSEA_all_plot)

GSEA_M_sig<-lapply(GSEA_M_sig,as.data.frame)
GSEA_M_sig<-bind_rows(GSEA_M_sig,.id="tejido")
GSEA_M_plot<-myGSEAplot(GSEA_M_sig,titlesize = 15,legsize = 15)
png(file=paste0(out_path,"GSEA_TE_male.png"),width = 10000, height = 8000, res = 1100)
GSEA_M_plot
dev.off()
ggsave(paste0(out_path,"GSEA_TE_M.pdf"),GSEA_M_plot,height = 8,width = 8)

GSEA_F_sig<-lapply(GSEA_F_sig,as.data.frame)
GSEA_F_sig<-bind_rows(GSEA_F_sig,.id="tejido")
GSEA_F_plot<-myGSEAplot(GSEA_F_sig,titlesize = 15,legsize = 15)
png(file=paste0(out_path,"GSEA_TE_female.png"),width = 10000, height = 8000, res = 1100)
GSEA_F_plot
dev.off()
ggsave(paste0(out_path,"GSEA_TE_F.pdf"),GSEA_F_plot,height = 8,width = 8)
