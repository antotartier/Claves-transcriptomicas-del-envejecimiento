#cargamos las librerias necesarias
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(ggpubr)
library(openxlsx)
library(readr)
library(tidyr)
library(rtracklayer)
library(clusterProfiler)
library(forcats)
library(fgsea)
library(gridExtra)

#definimos el out_path y el path de las matrices de conteo
out_path<-"~/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes"
counts_path<- Sys.glob(paste("~/data/transposones/fibrob_senescentes/TEtranscripts/TEcount", "/*cntTable", sep = ""))

#funciones del script

#selecciona los transposones
get_TE<-function(x){
  TE<-x[grepl(":",rownames(x)),]
  return(TE)
}

#selecciona los genes codificantes
get_codif<-function(x){
  cod<-x[!grepl(":",rownames(x)),]
  return(cod)
}

#permite ajustar los datos para que los valores altos no deformen la gráfica
square_data <- function(x,maxVal,modeFC=T) {
  if(modeFC == T) {
    ifelse(x > maxVal,maxVal,x)
    ifelse(abs(x) > maxVal,-maxVal,x)
  }
  if(modeFC == F){
    ifelse(x < maxVal,maxVal,x)
  }
}

#representación volcano plot
myVolcano<-function(data,biosig,sig,labsize,pointsize,axsize,legsize,axtextsize,topn){
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
                     max.overlaps = 4,
                     size=labsize) +
    scale_colour_manual(values = c("Sobreexpresado"="blue","Infraexpresado"="red","ns"="#3b3b3b")) + 
    scale_x_continuous(breaks = seq(round(-maxLFC),round(maxLFC),1),    
                       limits = c(-maxLFC,maxLFC))+
    theme(plot.title = element_text(hjust = 0.5),axis.title = element_text(size=axsize),axis.text = (element_text(size = axtextsize)),legend.title = element_blank(),legend.text = element_text(size=legsize),axis,legend.position = "top")
  return(p)
}

#obtiene los ID de los genes sobre expresados de forma significativa
get_ID_positivos<-function(x){
  res<-na.omit(x)
  res_positivos<-res[res$log2FoldChange>0 & res$padj<0.05,]
  ID_positivos<-rownames(res_positivos)
  return(ID_positivos)
}

#obtiene los ID de los genes infra expresados de forma significativa
get_ID_negativos<-function(x){
  res<-na.omit(x)
  res_negativos<-res[res$log2FoldChange<0 & res$padj<0.05,]
  ID_negativos<-rownames(res_negativos)
  return(ID_negativos)
}

#genera el universe del ORA
get_universe<-function(x){
  universe<-rownames(x)
  return(universe)
}

#ejecuta el ORA y calculea el enrichment score, la "x" es el ID de interes y la "y" el universe
get_ORA<-function(x,y){
  ORA<-fora(term_to_gene_list,x,y)
  ORA<-mutate(ORA, enrichment= overlap/((size/length(y))*length(x)))#enrichment
  ORA$pathway<-gsub(".gene_name", "", ORA$pathway, fixed = TRUE)#quitamos el ".gene_name"
  return(ORA)
}

#-------------------------------------------------------------------------------

#ANÁLISIS DE LOS CAMBIOS DE EXPRESIÓN ASOCIADOS A LA SENESCENCIA E INDEPENDIENTES A LA LÍNEA CELULAR

#Generamos el input de DESeqDataSetFromMatrix

#Fusionamos los archivos de count de todas las muestras
countData<-read_delim(counts_path, id = "sample", col_names = F, skip = 1,show_col_types = F)
countData<-pivot_wider(countData,names_from = sample, values_from = X2)
colnames(countData)<-c("gene_id", basename(colnames(countData[2:19])))
countData<-data.frame(countData, row.names = 1)

#Generamos el colData
sample_info<-read.xlsx("~/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/sample_info.xlsx")
sample_info<-sample_info[sample_info$ID==c(1:9)|sample_info$ID==c(19:27),1:3]
sample_info$Passage<-mgsub(sample_info$Passage, c("p14","p20","p24","p21","p27","p34"), c("Temprano","Medio","Tardío","Temprano","Medio","Tardío"))
colnames(sample_info)[colnames(sample_info)=="Passage"]<-"Pase"
colData<-sample_info[,2:3]
colnames(colData)[colnames(colData)=="Cells"]<-"Células"
numeric_part<-as.numeric(gsub("[^0-9]", "",counts_path))
colnames_cts<-counts_path[order(numeric_part)]
row.names(colData)<-basename(colnames_cts)
#Convertimos nuestras variables de inetrés en factores y determinamos el orden de los pases
colData$Células<-factor(colData$Células)
colData$Pase<-factor(colData$Pase,levels = c("Temprano","Medio","Tardío"))

#ordenamos las columnas del countData de forma que aparezcan en el mismo orden que las filas en el colData
countData=countData[,basename(colnames_cts)]


#Generamos el DESeqDataSet
DDS_dataset<-DESeqDataSetFromMatrix(countData = countData, colData = colData, design= ~Células + Pase)


#Eliminamos los genes con pocas lecturas, aquellos que no tienen al menos 10 lecturas en el grupo de análisis más pequeño (en este caso el pase)
smallestGroupSize <- 6 #6 muestras por cada pase
keep <- rowSums(counts(DDS_dataset) >= 10) >= smallestGroupSize
DDS_dataset <- DDS_dataset[keep,]


#Ejecutamos la función DESeq, usando el Likehood Ratio Test y el "reduce" para observar los cambios comunes a ambas líneas celulares en senescencia y guardamos el resultado.
DDS_exp<-DESeq(DDS_dataset, test = "LRT", reduced=~Células)
save(DDS_exp,file = paste0(out_path,"/DESeq_conj.obj"))


#Visualizamos las diferencias de expresión de los transposones a nivel general mediante un PCA
vsd <- vst(DDS_exp, blind=FALSE)
save(vsd,file = paste0(out_path,"vsd_conj.obj"))
vsd_TE<-vsd[grepl(":",rownames(vsd)),]

PCA_all_TE<-plotPCA(vsd_TE, intgroup=c("Pase","Células"),returnData =T) 
percentVar_TE<- round(100 * attr(PCA_all_TE, "percentVar"))
PCA_plot_all_TE<-ggplot(PCA_all_TE, aes(PC1, PC2, color=Pase, shape=Células))+
  geom_point(size=3)+
  xlab(paste0("PC1: ",percentVar_TE[1],"% Variabilidad")) +
  ylab(paste0("PC2: ",percentVar_TE[2],"% Variabilidad")) + 
  theme(legend.title = element_text(size = 20),legend.text = element_text(size=15))+
  coord_fixed()

png(filename = paste0(out_path,"/PCA_conj_TE.png"),height = 6000,width = 8000,res = 1100)
PCA_plot_all_TE
dev.off()
ggsave(paste0(out_path,"/PCA_conj_TE.pdf"),PCA_plot_all_TE)


#Genereamos los resultados y los guardamos
res<-results(DDS_exp)
write.csv2(res,file = paste0(out_path,"/DESeq_res_conj.csv"),quote = F)
save(res,file = paste0(out_path,"/DESeq_res_conj.obj"))


#Separamos los transposones y los genes codificantes
TE_only<-get_TE(res)
Codif<-get_codif(res)


#Análisis general de los resultados, comparando la expresión de los genes codificantes con la de los transposones
print("Genes codificantes")
DESeq2:::summary.DESeqResults(Codif,alpha = 0.05)
print("Transposones")
DESeq2:::summary.DESeqResults(TE_only,alpha = 0.05)


#Hacemos un volcano plot
TE_only_sq<-as.data.frame(TE_only) %>%
  mutate(padj, padj=square_data(padj,1e-25,F)) %>%
  mutate(log2FoldChange, aa=square_data(log2FoldChange,4,T))

Volc_conj<-myVolcano(TE_only_sq,biosig = 1,sig = 0.05,labsize = 6,pointsize = 5,axsize = 25,legsize = 30,axtextsize = 20,topn = 10)

png(filename = paste0(out_path,"/Volc_conj.png"),height = 10000,width = 11000,res = 1100)
Volc_conj
dev.off()
ggsave(paste0(out_path,"/Volc_conj.pdf"),Volc_conj,width =10,height = 10)

#HACEMOS UN ORA

#Cargamos los resultados
load(paste0(out_path,"/DESeq_res_conj.obj"))
res_df<-as.data.frame(res)


#Generamos la lista de gene_ID de interés a partir del df de los resultados, en este caso dos: genes sobre e infraexpresados
ID_positivos_conj<-get_ID_positivos(res_df)
ID_negativos_conj<-get_ID_negativos(res_df)


#Generamos el universea partir del df de los resultados, es decir una lista con todos los genes del experimento
universe_conj<-get_universe(res_df)


#creamos el term_to_gene a partir del rmsk
rmsk <- "/data/genomes/GRCh38_gencode/GRCh38_GENCODE_rmsk_TE.gtf"
rmsk <- import(rmsk, format = "gtf")
term_to_gene <- data.frame(
  term = rmsk$family_id,
  gene_name = paste(rmsk$gene_id, rmsk$family_id, rmsk$class_id, sep = ":")
) %>% unique()
#Term-to-gene agrupado en lista para fgsea
term_to_gene_list <-term_to_gene %>%
  group_by(term) %>%
  nest() %>%
  pull(data,name = term)
term_to_gene_list <- unlist(term_to_gene_list,recursive = F)


#Ejecutamos la función fora de fgsea, calculamos el enrichment score y filtramos los resultados significativos
ORA_positivos_conj<-get_ORA(ID_positivos_conj,universe_conj)
ORA_negativos_conj<-get_ORA(ID_negativos_conj,universe_conj)

ORA_positivos_conj_sig<-ORA_positivos_conj[ORA_positivos_conj$padj<0.05,]
ORA_negativos_conj_sig<-ORA_negativos_conj[ORA_negativos_conj$padj<0.05,]


#Representamos los resultados significaticos
ORA_conj<-ggplot(ORA_negativos_conj_sig,aes(x=enrichment,y=fct_reorder(pathway,enrichment),fill=padj)) +
  ylab("FAMILIAS")+
  xlab("ENRIQUECIMIENTO")+
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado") +
  theme(axis.title = element_text(size=16,face = "bold"),legend.title = element_text(size=15),axis.text = element_text(size=12,face = "bold"),legend.text = element_text(size=12))

png(filename = paste0(out_path,"/ORA_conj.png"),height = 6000,width = 7000,res = 1100)
ORA_conj
dev.off()
ggsave(paste0(out_path,"/ORA_conj.pdf"),ORA_conj)


#Finalmente, dado el caracter sorprendente de los resutlados los validamos estudiando la expresión de p21 (CDKN1A), un marcador de senescencia
p21<-plotCounts(DDS_exp,gene = "ENSG00000124762.14",intgroup = c("Pase","Células"),returnData = T)
p21_plot<-ggplot(p21,aes(x=Pase,y=count))+
  geom_point(aes(color=Células),size=3)+
  ylab("Lecturas normalizadas")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("EXPRESIÓN p21")

png(filename = paste0(out_path,"/p21_expr.png"),height = 6000,width = 8000,res = 1100)
p21_plot
dev.off()


#-------------------------------------------------------------------------------

#ANÁLISIS DE LOS CAMBIOS DE EXPRESIÓN ASOCIADOS A LA SENESCENCIA EN CADA LÍNEA CELULAR POR SEPARADO

#Separamos los counts en función de la línea celular
countData_GM05565<-countData[,c(1:9)]
countData_GM00038<-countData[,c(10:18)]

#Generamos el coldata de cada una de las líneas celulares
colData_GM05565<-data.frame(Pase=colData[1:9,2],row.names = colnames(countData_GM05565))
colData_GM00038<-data.frame(Pase=colData[10:18,2], row.names = colnames(countData_GM00038))

colData_GM05565$Pase<-factor(colData_GM05565$Pase)
colData_GM00038$Pase<-factor(colData_GM00038$Pase)


#generamos el DESeqDataSet
DESeqDataSet_GM05565<-DESeqDataSetFromMatrix(countData = countData_GM05565, colData = colData_GM05565, design = ~ Pase)
DESeqDataSet_GM00038<-DESeqDataSetFromMatrix(countData = countData_GM00038, colData = colData_GM00038, design = ~ Pase)


#Eliminamos los genes con pocas lecturas
smallestGroupSize <- 3 #en este caso el número mínimo de muestras por pase es 3 al analizar cada línea por separado
keep_GM05565 <- rowSums(counts(DESeqDataSet_GM05565) >= 10) >= smallestGroupSize
DESeqDataSet_GM05565 <- DESeqDataSet_GM05565[keep_GM05565,]
keep_GM00038 <- rowSums(counts(DESeqDataSet_GM00038) >= 10) >= smallestGroupSize
DESeqDataSet_GM00038 <- DESeqDataSet_GM00038[keep_GM00038,]


#Ejecutamos la función DESeq
DESeq_GM05565<-DESeq(DESeqDataSet_GM05565)
DESeq_GM00038<-DESeq(DESeqDataSet_GM00038)


#Visualizamos las diferencias de expresión de los transposones a nivel general mediante un PCA
vsd_GM05565 <- vst(DESeq_GM05565, blind=FALSE)
vsd_GM00038<- vst(DESeq_GM00038, blind = FALSE)
vsd_GM00038_TE<-vsd_GM00038[grepl(":",rownames(vsd_GM00038)),]
vsd_GM05565_TE<-vsd_GM05565[grepl(":",rownames(vsd_GM05565)),]

PCA_GM05565_TE<-plotPCA(vsd_GM05565_TE, intgroup=("Pase")) + labs(title="PCA GM05565 TE", colour = "Pase") + theme(plot.title = element_text(hjust = 0.5))
PCA_GM00038_TE<-plotPCA(vsd_GM00038_TE, intgroup=("Pase")) + labs(title="PCA GM00038 TE", colour = "Pase") + theme(plot.title = element_text(hjust = 0.5))

png(filename = paste0(out_path,"/PCA_indep_TE.png"),height = 10000,width = 8000,res = 1100)
grid.arrange(PCA_GM00038_TE,PCA_GM05565_TE,ncol=1)
dev.off()


#Ejecutamos los análisis comparativos
res_GM05565_TempranovsMedio<-results(DESeq_GM05565, contrast = c("Pase","Medio","Temprano"))
res_GM05565_MediovsTardío<-results(DESeq_GM05565, contrast = c("Pase","Tardío","Medio"))
res_GM05565_TempranovsTardío<-results(DESeq_GM05565, contrast = c("Pase","Tardío","Temprano"))

res_GM00038_TempranovsMedio<-results(DESeq_GM00038, contrast = c("Pase","Medio","Temprano"))
res_GM00038_MediovsTardío<-results(DESeq_GM00038, contrast = c("Pase","Tardío","Medio"))
res_GM00038_TempranovsTardío<-results(DESeq_GM00038, contrast = c("Pase","Tardío","Temprano"))


#Guardamos los resultados como csv2 y como objeto
res_list<-list(GM05565_Temprano_vs_Medio=res_GM05565_TempranovsMedio,GM05565_Medio_vs_Tardío=res_GM05565_MediovsTardío,GM05565_Temprano_vs_Tardío=res_GM05565_TempranovsTardío,GM00038_Temprano_vs_Medio=res_GM00038_TempranovsMedio,GM00038_Medio_vs_Tardío=res_GM00038_MediovsTardío,GM00038_Temprano_vs_Tardío=res_GM00038_TempranovsTardío)
res_list_df<-lapply(res_list,as.data.frame)
res_indpt<-bind_rows(res_list_df,.id = "experiment")

write.csv2(res_indpt, paste0(out_path,"/DESeq_res_indpte.csv"),quote = F)
save(res_list_df,file = paste0(out_path,"/DESeq_res_indpte.obj"))


#Separamos los transposones y los genes codificantes
TE_list<-lapply(res_list,get_TE)
codif_list<-lapply(res_list,get_codif)


#Análisis general de los resultados, comparando la expresión de los genes codificantes con la de los transposones
for (i in 1:6){
  print(paste(names(res_list[i]), "Genes codificantes",sep = " "))
  print(DESeq2:::summary.DESeqResults(res_list[[i]],alpha = 0.05))
  print(paste(names(TE_list[i]), "Transposones",sep = " "))
  print(DESeq2:::summary.DESeqResults(TE_list[[i]],alpha = 0.05))
}


#Hacemos los volcano plots
TE_squared<-list()
for (i in TE_list){
  TE<-as.data.frame(i)%>%
    mutate(log2FoldChange,aa = square_data(log2FoldChange,4,T)) %>%
    mutate(padj,padj = square_data(padj,1e-25,F))
  TE_squared<-append(TE_squared, list(TE))
}
names(TE_squared)<-names(TE_list)
TE_squared_df<-bind_rows(TE_squared, .id="experiment")

Volc_indep<-myVolcano(data = TE_squared_df,biosig = 1,sig = 0.05,labsize = 2,pointsize = 2,topn = 0,axsize = 10,legsize = 10,axtextsize = 10)+ facet_wrap(vars(factor(experiment, levels = c("GM05565_Temprano_vs_Medio","GM05565_Medio_vs_Tardío","GM05565_Temprano_vs_Tardío","GM00038_Temprano_vs_Medio","GM00038_Medio_vs_Tardío","GM00038_Temprano_vs_Tardío"))))

png(filename = paste0(out_path,"/Volc_indep.png"),height = 8000,width = 14000,res = 1100)
Volc_indep
dev.off()
ggsave(paste0(out_path,"/Volc_indep.pdf"),Volc_indep)


#HACEMOS EL ORA

#Cargamos los resultados de DESeq
load(paste0(out_path,"/DESeq_res_indpte.obj"))


#Obtenemos los conjuntos de interés: genes infra y sobreexpresados
ID_positivos_list<-lapply(res_list_df,get_ID_positivos)
ID_negativos_list<-lapply(res_list_df,get_ID_negativos)


#creamos el universe
universe_list<-lapply(res_list_df,get_universe)


#Ejecutamos la función fora de fgsea, calculamos el enrichment score y filtramos los resultados significativos (el term2gene es el mismo que antes)
ORA_positivos<-mapply(get_ORA, ID_positivos_list, universe_list, SIMPLIFY = FALSE)
ORA_positivos_signif<-lapply(ORA_positivos,function(x){
  filter(x,padj<0.05)
})

ORA_negativos<-mapply(get_ORA,ID_negativos_list,universe_list, SIMPLIFY = FALSE)
ORA_negativos_signif<-lapply(ORA_negativos,function(x){
  filter(x,padj<0.05)
})


#Representamos los resultados significaticos
ORA_negativos_signif_df <- dplyr::bind_rows(ORA_negativos_signif,.id="experiment")
ORA_GM00038<-ggplot(ORA_negativos_signif_df,aes(x=enrichment,y=fct_reorder(pathway,enrichment),fill=padj)) +
  ylab("Familia")+
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "blue", high = "red",name = "P.ajustado") +
  theme(plot.title = element_text(hjust = 0.5),strip.text = element_text(size = 12))+
  facet_wrap(~factor(experiment, c("GM00038_Temprano_vs_Medio","GM00038_Medio_vs_Tardío","GM00038_Temprano_vs_Tardío")))


png(filename = paste0(out_path,"/ORA_GM00038.png"),height = 10000,width = 14000,res = 1100)
ORA_GM00038
dev.off()
ggsave(paste0(out_path,"/ORA_GM00038.pdf"),ORA_GM00038,height = 8,width = 10)


