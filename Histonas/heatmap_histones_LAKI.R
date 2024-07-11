library(openxlsx)
library(dplyr)


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

#hacemos el heatmap
myHeatmap<-function(data,title){
heatmap(as.matrix(data),
        scale = "row",  
        margins = c(8, 8),      
        cexRow = 0.5,             # Row text size
        cexCol = 1
        )
}

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