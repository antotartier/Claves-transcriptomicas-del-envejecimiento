library(stringi)
library(ggplot2)
library(DESeq2)
library(ggpmisc)

#definimos el out_path
out_path<-"/home/antotartier/data/transposones/LAKI/resultados_anlaisis/"

#cargamos el dds dataset de todas las muestras, el gene set de histonas y el sample info
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/DDS_dataset_all.obj")
load("/home/antotartier/data/histones/gene_sets/mouse/gene_set_histonas_raton.obj")
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")

#eliminamos del sample info las muestras de FACE, las ausentes y el outlier
elim<-c("Tube_39","Tube_56","Tube_57","Tube_66","Tube_86","Tube_87","Tube_115","Tube_116","Tube_117","Tube_118","Tube_119","Tube_120")
sample_info<-sample_info[!grepl(paste0(elim,collapse = "|"),SampleName),]

#seaparmos las variables en columnas
sample_info<-sample_info %>% mutate(genotipo=str_extract(SampleID,"KO|WT"),
                            sexo=ifelse(grepl("M",SampleID),"M","H"))

#cambiamos el nombre del tejido a español
sample_info$Tissue<-stri_replace_all_regex(sample_info$Tissue, pattern = c("Ileon","Kidney","Muscle","Liver","Colon"  ,"Heart"),replacement=c("Íleon","Riñón","Músculo","Hígado","Colon","Corazón"),vectorize_all = F)

#obtenemos los count normalizaddos
DDS_exp_all<-DESeq(DDS_dataset)
vsd_all<-vst(DDS_exp_all)
cts_all<-assay(vsd_all)

#eliminamos la version del ensenmble ID
rownames(cts_all)<-gsub("\\.[0-9]+$","",rownames(cts_all))

#obtenemos los counts de los transposones y de las histonas
cts_TE<-cts_all[grepl(":",rownames(cts_all)),]
cts_hist<-cts_all[rownames(cts_all) %in% gene_set$Ensembl_Accession_ID,]

#calculamos los counts medios por muestra
cts_mean<-data.frame(TE_mean=colMeans(cts_TE),Hist_mean=colMeans(cts_hist),sexo=sample_info$sexo,genotipo=sample_info$genotipo,tejido=sample_info$Tissue,row.names = colnames(cts_TE))

#hacemos el scatter plot para ver si existe correlación entre los counts de histonas y los de transposones
correl<-ggplot(cts_mean,aes(x=TE_mean,y=Hist_mean))+
  geom_point(aes(color=genotipo))+
  geom_smooth(method = lm,)+
  stat_poly_eq(size=3)+
  labs(x="Lecturas normalizadas de transposones",y="Lecturas normalizadas de histonas")+
  facet_wrap(~tejido)
  
#---------------------------------------------------------------------------------------------------------------------
#hacemos lo mismo pero solo con los genes con lecturas consistentes (nZ)
smallestGroupSize <- min(c(sum(DDS_dataset$genotipo=="KO"),sum(DDS_dataset$genotipo=="WT")))
DDS_dataset_nZ<-DDS_dataset[rowSums(counts(DDS_dataset) >= 10) >= smallestGroupSize,]

#obtenemos los count normalizaddos
DDS_exp_all_nZ<-DESeq(DDS_dataset_nZ)
vsd_all_nZ<-vst(DDS_exp_all_nZ)
cts_all_nZ<-assay(vsd_all_nZ)

#eliminamos la version del ensenmble ID
rownames(cts_all_nZ)<-gsub("\\.[0-9]+$","",rownames(cts_all_nZ))

#obtenemos los counts de los transposones y de las histonas
cts_TE_nZ<-cts_all_nZ[grepl(":",rownames(cts_all_nZ)),]
cts_hist_nZ<-cts_all_nZ[rownames(cts_all_nZ) %in% gene_set$Ensembl_Accession_ID,]

#calculamos los counts medios por muestra
cts_mean_nZ<-data.frame(TE_mean=colMeans(cts_TE_nZ),Hist_mean=colMeans(cts_hist_nZ),sexo=sample_info$sexo,genotipo=sample_info$genotipo,tejido=sample_info$Tissue,row.names = colnames(cts_TE_nZ))

#hacemos el scatter plot para ver si existe correlación entre los counts de histonas y los de transposones
correl_nZ<-ggplot(cts_mean_nZ,aes(x=TE_mean,y=Hist_mean))+
  geom_point(aes(color=genotipo))+
  geom_smooth(method = lm,)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "left", 
               label.y = "top") +
  labs(x="Lecturas normalizadas de transposones",y="Lecturas normalizadas de histonas")+
  facet_wrap(~tejido,scales = "free")

#lo repetimos pero separando por genotipo
correl_nZ_genot<-ggplot(cts_mean_nZ,aes(x=TE_mean,y=Hist_mean))+
  geom_point(aes(color=genotipo))+
  geom_smooth(method = lm,aes(group=genotipo))+
  stat_poly_eq(aes(group=genotipo,label = paste(..rr.label.., sep = "~~~")), label.x.npc = 'left', label.y.npc = 'top', size = 3, color = "black", parse = TRUE) +
  labs(x="Lecturas normalizadas de transposones",y="Lecturas normalizadas de histonas")+
  facet_wrap(~tejido)

#guardamos los plots
png(file=paste0(out_path,"Coreel_TE_Hist.png"),width = 10000, height = 8000, res = 1100)
correl
dev.off()

png(file=paste0(out_path,"CoreelnZ_TE_Hist.png"),width = 10000, height = 8000, res = 1100)
correl_nZ
dev.off()

png(file=paste0(out_path,"CoreelnZ_genot_TE_Hist.png"),width = 10000, height = 8000, res = 1100)
correl_nZ_genot
dev.off()

ggsave(paste0(out_path,"CoreelnZ_TE_Hist.pdf"),correl_nZ)
