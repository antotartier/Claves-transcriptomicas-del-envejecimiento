#Loading functions
library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringi)

#cargamos el sample info y eliminamos las muestras ausentes, las de FACE y el outlier
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")
elim<-c("Tube_39","Tube_56","Tube_57","Tube_66","Tube_86","Tube_87","Tube_115","Tube_116","Tube_117","Tube_118","Tube_119","Tube_120")
sample_info<-sample_info[!grepl(paste0(elim,collapse = "|"),SampleName),]

#seaparmos las variables en columnas
sample_info<-sample_info %>% mutate(genotipo=str_extract(SampleID,"KO|WT"),
                                    sexo=ifelse(grepl("M",SampleID),"M","H"))

#función que extrae los resultados para cada tejido
res_speed<-function(path,sample_info,sex="ALL"){
#lista para guardar los dos resultados
results<-list()
#unimos los resultados 
path_res<-Sys.glob(path)
tis<-gsub("\\*","",gsub("_speed.tsv","",basename(path)))
res_fused<-read_tsv(path_res,id="SampleName",show_col_types = F)
res_fused$SampleName<-gsub(paste0("_",tis,".*$"),"",basename(res_fused$SampleName))
filt<-sample_info[sample_info$Tissue==tis,]
if(sex=="M"){
res_fused<-res_fused[res_fused$SampleName %in% filt$SampleName[filt$sexo=="M"],]
} else if (sex=="H"){
res_fused<-res_fused[res_fused$SampleName %in% filt$SampleName[filt$sexo=="H"],]
}
res_fused<-as.data.table(res_fused)

#separamos los datos donde vamos a trabajar con slope noZ
res_fused_noZ<-res_fused[,Change:=sum(slope_noZ<0), by=in_coords]
res_fused_noZ<-res_fused_noZ[Change==length(unique(res_fused$SampleName))]
res_fused_noZ[,velocidad_noZ:=-1/slope_noZ, by=in_coords]
res_fused_noZ<-merge(res_fused_noZ,sample_info,by="SampleName")

#separamos los datos en los que vamos a trabjar con la pendiente normal
res_fused[,Change:=sum(slope<0), by=in_coords]
res_fused<-res_fused[Change==length(unique(res_fused$SampleName))]
res_fused[,velocidad:=-1/slope, by=in_coords]
res_fused<-merge(res_fused,sample_info,by="SampleName")


# hacemos un subset de aquellos en los que la muestra joven con mayor velocidad tiene una velocidad menos que la senescente con menor velocidad
res_fused[,consist:=sum(max(velocidad[genotipo=="WT"])<min(velocidad[genotipo=="KO"])| max(velocidad[genotipo=="KO"])<min(velocidad[genotipo=="WT"])),by=in_coords]
res_consist<-res_fused[consist==1]

res_fused_noZ[,consist:=sum(max(velocidad_noZ[genotipo=="WT"])<min(velocidad_noZ[genotipo=="KO"])| max(velocidad_noZ[genotipo=="KO"])<min(velocidad_noZ[genotipo=="WT"])),by=in_coords]
res_consist_noZ<-res_fused_noZ[consist==1]

#aqui puede filtrarse por sexo
# res_fused<-res_fused[sexo=="M"]
# res_consist_noZ<-res_consist_noZ[sexo=="M"]

#calculamos el LFC de los intrones 
LFC_all<-res_fused[,.(LFC=log2(mean(velocidad[genotipo=="KO"])/mean(velocidad[genotipo=="WT"]))), by=.(in_coords)]
LFC_all[,Velocidad:=ifelse(LFC<0,"Disminuye","Aumenta")]

LFC_consist<-res_consist[,.(LFC=log2(mean(velocidad[genotipo=="KO"])/mean(velocidad[genotipo=="WT"]))), by=.(in_coords)]
LFC_consist[,Velocidad:=ifelse(LFC<0,"Disminuye","Aumenta")]


LFC_all_noZ<-res_fused_noZ[,.(LFC=log2(mean(velocidad_noZ[genotipo=="KO"])/mean(velocidad_noZ[genotipo=="WT"]))), by=.(in_coords)]
LFC_all_noZ[,velocidad_noZ:=ifelse(LFC<0,"Disminuye","Aumenta")]

LFC_consist_noZ<-res_consist_noZ[,.(LFC=log2(mean(velocidad_noZ[genotipo=="KO"])/mean(velocidad_noZ[genotipo=="WT"]))), by=.(in_coords)]
LFC_consist_noZ[,velocidad_noZ:=ifelse(LFC<0,"Disminuye","Aumenta")]

#calculamos la frecuencia de intrones en los que aumenta la velocidad
Freq_all<-LFC_all[,.(counts=.N), by=Velocidad]
Freq_all[,Frecuencia:=counts/nrow(LFC_all)]
Freq_all$grupo<-"Todos"
Freq_all$tejido<-tis
Freq_all$metodo<-"Normal"

Freq_consist<-LFC_consist[,.(counts=.N), by=Velocidad]
Freq_consist[,Frecuencia:=counts/nrow(LFC_consist)]
Freq_consist$grupo<-"Consistentes"
Freq_consist$tejido<-tis
Freq_consist$metodo<-"Normal"

Freq<-rbind(Freq_all,Freq_consist)

Freq_all_noZ<-LFC_all_noZ[,.(counts=.N), by=velocidad_noZ]
Freq_all_noZ[,Frecuencia:=counts/nrow(LFC_all_noZ)]
Freq_all_noZ$grupo<-"Todos"
Freq_all_noZ$tejido<-tis
Freq_all_noZ$metodo<-"NoZ"

Freq_consist_noZ<-LFC_consist_noZ[,.(counts=.N), by=velocidad_noZ]
Freq_consist_noZ[,Frecuencia:=counts/nrow(LFC_consist_noZ)]
Freq_consist_noZ$grupo<-"Consistentes"
Freq_consist_noZ$tejido<-tis
Freq_consist_noZ$metodo<-"NoZ"

Freq_noZ<-rbind(Freq_all_noZ,Freq_consist_noZ)
colnames(Freq_noZ)[colnames(Freq_noZ)=="velocidad_noZ"]<-"Velocidad"#cambiamos el nombre a la columna para evitar problemas al juntarlo con los resultados normales

Freq_both<-rbind(Freq,Freq_noZ)

results[["Frec"]]<-Freq_both

#hacemos un test de wilcoxon para calcular las diferencias de velocidad entre las dos condiciones
Wilc_all<-wilcox.test(LFC_all$LFC, conf.int=T)
Wilc_all_res<-data.frame(media=Wilc_all$estimate,row.names = NULL)
Wilc_all_res$plus<-Wilc_all$conf.int[[1]]
Wilc_all_res$minus<-Wilc_all$conf.int[[2]]
Wilc_all_res$grupo<-"Todos los intrones"
Wilc_all_res$metodo<-"Normal"
Wilc_all_res$tejido<-tis

Wilc_consist<-wilcox.test(LFC_consist$LFC, conf.int=T)
Wilc_consist_res<-data.frame(media=Wilc_consist$estimate,row.names = NULL)
Wilc_consist_res$plus<-Wilc_consist$conf.int[[1]]
Wilc_consist_res$minus<-Wilc_consist$conf.int[[2]]
Wilc_consist_res$grupo<-"Intrones consistentes"
Wilc_consist_res$metodo<-"Normal"
Wilc_consist_res$tejido<-tis

Wilc_res<-rbind(Wilc_all_res,Wilc_consist_res)

Wilc_all_noZ<-wilcox.test(LFC_all_noZ$LFC, conf.int=T)
Wilc_all_res_noZ<-data.frame(media=Wilc_all_noZ$estimate,row.names = NULL)
Wilc_all_res_noZ$plus<-Wilc_all_noZ$conf.int[[1]]
Wilc_all_res_noZ$minus<-Wilc_all_noZ$conf.int[[2]]
Wilc_all_res_noZ$grupo<-"Todos los intrones"
Wilc_all_res_noZ$metodo<-"NoZ"
Wilc_all_res_noZ$tejido<-tis

Wilc_consist_noZ<-wilcox.test(LFC_consist_noZ$LFC, conf.int=T)
Wilc_consist_res_noZ<-data.frame(media=Wilc_consist_noZ$estimate,row.names = NULL)
Wilc_consist_res_noZ$plus<-Wilc_consist_noZ$conf.int[[1]]
Wilc_consist_res_noZ$minus<-Wilc_consist_noZ$conf.int[[2]]
Wilc_consist_res_noZ$grupo<-"Intrones consistentes"
Wilc_consist_res_noZ$metodo<-"NoZ"
Wilc_consist_res_noZ$tejido<-tis

Wilc_res_noZ<-rbind(Wilc_all_res_noZ,Wilc_consist_res_noZ)

Wilc_both<-rbind(Wilc_res,Wilc_res_noZ)

results[["LFC"]]<-Wilc_both

return(results)
}

#ejecutamos la función para cada tejido
tissue_path<-list(Colon="/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/*Colon_speed.tsv",Íleon="/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/*Ileon_speed.tsv",Riñón="/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/*Kidney_speed.tsv",Hígado="/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/*Liver_speed.tsv",Corazón="/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/*Heart_speed.tsv",Músculo="/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/*Muscle_speed.tsv")

res_alltis<-lapply(tissue_path,res_speed,sample_info)

#juntamos los resultados de los diferentes tejidos
Freq_tissue<-sapply(res_alltis,"[[",1,simplify = F)
Freq_tissue<-do.call(rbind,Freq_tissue)
Freq_tissue$tejido<-stri_replace_all_regex(Freq_tissue$tejido,pattern = c("Ileon","Kidney","Liver","Heart" ,"Muscle"),replacement = c("Íleon","Riñón","Hígado","Corazón","Músculo"),vectorize_all = F)

LFC_tissue<-sapply(res_alltis,"[[",2,simplify = F)
LFC_tissue<-do.call(rbind,LFC_tissue)
LFC_tissue$tejido<-stri_replace_all_regex(LFC_tissue$tejido,pattern = c("Ileon","Kidney","Liver","Heart" ,"Muscle"),replacement = c("Íleon","Riñón","Hígado","Corazón","Músculo"),vectorize_all = F)

#nos quedamos con el método NoZ
Freq_tissue_NoZ<-Freq_tissue[Freq_tissue$metodo=="NoZ",]
LFC_tissue_Noz<-LFC_tissue[LFC_tissue$metodo=="NoZ",]

#representamos los resultados
Freq_plot<-ggplot(Freq_tissue_NoZ, aes(fill=Velocidad, y=Frecuencia, x=interaction(tejido,grupo,lex.order = T))) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Fracción de intrones")+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 12, angle = 90, hjust = 1))
  
LFC_plot<-ggplot(LFC_tissue_Noz,aes(x=tejido, y=media, fill=grupo))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.9))+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14,angle = 90))+
  labs(y=expression("Velocidad de la elongación (Log"[2]*"FC)"))+
  geom_errorbar(aes(x=tejido, ymin=minus, ymax=plus), width=0.2, colour="black", alpha=0.9, size=0.5,position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("Intrones consistentes"="#FFB347", "Todos los intrones"="#AEC6CF"))
 


#Guardamos los resultados
png(filename = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/Freq.png",height =15000 , width =17000 ,res =1100 )
Freq_plot
dev.off()
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/Freq.pdf",Freq_plot)

png(filename = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC.png",height =12000 , width =15000 ,res =1100 )
LFC_plot
dev.off()
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC.pdf",LFC_plot)

#----------------------------------------------------------------------------------------------------------------------------------------
#análisis por separado de machos y hembras

#machos
res_alltis_M<-lapply(tissue_path,res_speed,sample_info,sex="M")

#juntamos los resultados de los diferentes tejidos
Freq_tissue_M<-sapply(res_alltis_M,"[[",1,simplify = F)
Freq_tissue_M<-do.call(rbind,Freq_tissue_M)
Freq_tissue_M$tejido<-stri_replace_all_regex(Freq_tissue_M$tejido,pattern = c("Ileon","Kidney","Liver","Heart" ,"Muscle"),replacement = c("Íleon","Riñón","Hígado","Corazón","Músculo"),vectorize_all = F)

LFC_tissue_M<-sapply(res_alltis_M,"[[",2,simplify = F)
LFC_tissue_M<-do.call(rbind,LFC_tissue_M)
LFC_tissue_M$tejido<-stri_replace_all_regex(LFC_tissue_M$tejido,pattern = c("Ileon","Kidney","Liver","Heart" ,"Muscle"),replacement = c("Íleon","Riñón","Hígado","Corazón","Músculo"),vectorize_all = F)

#nos quedamos solo con el método noZ
Freq_tissue_NoZ_M<-Freq_tissue_M[Freq_tissue_M$metodo=="NoZ",]
LFC_tissue_Noz_M<-LFC_tissue_M[LFC_tissue_M$metodo=="NoZ",]

#representamos los resultados
Freq_plot_M<-ggplot(Freq_tissue_NoZ_M, aes(fill=Velocidad, y=Frecuencia, x=interaction(tejido,grupo,lex.order = T))) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Fracción de intrones")+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 12, angle = 90, hjust = 1))

LFC_plot_M<-ggplot(LFC_tissue_Noz_M,aes(x=tejido, y=media, fill=grupo))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.9))+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14,angle = 90))+
  labs(y=expression("Velocidad de la elongación (Log"[2]*"FC)"))+
  geom_errorbar(aes(x=tejido, ymin=minus, ymax=plus), width=0.2, colour="black", alpha=0.9, size=0.5,position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("Intrones consistentes"="#FFB347", "Todos los intrones"="#AEC6CF"))

#Guardamos los resultados
png(filename = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/Freq_machos.png",height =15000 , width =17000 ,res =1100 )
Freq_plot_M
dev.off()

png(filename = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_machos.png",height =12000 , width =15000 ,res =1100 )
LFC_plot_M
dev.off()

ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/Freq_machos.pdf",Freq_plot_M)
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_machos.pdf",LFC_plot_M)


#hembras
res_alltis_H<-lapply(tissue_path,res_speed,sample_info,sex="H")

#juntamos los resultados de los diferentes tejidos
Freq_tissue_H<-sapply(res_alltis_H,"[[",1,simplify = F)
Freq_tissue_H<-do.call(rbind,Freq_tissue_H)
Freq_tissue_H$tejido<-stri_replace_all_regex(Freq_tissue_H$tejido,pattern = c("Ileon","Kidney","Liver","Heart" ,"Muscle"),replacement = c("Íleon","Riñón","Hígado","Corazón","Músculo"),vectorize_all = F)

LFC_tissue_H<-sapply(res_alltis_H,"[[",2,simplify = F)
LFC_tissue_H<-do.call(rbind,LFC_tissue_H)
LFC_tissue_H$tejido<-stri_replace_all_regex(LFC_tissue_H$tejido,pattern = c("Ileon","Kidney","Liver","Heart" ,"Muscle"),replacement = c("Íleon","Riñón","Hígado","Corazón","Músculo"),vectorize_all = F)


#nos quedamos solo con el método noZ
Freq_tissue_NoZ_H<-Freq_tissue_H[Freq_tissue_H$metodo=="NoZ",]
LFC_tissue_Noz_H<-LFC_tissue_H[LFC_tissue_H$metodo=="NoZ",]

#representamos los resultados
Freq_plot_H<-ggplot(Freq_tissue_NoZ_H, aes(fill=Velocidad, y=Frecuencia, x=interaction(tejido,grupo,lex.order = T))) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Fracción de intrones")+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size = 12, angle = 90, hjust = 1))

LFC_plot_H<-ggplot(LFC_tissue_Noz_H,aes(x=tejido, y=media, fill=grupo))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.9))+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14,angle = 90))+
  labs(y=expression("Velocidad de la elongación (Log"[2]*"FC)"))+
  geom_errorbar(aes(x=tejido, ymin=minus, ymax=plus), width=0.2, colour="black", alpha=0.9, size=0.5,position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("Intrones consistentes"="#FFB347", "Todos los intrones"="#AEC6CF"))

#Guardamos los resultados
png(filename = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/Freq_hembras.png",height =15000 , width =17000 ,res =1100 )
Freq_plot_H
dev.off()

png(filename = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_hembras.png",height =12000 , width =15000 ,res =1100 )
LFC_plot_H
dev.off()

ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/Freq_hembras.pdf",Freq_plot_H)
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_hembras.pdf",LFC_plot_H)


#----------------------------------------------------------------------------------------------------------
#guardamos los resultados
save(LFC_tissue,file = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_all.obj")
save(LFC_tissue_H,file = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_hembras.obj")
save(LFC_tissue_M,file = "/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_machos.obj")
