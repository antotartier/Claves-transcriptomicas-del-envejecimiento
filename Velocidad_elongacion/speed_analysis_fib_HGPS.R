library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)

#cargamos el sample info
sample_info<-fread("/home/antotartier/data/fib_HGPS/FASTQ/Sra_meta.tsv")
names(sample_info)[names(sample_info)=="SRA"]<-"sample"

#unimos los resultados 
path_res<-Sys.glob("/home/antotartier/data/velocidad_transcrip/fib_HGPS/speed_res/*_speed.tsv")
res_fused<-read_tsv(path_res,id="sample",show_col_types = F)
res_fused$sample<-gsub("_Aligned.*$","",basename(res_fused$sample))
res_fused<-as.data.table(res_fused)

#separamos los datos donde vamos a trabajar con slope noZ
res_fused_noZ<-res_fused[,Change:=sum(slope_noZ<0), by=in_coords]
res_fused_noZ<-res_fused_noZ[Change==length(path_res)]
res_fused_noZ[,velocidad_noZ:=-1/slope_noZ, by=in_coords]
res_fused_noZ<-merge(res_fused_noZ,sample_info,by="sample")

#separamos los datos en los que vamos a trabjar con la pendiente normal
res_fused[,Change:=sum(slope<0), by=in_coords]
res_fused<-res_fused[Change==length(path_res)]
res_fused[,velocidad:=-1/slope, by=in_coords]
res_fused<-merge(res_fused,sample_info,by="sample")




# hacemos un subset de aquellos en los que la muestra joven con mayor velocidad tiene una velocidad menos que la senescente con menor velocidad
res_fused[,consist:=sum(max(velocidad[Genotipo=="WT"])<min(velocidad[Genotipo=="HGPS"])| max(velocidad[Genotipo=="HGPS"])<min(velocidad[Genotipo=="WT"])),by=in_coords]
res_consist<-res_fused[consist==1]

res_fused_noZ[,consist:=sum(max(velocidad_noZ[Genotipo=="WT"])<min(velocidad_noZ[Genotipo=="HGPS"])| max(velocidad_noZ[Genotipo=="HGPS"])<min(velocidad_noZ[Genotipo=="WT"])),by=in_coords]
res_consist_noZ<-res_fused_noZ[consist==1]

#calculamos el LFC de los intrones 
LFC_all<-res_fused[,.(LFC=log2(mean(velocidad[Genotipo=="HGPS"])/mean(velocidad[Genotipo=="WT"]))), by=.(in_coords)]
LFC_all[,Velocidad:=ifelse(LFC<0,"Disminuye","Aumenta")]

LFC_consist<-res_consist[,.(LFC=log2(mean(velocidad[Genotipo=="HGPS"])/mean(velocidad[Genotipo=="WT"]))), by=.(in_coords)]
LFC_consist[,Velocidad:=ifelse(LFC<0,"Disminuye","Aumenta")]


LFC_all_noZ<-res_fused_noZ[,.(LFC=log2(mean(velocidad_noZ[Genotipo=="HGPS"])/mean(velocidad_noZ[Genotipo=="WT"]))), by=.(in_coords)]
LFC_all_noZ[,velocidad_noZ:=ifelse(LFC<0,"Disminuye","Aumenta")]

LFC_consist_noZ<-res_consist_noZ[,.(LFC=log2(mean(velocidad_noZ[Genotipo=="HGPS"])/mean(velocidad_noZ[Genotipo=="WT"]))), by=.(in_coords)]
LFC_consist_noZ[,velocidad_noZ:=ifelse(LFC<0,"Disminuye","Aumenta")]

#calculamos la frecuencia de intrones en los que aumenta la velocidad
Freq_all<-LFC_all[,.(counts=.N), by=Velocidad]
Freq_all[,Frecuencia:=counts/nrow(LFC_all)]
Freq_all$group<-"Todos los intrones"

Freq_consist<-LFC_consist[,.(counts=.N), by=Velocidad]
Freq_consist[,Frecuencia:=counts/nrow(LFC_consist)]
Freq_consist$group<-"Intrones consistentes"

Freq<-rbind(Freq_all,Freq_consist)

Freq_all_noZ<-LFC_all_noZ[,.(counts=.N), by=velocidad_noZ]
Freq_all_noZ[,Frecuencia:=counts/nrow(LFC_all_noZ)]
Freq_all_noZ$group<-"Todos los intrones"

Freq_consist_noZ<-LFC_consist_noZ[,.(counts=.N), by=velocidad_noZ]
Freq_consist_noZ[,Frecuencia:=counts/nrow(LFC_consist_noZ)]
Freq_consist_noZ$group<-"Intrones consistentes"

Freq_noZ<-rbind(Freq_all_noZ,Freq_consist_noZ)
colnames(Freq_noZ)<-c("Velocidad","counts","Frecuencia","group")

Freq_both<-rbind(Freq,Freq_noZ)

#representamos los resultados
Freq_plot<-ggplot(Freq_noZ, aes(fill=Velocidad, y=Frecuencia, x=group)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 12))

#hacemos un test de wilcoxon para calcular las diferencias de velocidad entre las dos condiciones
Wilc_all<-wilcox.test(LFC_all$LFC, conf.int=T)
Wilc_all_res<-data.frame(media=Wilc_all$estimate,row.names = NULL)
Wilc_all_res$plus<-Wilc_all$conf.int[[1]]
Wilc_all_res$minus<-Wilc_all$conf.int[[2]]
Wilc_all_res$group<-"Todos los intrones"
Wilc_all_res$metodo<-"Normal"

Wilc_consist<-wilcox.test(LFC_consist$LFC, conf.int=T)
Wilc_consist_res<-data.frame(media=Wilc_consist$estimate,row.names = NULL)
Wilc_consist_res$plus<-Wilc_consist$conf.int[[1]]
Wilc_consist_res$minus<-Wilc_consist$conf.int[[2]]
Wilc_consist_res$group<-"Intrones consistentes"
Wilc_consist_res$metodo<-"Normal"

Wilc_res<-rbind(Wilc_all_res,Wilc_consist_res)

Wilc_all_noZ<-wilcox.test(LFC_all_noZ$LFC, conf.int=T)
Wilc_all_res_noZ<-data.frame(media=Wilc_all_noZ$estimate,row.names = NULL)
Wilc_all_res_noZ$plus<-Wilc_all_noZ$conf.int[[1]]
Wilc_all_res_noZ$minus<-Wilc_all_noZ$conf.int[[2]]
Wilc_all_res_noZ$group<-"Todos los intrones"
Wilc_all_res_noZ$metodo<-"NoZ"

Wilc_consist_noZ<-wilcox.test(LFC_consist_noZ$LFC, conf.int=T)
Wilc_consist_res_noZ<-data.frame(media=Wilc_consist_noZ$estimate,row.names = NULL)
Wilc_consist_res_noZ$plus<-Wilc_consist_noZ$conf.int[[1]]
Wilc_consist_res_noZ$minus<-Wilc_consist_noZ$conf.int[[2]]
Wilc_consist_res_noZ$group<-"Intrones consistentes"
Wilc_consist_res_noZ$metodo<-"NoZ"

Wilc_res_noZ<-rbind(Wilc_all_res_noZ,Wilc_consist_res_noZ)

Wilc_both<-rbind(Wilc_res,Wilc_res_noZ)

#representamos los resultados
LFC_plot<-ggplot(Wilc_res_noZ,aes(x=group, y=media, fill=group))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.9))+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14))+
  labs(y=expression("Log"[2]*"FC"))+
  geom_errorbar(aes(x=group, ymin=minus, ymax=plus), width=0.2, colour="black", alpha=0.9, size=0.5,position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("Intrones consistentes"="#FFB347", "Todos los intrones"="#AEC6CF"))
  



#Guardamos los resultados
png(filename = "/home/antotartier/data/velocidad_transcrip/fib_HGPS/speed_res/plots_fib_HGPS.png",height =6000 , width =18000 ,res =1100 )
grid.arrange(LFC_plot,Freq_plot,ncol=2)
dev.off()

ggsave("/home/antotartier/data/velocidad_transcrip/fib_HGPS/speed_res/Freq_fib_HGPS.pdf",Freq_plot)
ggsave("/home/antotartier/data/velocidad_transcrip/fib_HGPS/speed_res/LFC_fib_HGPS.pdf",LFC_plot)