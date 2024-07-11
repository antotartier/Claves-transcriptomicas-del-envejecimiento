
library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)


#función para analizar los resultados
speed_analysis<-function(path,sample_info){
res_fused<-read_tsv(path,id="sample")
res_fused$sample<-gsub("_Aligned.*$","",basename(res_fused$sample))
res_fused<-as.data.table(res_fused)

#separamos los datos donde vamos a trabajar con slope noZ
res_fused_noZ<-res_fused[,Change:=sum(slope_noZ<0), by=in_coords]
res_fused_noZ<-res_fused_noZ[Change==length(path)]
res_fused_noZ[,velocidad_noZ:=-1/slope_noZ, by=in_coords]
res_fused_noZ<-merge(res_fused_noZ,sample_info,by="sample")

#separamos los datos en los que vamos a trabjar con la pendiente normal
res_fused[,Change:=sum(slope<0), by=in_coords]
res_fused<-res_fused[Change==length(path)]
res_fused[,velocidad:=-1/slope, by=in_coords]
res_fused<-merge(res_fused,sample_info,by="sample")

# hacemos un subset de aquellos en los que la muestra joven con mayor velocidad tiene una velocidad menos que la senescente con menor velocidad
res_fused[,consist:=sum(max(velocidad[passage=="Temprano"])<min(velocidad[passage=="Tardío"])| max(velocidad[passage=="Tardío"])<min(velocidad[passage=="Temprano"])),by=in_coords]
res_consist<-res_fused[consist==1]

res_fused_noZ[,consist:=sum(max(velocidad_noZ[passage=="Temprano"])<min(velocidad_noZ[passage=="Tardío"])| max(velocidad_noZ[passage=="Tardío"])<min(velocidad_noZ[passage=="Temprano"])),by=in_coords]
res_consist_noZ<-res_fused_noZ[consist==1]

#calculamos el LFC de los intrones 
LFC_all<-res_fused[,.(LFC=log2(mean(velocidad[passage=="Tardío"])/mean(velocidad[passage=="Temprano"]))), by=.(in_coords)]
LFC_all[,Velocidad:=ifelse(LFC<0,"Disminuye","Aumenta")]

LFC_consist<-res_consist[,.(LFC=log2(mean(velocidad[passage=="Tardío"])/mean(velocidad[passage=="Temprano"]))), by=.(in_coords)]
LFC_consist[,Velocidad:=ifelse(LFC<0,"Disminuye","Aumenta")]


LFC_all_noZ<-res_fused_noZ[,.(LFC=log2(mean(velocidad_noZ[passage=="Tardío"])/mean(velocidad_noZ[passage=="Temprano"]))), by=.(in_coords)]
LFC_all_noZ[,velocidad_noZ:=ifelse(LFC<0,"Disminuye","Aumenta")]

LFC_consist_noZ<-res_consist_noZ[,.(LFC=log2(mean(velocidad_noZ[passage=="Tardío"])/mean(velocidad_noZ[passage=="Temprano"]))), by=.(in_coords)]
LFC_consist_noZ[,velocidad_noZ:=ifelse(LFC<0,"Disminuye","Aumenta")]

#calculamos la frecuencia de intrones en los que aumenta la velocidad
Freq_all<-LFC_all[,.(counts=.N), by=Velocidad]
Freq_all[,Frecuencia:=counts/nrow(LFC_all)]
Freq_all$group<-"Todos"

Freq_consist<-LFC_consist[,.(counts=.N), by=Velocidad]
Freq_consist[,Frecuencia:=counts/nrow(LFC_consist)]
Freq_consist$group<-"Consistentes"

Freq<-rbind(Freq_all,Freq_consist)

Freq_all_noZ<-LFC_all_noZ[,.(counts=.N), by=velocidad_noZ]
Freq_all_noZ[,Frecuencia:=counts/nrow(LFC_all_noZ)]
Freq_all_noZ$group<-"TodosnoZ"

Freq_consist_noZ<-LFC_consist_noZ[,.(counts=.N), by=velocidad_noZ]
Freq_consist_noZ[,Frecuencia:=counts/nrow(LFC_consist_noZ)]
Freq_consist_noZ$group<-"ConsistentesnoZ"

Freq_noZ<-rbind(Freq_all_noZ,Freq_consist_noZ)
colnames(Freq_noZ)<-c("Velocidad","counts","Frecuencia","group")

Freq_both<-rbind(Freq,Freq_noZ)

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
res<-list(Freq=Freq_noZ,LFC=Wilc_res_noZ)
return(res)
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------
#Figura TFM

sample_info_GM05565<-data.frame(sample=c("Tube_1", "Tube_2", "Tube_3", "Tube_4", "Tube_5", "Tube_6", "Tube_7", "Tube_8", "Tube_9"),passage=rep(c("Temprano","Medio","Tardío"),each=3))

sample_info_GM00038<-data.frame(sample=c("Tube_19", "Tube_20", "Tube_21", "Tube_22", "Tube_23", "Tube_24", "Tube_25" ,"Tube_26", "Tube_27"),passage=rep(c("Temprano","Medio","Tardío"),each=3))

# sample_info<-data.frame(sample=c("SRR5927738", "SRR5927739", "SRR5927740", "SRR5927741"),
#                         passage=rep(c("Temprano","Tardío"),each=2))

#path de los resultados
path_GM05565<-Sys.glob("/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/speed_res/allin1/*GM05565_speed.tsv")
path_GM00038<-Sys.glob("/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/speed_res/allin1/*GM00038_speed.tsv")

res_GM05565<-speed_analysis(path_GM05565,sample_info_GM05565)
res_GM00038<-speed_analysis(path_GM00038,sample_info_GM00038)

GM05565_Freq<-res_GM05565$Freq
GM00038_Freq<-res_GM00038$Freq

GM05565_Freq$group<-gsub("noZ","",GM05565_Freq$group)
GM05565_Freq$group<- paste0("GM05565","(",GM05565_Freq$group,")")
GM00038_Freq$group<-gsub("noZ","",GM00038_Freq$group)
GM00038_Freq$group<- paste0("GM00038","(",GM00038_Freq$group,")")

Freq_TFM<-rbind(GM05565_Freq,GM00038_Freq)

GM05565_LFC<-res_GM05565$LFC
GM00038_LFC<-res_GM00038$LFC

GM05565_LFC$línea<-"GM05565"
GM00038_LFC$línea<-"GM00038"

LFC_TFM<-rbind(GM05565_LFC,GM00038_LFC)

#representamos los resultados
Freq_plot<-ggplot(Freq_TFM, aes(fill=Velocidad, y=Frecuencia, x=group)) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Fracción de intrones")+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 12, angle = 90, hjust = 1),axis.title.y = element_text(size = 15),legend.text = element_text(size=12),axis.text.y = element_text(size = 10),legend.title = element_text(size = 12))+
  scale_fill_discrete(name="Cambio en la velocidad")
  
LFC_plot<-ggplot(LFC_TFM,aes(x=línea, y=media, fill=group))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.9))+
  theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14,face = "bold"),legend.title = element_blank(),axis.title.y = element_text(size = 18),legend.text = element_text(size=15),axis.text.y = element_text(size = 12))+
  labs(y=expression("Velocidad de la elongación (log"[2]*"FC)"))+
  geom_errorbar(aes(x=línea, ymin=minus, ymax=plus), width=0.2, colour="black", alpha=0.9, size=0.5,position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("Intrones consistentes"="#FFB347", "Todos los intrones"="#AEC6CF"))


png(filename = "/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/speed_res/allin1/Freq_both.png",height =8000 , width =8000 ,res =1100 )
Freq_plot
dev.off()

png(filename = "/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/speed_res/allin1/LFC_both.png",height =6000 , width =8000 ,res =1100 )
LFC_plot
dev.off()

ggsave("/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/speed_res/allin1/Freq_both.pdf",Freq_plot)
ggsave("/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/speed_res/allin1/LFC_both.pdf",LFC_plot)
