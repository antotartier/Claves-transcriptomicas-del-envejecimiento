library(ggplot2)
library(ggpmisc)


#cargamos los resultados conjuntos, los de machos y los de hemrbas; y los de la velocidad (generados en el script "speed_analysis_LAKI.R")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_all.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_M.obj")
load("/home/antotartier/data/transposones/LAKI/resultados_anlaisis/res_F.obj")

load("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_all.obj")
load("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_machos.obj")
load("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_hembras.obj")

#nos quedamos con los LFC de la velocidad de todos los intrones porque son los que tienen menos error y con el metodo noZ
LFC_tissue<-LFC_tissue[LFC_tissue$grupo=="Todos los intrones"&LFC_tissue$metodo=="NoZ",c(1,6)]
LFC_tissue_M<-LFC_tissue_M[LFC_tissue_M$grupo=="Todos los intrones"&LFC_tissue_M$metodo=="NoZ",c(1,6)]
LFC_tissue_H<-LFC_tissue_H[LFC_tissue_H$grupo=="Todos los intrones"&LFC_tissue_H$metodo=="NoZ",c(1,6)]

#nos quedamos con los transposones
get_TE<-function(x){
  TE<-x[grepl(":",rownames(x)),]
  return(TE)
}
TE_all<-lapply(res_all,get_TE)
TE_M<-lapply(res_M,get_TE)
TE_F<-lapply(res_F,get_TE)

#calculamos la media del LFC
tis<-c("Íleon","Riñón","Músculo","Hígado","Colon","Corazón")
LFC_mean<-function(x,tis){
y<-mean(x$log2FoldChange)
LFC<-data.frame(LFC_TE=y,tejido= tis)
}

LFC_TE_all<-mapply(LFC_mean,TE_all,tis,SIMPLIFY = F)
LFC_TE_all<-do.call(rbind,LFC_TE_all)

LFC_TE_M<-mapply(LFC_mean,TE_M,tis,SIMPLIFY = F)
LFC_TE_M<-do.call(rbind,LFC_TE_M)

LFC_TE_F<-mapply(LFC_mean,TE_F,tis,SIMPLIFY = F)
LFC_TE_F<-do.call(rbind,LFC_TE_F)

#juntamos los datos de histonas y velocidad
correl<-merge(LFC_tissue,LFC_TE_all,by="tejido")
correl_M<-merge(LFC_tissue_M,LFC_TE_M,by="tejido")
correl_F<-merge(LFC_tissue_H,LFC_TE_F,by="tejido")

#representamos los datos

#todos
plot_correl<-ggplot(correl, aes(x = LFC_TE, y = media)) +
  geom_point(aes(color = tejido)) +
  geom_smooth(method = "lm") +
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "right", 
               label.y = "top") + 
  labs(x = "LFC medio de los transposones", y = "Velocidad de la elongación media") +
  ggtitle("Correlación expresión de los transposones y velocidad de la elongación LAKI") +
  theme(plot.title = element_text(hjust = 0.5))

#machos
plot_correl_M<-ggplot(correl_M,aes(x=LFC_TE,y=media))+
  geom_point(aes(color=tejido))+
  geom_smooth(method = lm,)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "right", 
               label.y = "top") +
  labs(x="LFC medio de los transposones",y="Velocidad de la elongación media")+
  ggtitle("Correlación expresión de los transposones y 
  velocidad de la elongación LAKI machos")+
  theme(plot.title = element_text(hjust = 0.5))

#hembras
plot_correl_F<-ggplot(correl_F,aes(x=LFC_TE,y=media))+
  geom_point(aes(color=tejido))+
  geom_smooth(method = lm,)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "right", 
               label.y = "top") +
  labs(x="LFC medio de los transposones",y="Velocidad de la elongación media")+
  ggtitle("Correlación expresión de  transposones y 
  velocidad de la elongación LAKI hembras")+
  theme(plot.title = element_text(hjust = 0.5))

#guardamos el resultado
png(file="/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_TE.png",width = 20000, height = 5000, res = 1100)
grid.arrange(plot_correl,plot_correl_M,plot_correl_F,ncol=3)
dev.off()

ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_TE.pdf",plot_correl)
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_TE_M.pdf",plot_correl_M)
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_TE_F.pdf",plot_correl_F)
