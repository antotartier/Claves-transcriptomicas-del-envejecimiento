#cargamos las librerías necesarias
library(ggplot2)
library(ggpmisc)

#cargamos los datos del cambio en la velocidad y el LFC medio de las histonas en los diferentes tejidos (generado al final del script "histone_LAKI.R")
load("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_all.obj")
load("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_machos.obj")
load("/home/antotartier/data/velocidad_transcrip/LAKI/speed_res/LFC_hembras.obj")
load("/home/antotartier/data/histones/LAKI/LFC_hist.obj")
load("/home/antotartier/data/histones/LAKI/LFC_hist_M.obj")
load("/home/antotartier/data/histones/LAKI/LFC_hist_F.obj")

#nos quedamos con los LFC de la velocidad de todos los intrones porque son los que tienen menos error y con el metodo noZ
LFC_tissue<-LFC_tissue[LFC_tissue$grupo=="Todos los intrones"&LFC_tissue$metodo=="NoZ",c(1,6)]
LFC_tissue_M<-LFC_tissue_M[LFC_tissue_M$grupo=="Todos los intrones"&LFC_tissue_M$metodo=="NoZ",c(1,6)]
LFC_tissue_H<-LFC_tissue_H[LFC_tissue_H$grupo=="Todos los intrones"&LFC_tissue_H$metodo=="NoZ",c(1,6)]

#juntamos los datos de histonas y velocidad
correl<-merge(LFC_tissue,LFC_hist,by="tejido")
correl_M<-merge(LFC_tissue_M,LFC_hist_M,by="tejido")
correl_H<-merge(LFC_tissue_H,LFC_hist_H,by="tejido")

#representamos los datos
plot_correl<-ggplot(correl,aes(x=LFC,y=media))+
  geom_point(aes(color=tejido))+
  geom_smooth(method = lm,)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "left", 
               label.y = "top") +
  labs(x="LFC medio de las histonas",y="Velocidad de la elongación media")+
  ggtitle("Correlación expresión de las histonas y 
  velocidad de la elongación LAKI")+
  theme(plot.title = element_text(hjust = 0.5))

plot_correl_M<-ggplot(correl_M,aes(x=LFC,y=media))+
  geom_point(aes(color=tejido))+
  geom_smooth(method = lm,)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "left", 
               label.y = "top") +
  labs(x="LFC medio de las histonas",y="Velocidad de la elongación media")+
  ggtitle("Correlación expresión de las histonas y 
  velocidad de la elongación LAKI machos")+
  theme(plot.title = element_text(hjust = 0.5))

plot_correl_H<-ggplot(correl_H,aes(x=LFC,y=media))+
  geom_point(aes(color=tejido))+
  geom_smooth(method = lm,)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, 
               size = 4,
               label.x = "left", 
               label.y = "top") +
  labs(x="LFC medio de las histonas",y="Velocidad de la elongación media")+
  ggtitle("Correlación expresión de las histonas y 
  velocidad de la elongación LAKI hembras")+
  theme(plot.title = element_text(hjust = 0.5))

#guardamos el resultado
png(file="/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_Hist.png",width = 20000, height = 5000, res = 1100)
grid.arrange(plot_correl,plot_correl_M,plot_correl_H,ncol=3)
dev.off()

ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_Hist.pdf",plot_correl)
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_Hist_M.pdf",plot_correl_M)
ggsave("/home/antotartier/data/velocidad_transcrip/LAKI/Correl_speed_Hist_H.pdf",plot_correl_H)
