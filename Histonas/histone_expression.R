library(ggplot2)
library(ggsignif)

#cargamos el resultado de DESeq y los id y los nombres de las histonas, los genes y los pseudogenes
load("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/DESeq_res_conj.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_genes.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_pseudogenes.obj")


#lo convertimos en un dataframe
res_df<-as.data.frame(res)
rownames(res_df)<-gsub("\\.[0-9]+$","",rownames(res_df))#eliminamos el punto de los ensemble IDs

#nos quedamos con las histonas
res_hist<-res_df[rownames(res_df)%in%histone_name_id$ensembl_gene_id,]
#sustituimos los id por los nombres
id_to_name<-setNames(histone_name_id$symbol,histone_name_id$ensembl_gene_id)
rownames(res_hist)<-id_to_name[rownames(res_hist)]

#hacemos una lista de los nombres de las histonas cluster y no cluster
cluster<-histone_name_id$symbol[grepl("cluster",histone_name_id$symbol)]
cluster<-cluster[cluster %in% histone_genes$symbol]
no_cluster<-histone_name_id$symbol[!grepl("cluster",histone_name_id$symbol)]
no_cluster<-no_cluster[!(no_cluster %in% histone_pseudogenes$symbol)]

#hacemos una lista de los nombres de las histonas por familia
H1<-histone_genes$symbol[grepl("H1",histone_genes$symbol)]
H2A<-histone_genes$symbol[grepl("H2A",histone_genes$symbol)]
H2B<-histone_genes$symbol[grepl("H2B",histone_genes$symbol)]
H3<-histone_genes$symbol[grepl("H3",histone_genes$symbol)]
H4<-histone_genes$symbol[grepl("H4",histone_genes$symbol)]

#extraemos los genes up-regulados
up_reg<-res_hist[res_hist$log2FoldChange>0 & res_hist$padj<0.05,]
#definimos un patrón para seleccionar los genes dependientes de la replicación
dpts<-c("cluster","H2AX","H4C16")


#clasificamos los genes
res_hist<-res_hist %>%
  mutate(type=ifelse(rownames(res_hist)%in%histone_pseudogenes$symbol, "pseudogen",
         ifelse(rownames(res_hist)%in%histone_genes$symbol,"gen", "unknown")),
         clustering=ifelse(rownames(res_hist)%in%cluster, "cluster",
          ifelse(rownames(res_hist)%in%no_cluster,"No cluster", "unknown")),
         familia=ifelse(rownames(res_hist)%in%H1, "H1",
                ifelse(rownames(res_hist)%in%H2A,"H2A",
                ifelse(rownames(res_hist)%in%H2B,"H2B",
                ifelse(rownames(res_hist)%in%H3,"H3",
                ifelse(rownames(res_hist)%in%H4,"H4", "Unknown"))))),
         reg=ifelse(grepl(paste(dpts, collapse = "|"),rownames(res_hist)) & !grepl("H3-4",rownames(res_hist)),"HDR",
            ifelse(!grepl(paste(dpts, collapse = "|"),rownames(res_hist)) | grepl("H3-4",rownames(res_hist)),"HIR"
                   ,"unknown")))
#nos quedamos con los resultados estadísticamente significativos
res_hist_sig<-res_hist[res_hist$padj<0.05,]
#eliminamos los pseudogenes
res_hist_sig_nopseudo<-res_hist_sig[res_hist_sig$type=="gen",]
res_hist_nopseudo<-res_hist[res_hist$type=="gen",]

#creamos una función que haga los boxplots
myBoxplot<-function(data,group,titulo,sizesig=0,heightsig=0,axsize,axtextsize,filter=F,comp_vec){
  group <- enquo(group)
  group_col <- ensym(group)  # Capture the group column symbolically
  
  generate_comparisons <- function(levels) {
    combn(levels, 2, simplify = FALSE)
  }#genera todas las posibles comparaciones por pares a partir de los niveles
  
  levels_list <- unique(data[[as.character(group_col)]])  # Extract unique levels
  comparisons <- generate_comparisons(levels_list)
  
  plot<-ggplot(data, aes(x=!!group, y=log2FoldChange, fill=!!group)) +
  geom_boxplot(alpha = 0.80, coef=0,outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width=0.25)+
  geom_point(aes(fill = !!group), size = 2, shape = 21, position = position_jitterdodge()) +
  ggtitle(titulo)+
  ylab(expression("log"[2]*" Fold Change"))+
  theme(axis.title = element_text(size = axsize),text = element_text(size = axtextsize,face = "bold"),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5))
  if(!filter){
    plot<-plot + geom_signif(comparisons = comparisons, test = "t.test",y_position = heightsig,step_increase = 0.1,
                       textsize = sizesig)}
  else{
    plot<-plot +geom_signif(comparisons = list(comp_vec), test = "t.test",y_position = heightsig,step_increase = 0.1,
                   textsize = sizesig)}#este condicional permite quedarse solo con ciertas comparaciones
  
  return(plot)
}

#hacemos los boxplots
type_sig<-myBoxplot(res_hist_sig,type,"Boxplot LFC pseudogenes vs genes de histonas",4,0.4,axtextsize = 10,axsize = 10)
culstering_sig<-myBoxplot(res_hist_sig_nopseudo,clustering,"Boxplot LFC genes de histonas en cluster vs 'solitarias'",4,0.4,axtextsize = 10,axsize = 10)
family_sig<-myBoxplot(res_hist_sig_nopseudo,familia,"Boxplot LFC familias histonas",3,1,filter = T,comp_vec = NULL,axtextsize = 10,axsize = 10)
reg_sig<-myBoxplot(res_hist_sig_nopseudo,reg,"Boxplot LFC genes de histonas dependientes vs independientes de la replicación",4,0.4,axtextsize = 10,axsize = 10)

#guardamos todos los plots
layoutmatrix<-rbind(
    c(1, 2),
    c(3,4)
   )
png(filename = "/home/antotartier/data/histones/senescencia/boxplots_senescencia.png", height =9000 , width =18000 ,res =1100 )
grid.arrange(type_sig,culstering_sig,family_sig,reg_sig,layout_matrix=layoutmatrix)
dev.off()

#el principal regulador de la expresion de las histonas dependientes de la replicación es NPAT, un activador de la transcripción de las mismas que se activa durante la fase S del ciclo celular gracias a la fosforilación del mismo por parte de CDK2/CycE. Dados los cambios en la expresión de las histonas observados procedo a analizar la expresión de estos genes en ambas líneas celulares en busca de una posible causa para la down-regulación observada
load("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/DESeq_conj.obj")
rownames(DESeq)<-gsub("\\.[0-9]+$","",rownames(DESeq))
NPAT<-plotCounts(DESeq, gene="ENSG00000149308", intgroup=c("Passage","Cells"),returnData = T)
NPAT_plot<-ggplot(NPAT, aes(x=Passage,y=count,colour=Cells))+geom_point()+ggtitle("Count plot NPAT")
CDK2<-plotCounts(DESeq, gene="ENSG00000123374", intgroup=c("Passage","Cells"),returnData = T)
CDK2_plot<-ggplot(CDK2, aes(x=Passage,y=count,colour=Cells))+geom_point()+ggtitle("Count plot CDK2")
CycE<-plotCounts(DESeq, gene="ENSG00000105173", intgroup=c("Passage","Cells"),returnData = T)
CycE_plot<-ggplot(CycE, aes(x=Passage,y=count,colour=Cells))+geom_point()+ggtitle("Count plot CycE")

#guardamos todos los plots
png(filename = "~/histones/senescencia/NPAT_CDK2_CycE_cts.png", height =5000 , width =18000 ,res =1100 )
grid.arrange(NPAT_plot,CDK2_plot,CycE_plot,ncol=3)
dev.off()

#los resultados muestran cambios consistentes en ambas líneas celulares en la CDK2, concretamente la downregulación de la misma con un LFC de -1,4. Dado que si está disminye también disminuye la fosforilación de NPAT, este cambio podría ser uno de los responsables de la disminución de la expresión de las histonas. Además, teniendo en cuenta que la disminución es significativamente mayor en las histonas en cluster, que en su mayoría son dependientes de la replicación (aquellas reguladas por NPAT, CDK2 y CycE), esta causa resulta verosimil. No obstante, pese a que la expresión de las CDKs teóricamente permanece constante a lo largo del ciclo celular, para confirmar que este cambio está asociado a la senescencia es necesario descartar que este sea consecuencia del arresto de la proliferación que se produce durante la senescencia. Para ello realizaremos este análisis en datos de arresto por falta de suero en el que se produce un arresto de la proliferación independiente a la senscencia. 


#----------------------------------------------------------------------------------------------------------------------
#Plot TFM
reg_TFM<-myBoxplot(data = res_hist_nopseudo,group = reg,titulo = NULL,sizesig = 4,heightsig = 0.4,axsize = 18,axtextsize = 18)

png(filename = "/home/antotartier/data/histones/senescencia/boxplot_reg.png", height =5000 , width =5000 ,res =1100 )
reg_TFM
dev.off()

fam_TFM<-myBoxplot(data = res_hist_nopseudo,group = familia,titulo = NULL,axsize = 18,axtextsize = 18,filter = T,comp_vec = NULL)

png(filename = "/home/antotartier/data/histones/senescencia/boxplot_fam.png", height =5000 , width =8000 ,res =1100 )
fam_TFM
dev.off()

ggsave("/home/antotartier/data/histones/senescencia/boxplot_reg.pdf",reg_TFM)
ggsave("/home/antotartier/data/histones/senescencia/boxplot_fam.pdf",fam_TFM)
        