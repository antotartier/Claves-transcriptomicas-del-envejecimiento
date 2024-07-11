library(openxlsx)
library(dplyr)

#cargamos el vsd de DESeq y los id de las histonas y los nombres
load("~/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/vsd_conj.obj")
load("/home/antotartier/data/histones/gene_sets/human/histone_name_id.obj")
sample_info<-read.xlsx("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/sample_info.xlsx")

#sacamos la matriz de counts normalizados y nos quedamos con los de las histonas
counts<-assay(vsd)
rownames(counts)<-gsub("\\.[0-9]+$","",rownames(counts))#eliminamos el punto de los ensemble IDs
counts_hist<-counts[rownames(counts) %in% histone_name_id$ensembl_gene_id,c(1:3,7:12,16:18)]

#sustituimos los id por los nombres
id_to_name<-setNames(histone_name_id$symbol,histone_name_id$ensembl_gene_id)
rownames(counts_hist)<-id_to_name[rownames(counts_hist)]

#sustituimos en las columnas los tubos por la condición
sample_info<-sample_info[c(1:3,7:9,19:21,25:27),]
sample_info <- sample_info %>%
  mutate(
    condition = case_when(
      Passage == "p14" ~ paste0(Cells, "_Temprano"),
      Passage == "p20" ~ paste0(Cells, "_Medio"),
      Passage == "p24" ~ paste0(Cells, "_Tardío"),
      Passage == "p21" ~ paste0(Cells, "_Temprano"),
      Passage == "p27" ~ paste0(Cells, "_Medio"),
      Passage == "p34" ~ paste0(Cells, "_Tardío"),
      TRUE ~ NA_character_
    ),
    ID = paste0("Tube_",ID,".cntTable")
  )#creamos la columna condition
tube_to_condition<-setNames(sample_info$condition,sample_info$ID)
colnames(counts_hist)<-tube_to_condition[colnames(counts_hist)]

#hacemos el heatmap
heatplot<-heatmap(counts_hist,
        scale = "row",  
        margins = c(8, 8),      
        cexRow = 0.5,             # Row text size
        cexCol = 1.5)

png(filename = "/home/antotartier/data/histones/senescencia/heatmap_senescencia.png", width = 9600, height = 12800, res = 1100)
heatplot
dev.off()

pdf(file =  "/home/antotartier/data/histones/senescencia/heatmap_senescencia.pdf",height = 10,width = 10)
heatmap(counts_hist,
                  scale = "row",  
                  margins = c(8, 8),      
                  cexRow = 0.5,             # Row text size
                  cexCol = 1.5)
dev.off()


#extraemos las histonas no detectadas en el experimento, que no se expresan
 hist_noexpres<-histone_name_id[!(histone_name_id$symbol %in% rownames(counts_hist)),]
  
