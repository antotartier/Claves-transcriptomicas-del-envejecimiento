require(rtracklayer)
require(GenomicFeatures)


# Load seqinfo (to avoid warnings when doing intersections because some chrom names are not present)
seqin <- read.table('/data/genomes/GRCh38_gencode/STAR_index/chrNameLength.txt')
seqin <- Seqinfo(seqin$V1,seqin$V2)

# Generate introns the way they do but using built-in genomicranges functions
require(rtracklayer)
gtf <- import.gff('/data/genomes/GRCh38_gencode/gencode.v44.primary_assembly.annotation.gtf')
# Remove exons and UTRs from the gene range
genes <- gtf[gtf$type == 'gene']
names(genes) <- genes$gene_id
introns <- subtract(genes,gtf[gtf$type %in% c('exon','UTR')],ignore.strand=F)
# Reduce (like betools merge)
introns <- GenomicRanges::reduce(unlist(introns))
saveRDS(introns,'/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/rec.anot_res/gencode.v44.introns.rds')


# Load the sample info
sample_info<-read.xlsx("/home/antotartier/data/transposones/fibrob_senescentes/TEtranscripts/DESeq_fib_senescentes/sample_info.xlsx")


require(data.table)
require(GenomicRanges)
# Load the junctions into a data.frame
star.files <- list.files('/home/antotartier/data/fib_senescentes_data/STAR_aligns/',pattern='*SJ.out.tab',full.names=T)

#aplaying filters to select the samples from the GM00038 cell line
star.files<-star.files[grepl("Tube_[0-9]{2}_",star.files)]
#star.files<-star.files[grepl("Tube_[0-9]{1}_",star.files)] #this would be executed instead for the GM05565 cell line

star.jun <- lapply(star.files, fread)
star.jun <- do.call('rbind',star.jun)
# Filter by number of uniquely mapped reads (V7) and overhang length (V9)
star.jun <- star.jun[V9 > 20 & V7 > 5]
# Collapse equal ranges and sum up reads. Get also in how many samples the junction is detected (.N)
star.jun <- star.jun[,.(.N, reads=sum(V7)),by=.(V1,V2,V3,V4)]
# Arrange columns
colnames(star.jun)<-c('chr','start','end','strand','N','reads')
# Transform strand to +/- instead of 0,1,2
star.jun[,strand:=factor(strand,levels=c(0,1,2),labels=c('*','+','-'))]
# Keep junctions supported in all samples
star.jun <- star.jun[N == length(star.files)]
# Convert to GRanges
star.jun <- makeGRangesFromDataFrame(star.jun,keep.extra.columns = T,seqinfo = seqin)

# Load gencode junctions
known.jun <- readRDS('/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/rec.anot_res/gencode.v44.introns.rds')

# Extract STAR junctions that are equal or within the predicted introns
eq_se <- subsetByOverlaps(known.jun,star.jun,type='equal',maxgap = 2)
eq_s <- subsetByOverlaps(known.jun,star.jun,type='start',maxgap = 2)
eq_e <- subsetByOverlaps(known.jun,star.jun,type='end',maxgap = 2)
dif <- subsetByOverlaps(known.jun,star.jun,type='within')
sel_in <- disjoin(c(dif,eq_se,eq_s,eq_e))
# Filter out short introns
sel_in <- sel_in[width(sel_in) >= 1000]
# Export the introns
rtracklayer::export.bed(sel_in,'/home/antotartier/data/velocidad_transcrip/fibrob_senescentes/rec.anot_res/selected_introns_GM00038.bed')




#this is just an example to facilitate the comprehension of the code
#plot of the result for actin, use sel_in before filtering
library(Gviz)
options(ucscChromosomeNames=FALSE)
actb <- GRanges('chr7',IRanges(5592816,5606655))
exns <- gtf[gtf$type == 'exon']
exnsactb <- subsetByOverlaps(exns, actb)

pdf(file = "~/data/velocidad_transcrip/fibrob_senescentes/Tracks.pdf")
plotTracks(list(
  AnnotationTrack(subsetByOverlaps(exnsactb, actb),name='Isoformas actina',group=exnsactb$transcript_id,cex.title=1),
  AnnotationTrack(subsetByOverlaps(known.jun,actb,type = 'any'),name='Intrones',fill='orange',cex.title=1),
  AnnotationTrack(subsetByOverlaps(star.jun, actb),name='STAR junctions',fill="green",cex.title=1),
  AnnotationTrack(subsetByOverlaps(c(dif,eq_se,eq_s,eq_e), actb),name='Seleccionados',fill='blue',cex.title=1,ttitleJustify="center"),
  AnnotationTrack(subsetByOverlaps(sel_in, actb),name='Fusionados',fill='darkred',cex.title=1,just.group="left")))
dev.off()
