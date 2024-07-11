require(rtracklayer)
require(GenomicFeatures)
#set outpath
outpath<-'/home/antotartier/data/velocidad_transcrip/LAKI/rec.anot_res/'

# Load seqinfo (to avoid warnings when doing intersections because some chrom names are not present)
seqin <- read.table('/data/genomes/GRCm39_gencode/STAR_gencodeM34_index/chrNameLength.txt')
seqin <- Seqinfo(seqin$V1,seqin$V2)

# Generate introns the way they do but using built-in genomicranges functions
require(rtracklayer)
gtf <- import.gff('/data/genomes/GRCm39_gencode/gencode.vM32.primary_assembly.annotation.gtf')
# Remove exons and UTRs from the gene range
genes <- gtf[gtf$type == 'gene']
names(genes) <- genes$gene_id
introns <- subtract(genes,gtf[gtf$type %in% c('exon','UTR')],ignore.strand=F)
# Reduce (like betools merge)
introns <- GenomicRanges::reduce(unlist(introns))
saveRDS(introns,paste0(outpath,'gencode.vM32.introns.rds'))


# Load the sample info
sample_info<-fread("/home/antotartier/data/20231017_RNASeq_LAKI/Sample_info.tsv")
sample_info<-sample_info %>% mutate(genotipo=str_extract(SampleID,"KO|WT"),
                            sexo=ifelse(grepl("M",SampleID),"M","H"))
#eliminate ausent samples and FACE samples
ausent <- c("Tube_39","Tube_56","Tube_57","Tube_66","Tube_87") #error while transfering data
FACE <- c("Tube_115","Tube_116","Tube_117","Tube_118","Tube_119","Tube_120")
sample_info<-sample_info[!grepl(paste(ausent,collapse = "|"),sample_info$SampleName)]
sample_info<-sample_info[!grepl(paste(FACE,collapse = "|"),sample_info$SampleName)]



require(data.table)
require(GenomicRanges)
# Load the junctions into a data.frame
star.files <- list.files('/home/antotartier/data/20231017_RNASeq_LAKI/newSTAR_aligns/',pattern='*SJ.out.tab',full.names=T)
#eliminate the erroneous files and the FACE files
star.files <- star.files[!grepl(paste(ausent,collapse = "|"),star.files)]
star.files <- star.files[!grepl(paste(FACE,collapse = "|"),star.files)]

#Here you can apply filters to select certain samples
# star.files<-star.files[grepl(paste0(unlist(as.vector(sample_info[sample_info$Tissue=="Ileon",1])),"_",collapse = "|"),star.files)]
# star.files<-star.files[grepl(paste0(unlist(as.vector(sample_info[sample_info$Tissue=="Kidney",1])),"_",collapse = "|"),star.files)]
# star.files<-star.files[grepl(paste0(unlist(as.vector(sample_info[sample_info$Tissue=="Muscle",1])),"_",collapse = "|"),star.files)]
# star.files<-star.files[grepl(paste0(unlist(as.vector(sample_info[sample_info$Tissue=="Liver",1])),"_",collapse = "|"),star.files)]
# star.files<-star.files[grepl(paste0(unlist(as.vector(sample_info[sample_info$Tissue=="Colon",1])),"_",collapse = "|"),star.files)]
star.files<-star.files[grepl(paste0(unlist(as.vector(sample_info[sample_info$Tissue=="Heart",1])),"_",collapse = "|"),star.files)]

#reading the junctions
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
known.jun <- readRDS('/home/antotartier/data/velocidad_transcrip/LAKI/rec.anot_res/gencode.vM32.introns.rds')

# Extract STAR junctions that are equal or within the predicted introns
eq_se <- subsetByOverlaps(known.jun,star.jun,type='equal',maxgap = 2)
eq_s <- subsetByOverlaps(known.jun,star.jun,type='start',maxgap = 2)
eq_e <- subsetByOverlaps(known.jun,star.jun,type='end',maxgap = 2)
dif <- subsetByOverlaps(known.jun,star.jun,type='within')
sel_in <- disjoin(c(dif,eq_se,eq_s,eq_e))
# Filter out short introns
sel_in <- sel_in[width(sel_in) >= 1000]
# Export the introns
rtracklayer::export.bed(sel_in,'/home/antotartier/data/velocidad_transcrip/LAKI/rec.anot_res/selected_introns_Heart.bed')





#plot of the result for actin, use sel_in before filtering
options(ucscChromosomeNames=FALSE)
actb <- GRanges('chr7',IRanges(5592816,5606655))
exns <- gtf[gtf$type == 'exon']
exnsactb <- subsetByOverlaps(exns, actb)

plotTracks(list(
  AnnotationTrack(subsetByOverlaps(exnsactb, actb),name='exonsbytx',group=exnsactb$transcript_id),
  AnnotationTrack(subsetByOverlaps(known.jun,actb,type = 'any'),name='known',fill='orange'),
  AnnotationTrack(subsetByOverlaps(star.jun, actb),name='STAR junctions'),
  AnnotationTrack(subsetByOverlaps(c(dif,eq_se,eq_s,eq_e), actb),name='Selected',fill='blue'),
  AnnotationTrack(subsetByOverlaps(sel_in, actb),name='Disjoin',fill='darkred')))