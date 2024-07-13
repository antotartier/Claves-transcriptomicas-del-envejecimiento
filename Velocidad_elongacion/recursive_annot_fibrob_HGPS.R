require(rtracklayer)
require(GenomicFeatures)


# Load seqinfo (to avoid warnings when doing intersections because some chrom names are not present)
seqin <- read.table('/data/genomes/GRCh38_gencode/STAR_index/chrNameLength.txt')
seqin <- Seqinfo(seqin$V1,seqin$V2)

require(data.table)
require(GenomicRanges)

# Load the sample info
sample_info<-fread("/home/antotartier/data/fib_HGPS/FASTQ/Sra_meta.tsv")

# Load the junctions into a data.frame
star.files <- list.files('/home/antotartier/data/fib_HGPS/STAR_aligns/',pattern='*SJ.out.tab',full.names=T)
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

# Load gencode junctions, they were previously generated in the senescent fibroblasts script
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
rtracklayer::export.bed(sel_in,'/home/antotartier/data/velocidad_transcrip/fib_HGPS/rec.anot_res/selected_introns_fib_HGPS.bed')
