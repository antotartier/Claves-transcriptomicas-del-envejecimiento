#!/usr/bin/env python

import HTSeq
import numpy as np
import os.path
from sklearn import linear_model
import argparse
import sys


#handles arguments and saves them in the args object and then calls them along the script
parser = argparse.ArgumentParser()
parser.add_argument('-b','--inbam',required=True,type=str,help='(sorted) bam file with alignments')
parser.add_argument('-i','--introns',required=True,type=str,help='bed file with intron coordinates')
parser.add_argument('-n','--narrow', type=int,default=0, help='number of nucleotides to trim from the ends of the intron (optional, default = 0)')
parser.add_argument('-o','--outf',default=None, help="output file name (optional). Otherwise the .bam extension of input file is replaced with txspeed.tsv")
parser.add_argument('-m','--imin',default=1000,type=int, help="minimum intron size (optional, default 1000 nt)")
args = parser.parse_args()

#function that calculates the coverage, int= intron, which is an interval genomic (hTSEQ, COORDINATES), bamgen (is a generator, gives information from little to little without the need to load the entire BAM, that is it gives one read per iteration until the content is finished)
def get_cov(int, bamgen, narrow=0):
    ints = int.iv.start #start of the intron
    intnp = np.zeros(int.iv.length, dtype='int') #object where the coverage will be stored, a 1D empty numpy array of the length of the intron
    for r in bamgen[int.iv]: #extract the lectures of the intron, all overlapping reads including that ones spliced because the overlap is based only on the coordinates
        if not (r.aligned and r.proper_pair and r.mate_aligned): #of the read given by the generator verify if it is aligned, and a proper pair, if not give me another read and skip the next code chunks
            continue
        if 'N' in [c.type for c in r.cigar]: #eliminate splice reads inspecting cigar, if there is an N in cigar type, the read is spliced
            continue
        #these part of the code corrects the sign of the first pair, an issue due to the nature of pair-end sequencing
        if r.pe_which == 'first': #if the read is pair 1
            r.iv.strand = "+" if r.iv.strand == '-' else '-' #change the sign of the strand, because our library is first stranded
        if r.iv.strand == int.iv.strand: #check if the strand of the read and the intron is the same
            for cig in r.cigar:
                if cig.type == 'M': #select the part of the read aligned with the intron, cyg.type == M means aligned
                    s = max([0,cig.ref_iv.start - ints]) #the start is the start of the cigar block if the substraction of the star of the intron gives a positive number, if not is set to 0 
                    e = max([0,cig.ref_iv.end - ints]) #the end is the end of the cigar block if the substraction of the star of the intron gives a positive number, if not is set to 0. It does not matter if the end of the cygar block is bigger than the end of the intron because in that case numpy automatically interprets it as if it were the end of the intron
                    intnp[s:e] += 1 #add 1 in the positions where the read overlapps with the intron
    if narrow > 0: #if you set narrow  to more than zero, tha narrow is applied
        return(intnp[narrow:-narrow])
    else:
        return(intnp)

#function to calculate the slope of the intron using a lineal regression
def runlm(model,intnp,strand):
    strand = -1 if strand == '-' else 1 #give the value -1 to "-" strand and 1 to "+" strand
    x=np.arange(0,intnp.shape[0],dtype='int') #crate an array of the dimension of the intron, [0,1,2,3,4,...,end of intron], this will be the "x" of the lineal model
    sel = intnp > 0 #an array of true or false, to determine the positions in the coverage which are not zero
    uniqN = np.unique(intnp[sel]).size #extract the unnique values of the coverage which are greater than zero
    if  uniqN < 2:
        res = ['NA','NA','NA','NA',uniqN] #if the coverage array has less than 2 unique values it makes no sanse to calculate the slope so we return NA, if not  we perform the lienal regression 
    else:
        model.fit(x.reshape(-1,1),intnp) # the model function waits for a vertical array so we transform our horizontal array to the vertical one
        res = [model.coef_[0]*strand,model.intercept_] #if strand is positive we multiply per 1 and if not per -1 so we change the sign of the slope to get the real one, store the results in a vector
        model.fit(x[sel].reshape(-1,1),intnp[sel])
        res = res + [model.coef_[0]*strand,model.intercept_] + [uniqN] #we add to results the linal model results only for values greater than zero and the number of unique values
    return("\t".join([str(i) for i in res])) #convert results in a tsv file, we do str(i) because koin functions only with text so we have to transform the numbers in strings
    

def main():
    #Check bam and create generator
    bamgen=HTSeq.BAM_Reader(args.inbam)
    ermsg = "bam file must be sorted by coordinate and indexed (samtools sort / index)"
    try:
        assert bamgen.get_header_dict()['HD']['SO'] == 'coordinate', ermsg #check if the bam is ordered and indexed if not stop and give error message, is imortant to use try because if not the proper check gives error
    except:
        print(ermsg)
        sys.exit(1)
    # Check intron bed file
    introns = HTSeq.BED_Reader(args.introns) #HTSeq it self gives an error if the file is not OK
    
    # output file prefix
    outf = args.inbam
    outf = outf.replace('bam','txspeed.tsv') if args.outf is None else args.outf
    outf = open(outf,'w+') #open the outputfile in write mode (w+)
    outf.write("\t".join(['in_coords','slope','inter','slope_noZ','inter_noZ','unique_noZ','r_mean','r_max']) + '\n') #writing header, "\n" means escape to the next line after writing the line
    
    # Model
    model = linear_model.LinearRegression() #generating the object where the liner regression result will be stored

    # Main loop
    for int in introns:
        if int.iv.length < args.imin: #filter for minimum length of the intron
            continue
        # Calculate coverage
        intnp = get_cov(int,bamgen,args.narrow)
        # Calculate slope
        res = runlm(model,intnp,int.iv.strand)
        #write result
        outf.write("\t".join([int.iv.__str__(),#the coordinates of the intron
                            res,
                            str(intnp.mean()),
                            str(intnp.max())]) +
                            '\n')
    outf.close()

#this line enables you to use this script as a library in python
if __name__ == "__main__":
    main()
