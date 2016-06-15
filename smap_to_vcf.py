
description='''Stand-alone script to convert BioNano smap file format to vcf.'''

vcfheader_orig='''##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">   
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant"> 
##ALT=<ID=DEL,Description="Deletion"> 
##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">  
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">'''

vcfheader_1='''##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant"> 
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">   
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Breakpoint uncertainty for POS"> 
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Breakpoint uncertainty for END"> 
##ALT=<ID=DEL,Description="Deletion"> 
##ALT=<ID=INS,Description="Insertion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'''


import os
import argparse
import math


#create .vcf at smappath (replacing suffix '.smap' with '.vcf')
#currently only processing insertions and deletions
def smap_to_vcf(smappath, refcmap, sample, vcfh) :
    doconf = True #put Confidence in QUAL field
    defaultconf = "." #missing value in vcf format: used if !doconf
    maxconf = "20" #if doconf, this is maximum reported: equivalent to ppv of 0.99
    
    si = ";" #separator for 'INFO' field
    colhead = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" #this is complete if no genotype
    colform = "FORMAT\t%s\n" % sample #if genotype, add this

    outpath = smappath.replace(".smap", ".vcf")
    fout = open(outpath,"w")
    fout.write(vcfh+"\n")

    f1 = open(smappath)
    dogeno = False #smap has genotype info: add to vcf
    nent = 0
    conf = defaultconf #if !doconf, this is used
    for line in f1 :
        if line[0] == "#" :
            continue
        tokens = line.split()
        if len(tokens) < 17 :
            print "ERROR: line incomplete, terminating:\n%s\n" % line
            break
        svtype = tokens[9] 
        if svtype != "insertion" and svtype != "deletion" :
            continue
        if nent == 0 :
            if len(tokens) > 19 : #at least 20 columns for genotype
                dogeno = True
                fout.write(colhead+colform)
            else : #no genotype
                fout.write(colhead+"\n")
        vcftype = "INS" if svtype == "insertion" else "DEL"
        # 0 : SmapEntryID
        #qry = tokens[1] #QryContigID
        ref = int(tokens[2]) #RefcontigID1 -- must be int bc refcmap keys are ints
        # 3 : RefcontigID2 
        qrystart = int(float(tokens[4])) #QryStartPos
        qrystop  = int(float(tokens[5])) #QryEndPos
        refstart = float(tokens[6]) #RefStartPos
        refstop  = float(tokens[7]) #RefEndPos
        if doconf :
            conf = tokens[8] #Confidence
            try :
                conf = float(conf)
                #I was going to put this if block in an else of the try, but this can also catch exceptions raised in log10
                if conf >= 1.0 :
                    conf = maxconf
                else :
                    conf = str(round(-10*math.log10(1-conf),2))
            except :
                conf = defaultconf
        #9: type (above); 10, 11: xmapids; 12: Linkid; 13, 14: Qry Idxs
        refstartidx = int(tokens[15]) #RefStartIdx
        refstopidx  = int(tokens[16]) #RefEndIdx
        if dogeno :
            #tokens[18] = Genotype
            if int(tokens[19]) > 0 : #GenotypeGroup
                if tokens[17] == "homozygous" : #zygosity
                    gt = "1/1"
                else : #note: treating 'unknown' same as 'heterozygous'
                    gt = "1/2"
            else : #note that GenotypeGroup should always be -1 or > 0
                gt = "0/1"
        refstartpos = refcmap[ref][refstartidx-1] #if start site is misresolved
        refstoppos  = refcmap[ref][refstopidx-1] #same for end
        #print refstartpos, refstart, #debug
        #check for the wrong reference:
        if refstartpos > refstop or refstoppos < refstart :
            print "WARNING: incompatible site positons from coordinates (%.1f - %.1f) and indices (%.1f - %.1f): please check reference cmap" % (refstart, refstop, refstartpos, refstoppos)
        if refstart == refstartpos : #no misresolved site
            nextrefstart = refcmap[ref][refstartidx] #list is 0-indexed while the ID in smap is 1-indexed, so +1 is implied
        else : #misresolved, so go to next label
            nextrefstart = refcmap[ref][refstartidx+1]
        #print nextrefstart, ";", refstoppos, refstop, #debug
        pos   = round( (refstart + nextrefstart)/2 )
        cipos = round( (nextrefstart - refstart)/2 )
        if refstoppos == refstop : #no misresolved site on stop
            prevrefstop = refcmap[ref][refstopidx-2] #here, need to take away 2, one for index, one for prev label
        else :
            prevrefstop = refcmap[ref][refstopidx-3]
        #print prevrefstop #debug
        end = round( (prevrefstop + refstop)/2 )
        ciend = round( (refstop - prevrefstop)/2 )
        #note: two inner abs should not be necessary
        svlen = int(round( abs( abs(refstop-refstart) - abs(qrystop-qrystart) ) )) * (1. if svtype == "insertion" else -1.)

        nent += 1
        fout.write(("chr%i\t%i\t.\tN\t<%s>\t%s\tPASS\tSVTYPE=%s"+si+"END=%i"+si+"SVLEN=%i"+si+"CIPOS=-%i,%i"+si+"CIEND=-%i,%i") % (ref, pos, vcftype, conf, vcftype, end, svlen, cipos, cipos, ciend, ciend))
        if dogeno :
            fout.write("\tGT\t%s" % gt)
        fout.write("\n")
    #end loop on input smap
    f1.close()
    fout.close()
    print "Created vcf with %i entries: %s" % (nent, outpath)
#end smap_to_vcf


#return cmap represented by a dict whose keys are integers for refIDs and values are list of positions
def load_cmap(cmappath) :
    cmap = {}
    f1 = open(cmappath)
    refid = -1
    for line in f1 :
        if line[0] == "#" :
            continue
        tokens = line.split()
        if int(tokens[0]) != refid :
            refid = int(tokens[0])
            cmap[refid] = []
        cmap[refid].append( float(tokens[5]) )

    f1.close()
    #print "load_cmap:", len(cmap), ":", " ".join([str(len(cmap[x])) for x in cmap])
    return cmap


#copied from utilities.py
def checkFile(filepath, filesuff="", checkReadable=True) :
    try :
        valid = os.path.isfile(filepath)
        if filesuff != "" :
            valid = (valid and filepath.endswith(filesuff))
        if checkReadable :
            valid = (valid and os.access(filepath, os.R_OK))
    except :
        valid = False
    #print "checkFile:", filepath, valid #debug
    return valid



def getArgs() :    
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-s', dest='smappath'   , help='Path to smap file to convert (required)', type=str) 
    parser.add_argument('-r', dest='refcmappath', help='Path to reference cmap which is reffered to in the smap supplied as -s argument (required)', type=str)
    defsamp = "Sample1"
    parser.add_argument('-n', dest='sample'     , help='Sample ID name for genotype data (optional, default "%s")'%defsamp, type=str, default=defsamp)
    result = parser.parse_args()

    smappath    = result.smappath
    refcmappath = result.refcmappath
    sample      = result.sample

    if not checkFile(smappath,".smap") :
        print "ERROR: smap does not exist, is not readable, or does not end with '.smap':", smappath
        return None

    if not checkFile(refcmappath,".cmap") :
        print "ERROR: reference cmap does not exist, is not readable, or does not end with '.cmap':", refcmappath
        return None

    return smappath, load_cmap(refcmappath), sample


def run_smap_to_vcf():
    getargs = getArgs()
    if getargs != None :
        #smappath, refcmap, sample
        #smap_to_vcf(smappath, refcmap, sample, vcfh=vcfheader_1)
        smap_to_vcf(getargs[0], getargs[1], getargs[2], vcfh=vcfheader_1)


if __name__ == "__main__" :
    run_smap_to_vcf()
