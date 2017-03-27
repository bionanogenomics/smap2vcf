#VCF converter version 2

#Version 1 supports conversion of insertions and deletions
#Version 2 supports conversion of translocation breakpoints

description='''Stand-alone script to convert BioNano smap file format to vcf.'''

#Template VCF headers
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

#Version 1 headers
vcfheader_1='''##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Breakpoint uncertainty for POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Breakpoint uncertainty for END">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'''

#Version 2 headers
vcfheader_2='''##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Breakpoint uncertainty for POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Breakpoint uncertainty for END">
##INFO=<ID=CT,Number=1,Type=String,Description="Breakpoint connection type">
##INFO=<ID=EXPERIMENT,Number=1,Type=Integer,Description="experiment_id from dbVar submission of the experiment that generated this call">
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="sample_id from dbVar submission. Each call must have only one of either SAMPLE or SAMPLESET">
##INFO=<ID=SAMPLESET,Number=1,Type=Integer,Description="sampleset_id from dbVar submission. Each call must have only one of either SAMPLE or SAMPLESET">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=TRA,Description="Translocation breakpoint">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'''

import os
import argparse
import datetime
import math

#define accepted SV types and their corresponding converted types in VCF
class SV(object):
    def __init__(self, smap_type, vcf_type, svlen=None):
        self.smap_type = smap_type
        self.vcf_type = vcf_type
        self.svlen = svlen #default SV length

sv_types = [SV("insertion", "INS"), SV("deletion","DEL"), SV("translocation_intrachr","TRA",0),SV("translocation_interchr","TRA",0)] #conversion table
accepted_svtypes = dict([(sv.smap_type,sv.vcf_type) for sv in sv_types])
default_svlen = dict([(sv.smap_type,sv.svlen) for sv in sv_types])

#utility functions

#read xmap; returns xmap entries in a list
def readxmap(inxmap):
    xmap_fh = open(inxmap,'r')

    xmap_list = []

    for line in xmap_fh:
        if line.startswith("#"):
            continue

        xmap_list.append(line.strip())

    xmap_fh.close()

    return xmap_list

#read cmap; returns cmap represented by a dict whose keys are integers for refIDs and values are list of positions
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

#check file; return bool (copied from utilities.py from assembly pipeline)
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

#main parser function
def smap_to_vcf(smappath, refcmap, xmappath, sample, output_prefix, ref_accession, human_bool, exp_id, vcfh) :

    doconf = True #put Confidence in QUAL field
    defaultconf = "." #missing value in vcf format: used if !doconf
    maxconf = "20" #if doconf, this is maximum reported: equivalent to ppv of 0.99

    si = ";" #separator for 'INFO' field
    colhead = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" #this is complete if no genotype
    colform = "FORMAT\t%s\n" % sample #if genotype, add this

    outpath = smappath.replace(".smap", ".vcf") #create default .vcf at smappath (replacing suffix '.smap' with '.vcf')
    if output_prefix is not None:
        outpath = output_prefix+".vcf"

    today = datetime.date.today() #prepare fileDate for header

    fileDate = "##fileDate=%s" % today

    fout = open(outpath,"w")
    fout.write(fileDate+"\n")
    fout.write(vcfh+"\n")

    refcmappath = refcmap #keep path
    refcmap = load_cmap(refcmap) #read cmap

    xmap_list = readxmap(xmappath) #read xmap

    reference_header = "##reference=%s\n" % ref_accession
    reference_contigids = sorted([key for key in refcmap])

    if human_bool:
        reference_contigids[22] = 'X'
        reference_contigids[23] = 'Y'
        contig_header = ["##contig=<ID=chr%s>\n" % elem for elem in reference_contigids]
    else:
        contig_header = ["##contig=<ID=chr%s>\n" % elem for elem in reference_contigids]

    f1 = open(smappath)
    dogeno = False #smap has genotype info: add to vcf
    nent = 0 #number of processed entries
    conf = defaultconf #if !doconf, this is used
    for line in f1 :
        if line[0] == "#" :
            if '# Reference Maps From:' in line:
                tokens = line.split("\t")
                ref = tokens[1] #get path to reference cmap
            continue
        tokens = line.split("\t")
        if len(tokens) < 17 :
            print "ERROR: line incomplete, terminating:\n%s\n" % line
            break

        svtype = tokens[9]
        if svtype not in accepted_svtypes: #only process certain SV types
            continue
        vcftype = accepted_svtypes[svtype]

        if nent == 0 :
            if len(tokens) > 19 : #at least 20 columns for genotype
                dogeno = True
                fout.write(reference_header)
                fout.write("".join(contig_header))
                fout.write(colhead+colform)
            else : #no genotype
                fout.write(colhead+"\n")

        smapid = tokens[0] #SmapEntryID
        #qry = tokens[1] #QryContigID
        ref = int(tokens[2]) #RefcontigID1 -- must be int bc refcmap keys are ints
        ref1 = int(tokens[3]) #RefcontigID2 -- must be int bc refcmap keys are ints
        qrystart = int(float(tokens[4])) #QryStartPos
        qrystop  = int(float(tokens[5])) #QryEndPos
        refstart = float(tokens[6]) #RefStartPos
        refstop  = float(tokens[7]) #RefEndPos

        if doconf :
            conf = tokens[8] #Confidence
            try :
                conf = float(conf)
                #I was going to put this if block in an else of the try, but this can also catch exceptions raised in log10

                if conf==-1.0 :
                    conf = defaultconf
                elif conf >= 1.0 :
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
                    gt = "." #"1/2"
            else : #note that GenotypeGroup should always be -1 or > 0
                gt = "0/1"

        refstartpos = refcmap[ref][refstartidx-1] #if start site is misresolved
        refstoppos  = refcmap[ref1][refstopidx-1] #same for end
        #print refstartpos, refstart, #debug

        #check for the wrong reference for insertions and deletions
        if svtype=="insertion" or svtype=="deletion":
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
                prevrefstop = refcmap[ref1][refstopidx-2] #here, need to take away 2, one for index, one for prev label
            else :
                prevrefstop = refcmap[ref1][refstopidx-3]
            #print prevrefstop #debug
            end = round( (prevrefstop + refstop)/2 )+1
            ciend = round( (refstop - prevrefstop)/2 )

            if svtype=="insertion": #special case where ref start and end are the same for insertions
                if refstart==refstop:
                    pos,end = end,pos
                    cipos,ciend = ciend,cipos

            if pos>end:
                print "WARNING: start position (%.1f) larger than end position (%.1f)" % (pos, end)
                print line.strip()
                print refstartpos, refstart, refstoppos, refstop
                print refstart, nextrefstart, prevrefstop, refstop
            if pos==end:
                print "WARNING: start position (%.1f) is equal to end position (%.1f)" % (pos, end)

        elif svtype=="translocation_intrachr" or svtype=="translocation_interchr":

            #get orientation information
            xmapid1 = int(tokens[10]) #XmapID1
            xmapid2 = int(tokens[11]) #XmapID2

            xmap1 = xmap_list[xmapid1-1].split('\t')
            xmap2 = xmap_list[xmapid2-1].split('\t')

            if (xmapid1!=int(xmap1[0])) or (xmapid2!=int(xmap2[0])):
                print "WARNING: xmap ids are not consistent"

            xmap1_refstart = float(xmap1[5])
            xmap2_refstart = float(xmap2[5])
            xmap1_refend = float(xmap1[6])
            xmap2_refend = float(xmap2[6])

            if refstart==xmap1_refstart:
                orientation1 = "+"
            else:
                if refstart==xmap1_refend:
                    orientation1 = "-"
                else:
                    print "WARNING: breakpoint start position not consistent with xmap data"

            if refstop==xmap2_refstart:
                orientation2 = "+"
            else:
                if refstop==xmap2_refend:
                    orientation2 = "-"
                else:
                    print "WARNING: breakpoint start position not consistent with xmap data"

            #assign orientation to breakpoint
            if orientation1=="+" and orientation2=="+":
                orientation="3to3"
            elif orientation1=="+" and orientation2=="-":
                orientation="3to5"
            elif orientation1=="-" and orientation2=="+":
                orientation="5to3"
            elif orientation1=="-" and orientation2=="-":
                orientation="5to5"
            else:
                print "WARNING: orientation combination not recognized"

            if refstart == refstartpos : #no misresolved site
                prevrefstart = refcmap[ref][refstartidx-2]
                nextrefstart = refcmap[ref][refstartidx] #list is 0-indexed while the ID in smap is 1-indexed, so +1 is implied
            else : #misresolved, so go to next label
                prevrefstart = refcmap[ref][refstartidx-3]
                nextrefstart = refcmap[ref][refstartidx+1]
            #print nextrefstart, ";", refstoppos, refstop, #debug
            #pos   = round( (prevrefstart + nextrefstart)/2 )
            pos = refstart
            cipos = round( (nextrefstart - prevrefstart)/2 )

            if refstoppos == refstop : #no misresolved site on stop
                prevrefstop = refcmap[ref1][refstopidx-2] #here, need to take away 2, one for index, one for prev label
                nextrefstop = refcmap[ref1][refstopidx]
            else :
                prevrefstop = refcmap[ref1][refstopidx-3]
                nextrefstop = refcmap[ref1][refstopidx+1]
            #print prevrefstop #debug
            #end = round( (prevrefstop + nextrefstop)/2 )
            end = refstop
            ciend = round( (nextrefstop - prevrefstop)/2 )

        svlen = default_svlen[svtype]
        if svlen is None:
            #note: two inner abs should not be necessary
            svlen = int(round( abs( abs(refstop-refstart) - abs(qrystop-qrystart) ) )) * (1. if svtype == "insertion" else -1.)

        nent += 1

        if human_bool:
            if ref==23:
                ref='X'
            elif ref==24:
                ref='Y'

            if ref1==23:
                ref1='X'
            elif ref1==24:
                ref1='Y'

        if svtype=="insertion" or svtype=="deletion":
            fout.write(("chr%s\t%i\tSMAP%s\tN\t<%s>\t%s\tPASS\tSVTYPE=%s"+si+"END=%i"+si+"SVLEN=%i"+si+"CIPOS=-%i,%i"+si+"CIEND=-%i,%i"+si+"EXPERIMENT=%s") % (ref, pos, smapid,vcftype, conf, vcftype, end, svlen, cipos, cipos, ciend, ciend, exp_id))
        elif svtype=="translocation_intrachr" or svtype=="translocation_interchr":
            fout.write(("chr%s\t%i\tSMAP%s\tN\t<%s>\t%s\tPASS\tSVTYPE=%s"+si+"CHR2=chr%s"+si+"END=%i"+si+"SVLEN=%i"+si+"CIPOS=-%i,%i"+si+"CIEND=-%i,%i"+si+"CT=%s"+si+"EXPERIMENT=%s") % (ref, pos, smapid,vcftype, conf, vcftype, ref1, end, svlen, cipos, cipos, ciend, ciend, orientation, exp_id))
        if dogeno :
            fout.write("\tGT\t%s" % gt)
        fout.write("\n")

    #end loop on input smap
    f1.close()
    fout.close()

    print "Created vcf with %i entries: %s" % (nent, outpath)
#end smap_to_vcf

def getArgs() :
    parser = argparse.ArgumentParser(description=description)

	#required input
    parser.add_argument('-s', dest='smappath', help='Path to smap file to convert (required)', type=str)
    parser.add_argument('-r', dest='refcmappath', help='Path to reference cmap which is refered to in the smap supplied as -s argument (required)', type=str)
    parser.add_argument('-x', dest='xmappath', help='Path to xmap which is refered to in the smap supplied as -s argument (required)', type=str)

    #optional input
    defsamp = "Sample1"
    parser.add_argument('-n', dest='sample', help='Sample ID name for genotype data (optional, default "%s")'%defsamp, type=str, default=defsamp)
    output_prefix = None
    parser.add_argument('-o', dest='output_prefix', help='Prefix for output vcf (optional, default to be same as input smap)', type=str, default=output_prefix)
    ref_accession = "GCA_000001405.1"
    parser.add_argument('-a', dest='ref_accession', help='RefSeq assembly accession version (optional, default "%s")'%ref_accession, type=str, default=ref_accession)
    human_bool = True
    parser.add_argument('-b', dest='human_bool', help='Whether sample is human (optional, default "%s")'%human_bool, type=str, default=human_bool)
    exp_id = '1'
    parser.add_argument('-i', dest='exp_id', help='dbVar-required experiemnt ID (optional, default "%s")'%exp_id, type=str, default=exp_id)

    result = parser.parse_args()

    #unpack
    smappath = result.smappath #required
    refcmappath = result.refcmappath
    xmappath = result.xmappath

    sample = result.sample #optional
    output_prefix = result.output_prefix
    ref_accession = result.ref_accession
    human_bool = result.human_bool
    exp_id = result.exp_id

    #checking files
    if not checkFile(smappath,".smap") :
        print "ERROR: smap does not exist, is not readable, or does not end with '.smap':", smappath
        return None

    if not checkFile(refcmappath,".cmap") :
        print "ERROR: reference cmap does not exist, is not readable, or does not end with '.cmap':", refcmappath
        return None

    if not checkFile(xmappath,".xmap") :
        print "ERROR: xmap does not exist, is not readable, or does not end with '.xmap':", xmappath
        return None

    return smappath, refcmappath, xmappath, sample, output_prefix, ref_accession, human_bool, exp_id

def run_smap_to_vcf():
    getargs = getArgs()
    if getargs != None :
        smap_to_vcf(getargs[0], getargs[1], getargs[2], getargs[3], getargs[4], getargs[5], getargs[6], getargs[7], vcfh=vcfheader_2)

if __name__ == "__main__" :
    run_smap_to_vcf()
