import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    """ import data either from a gzipped or or uncrompessed file or from STDIN"""
    import gzip

    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y

def GetIndexOfMostCommonAltAllel(pops,formID):
    AFs=d(int)
    ALLCOUNT=0
    ## loop through all populations
    for n in pops:
        pophash = dict(zip(formID, n.split(":")))
        if pophash["GT"] == "./.":
            continue
        ## append the allele frequencies to dictionary
        for i in range(len(pophash["GT"].split("/"))):
            ## get ID of allele
            ALLELE = pophash["GT"].split("/")[i]
            ## ALLELE 0 == Ref
            if ALLELE=="0":
                ACOUNT = pophash["RD"]
            else:
                ## else loop through alternative alleles 
                ACOUNT = pophash["AD"].split(",")[i-1]

            AFs[int(ALLELE)]+=float(ACOUNT)
        ## sum up counts accross all pops
        ALLCOUNT+=int(pophash["DP"])

    ## make list of average AFs
    AFL = []
    for k, v in sorted(AFs.items()):
        if k==0:
            continue
        AFL.append(v/ALLCOUNT)
    
    ## return the index of the most common allele (+1 for the alternative allel code in the GT)
    return str(AFL.index(max(AFL))+1)

for l in load_data(options.IN):
    ## skip primary header
    if l.startswith("#"):
        continue
    ## check if VCF file correct
    if len(l.split("\t")) < 9:
        print("ERRROR IN VCF FILE, CHECK CAREFULLY!!")
        sys.exit()
    ## split and name first 9 columns
    chr, pos, ID, ref, alt, qual, Filter, Info, Format = l.split()[:9]
    ## get Fields in Format field
    formID = Format.split(":")
    ## ignore positions with In/Dels
    if "DELETION" in Info or "INSERTION" in Info:
        continue
    ## get alternative alleles 
    ALT = alt.split(",")
    ## get population fields
    pops = l.split()[9:]
    ## get index of most common alternative allel if there are more than two alleles
    if len(ALT)>1:
        ID=GetIndexOfMostCommonAltAllel(pops,formID)
    else:
        ID="1"
    totalsync=[]
    ## loop through all populations and get 
    for n in pops:
        pophash = dict(zip(formID, n.split(":")))
        ## account for missing data
        if pophash["GT"] == "./.":
            totalsync.append(".,.")
            continue
        ## Get allele counts of REF and most common alternative allele
        GTcount=d(int)
        for i in range(len(pophash["GT"].split("/"))):
            ## get ID of allele
            ALLELE = pophash["GT"].split("/")[i]
            ## ALLELE 0 == Ref
            if ALLELE =="0":
                continue
            else:
                GTcount[ALLELE]= pophash["AD"].split(",")[int(i)-1]
        if ID in GTcount:
            totalsync.append(pophash["RD"]+","+GTcount[ID])
        else:
            totalsync.append(pophash["RD"]+",0")
    ## print the allele count matrix 
    print(chr+"\t"+pos+"\t"+ref+","+ALT[int(ID)-1]+"\t"+"\t".join(totalsync))