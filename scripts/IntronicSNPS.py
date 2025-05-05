import sys
import collections
import gzip
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --sync input.sync --gff transcriptome.gff --target-length 60 > intron60.sync  "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script identifies SNPs (--sync) located in introns (-gff) of a given length (--target-length)
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--gff", dest="gff", help="A GFF file for the transcriptome")
parser.add_option("--vcf", dest="input", help="A vcf file")
parser.add_option("--target-length", dest="intron", help="The target intron length")
parser.add_option("--output", dest="output", help="OutputFile")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_data(x):
	''' import data either from a gzipped or or uncrompessed file or from STDIN'''
	import gzip         
	if x=="-":
		y=sys.stdin
	elif x.endswith(".gz"):
		y=gzip.open(x,"r")
	else: 
		y=open(x,"r")
	return y


gff=load_data(options.gff)
#gff=load_data("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/Data/Drosophila/dmel-all-r6.57.gff.gz")

vcf=load_data(options.input)
#vcf=load_data("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData/data/PoolSeq2024.vcf.gz")
intron_length_threshold=int(options.intron)
#intron_length_threshold=int(50)
outfile=options.output

intronhash=collections.defaultdict(lambda: collections.defaultdict(list))

for l in gff:
    if l.decode('utf-8').startswith("#"):
        continue
    if len(l.decode('utf-8').split())<8:
        continue
    chrom,source,typ,start,end=l.decode('utf-8').split()[:5]
    if source!="FlyBase":
        continue
    if typ!="intron":
        continue
    if abs(int(end)-int(start))>intron_length_threshold:
        continue
    if int(end)>int(start):
        for i in range(int(start),int(end)+1):
            intronhash[chrom][i]
    else:
        for i in range(int(end),int(start)+1):
            intronhash[chrom][i]

with open(outfile, "w") as o:
    for l in vcf:
        if len(l.decode('utf-8').split("\t"))<2 or l.decode('utf-8').startswith("#"):
            continue
        chrom,pos=l.decode('utf-8').split("\t")[:2]
        if int(pos) in intronhash[chrom]:
            o.write(f"{chrom}\t{pos}\n")
        #print(l.rstrip())