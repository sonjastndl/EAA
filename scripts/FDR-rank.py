import sys
from collections import defaultdict
from bisect import bisect
import math
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input input.fet --true -1 -permuted -2,-3,-4 > output.fdr"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

Calculate False Discovery Rate (FDR) based on the approach of Jha et al. 2015, which requires a "true" pvalue and several "permuted" pvalues based on permuations of the labels in the statistical model. The average value of these permutations per SNP are then used to estimate the distribution of permuted p-values. The parameter (--true) indicates the position of the column in the input that contains the "true" pvalue. Since this is pythonic, this index is 0-based!! Likewise (--permuted) is a comma-separated list of the positions of "permuted" p-values.
""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="input file with true and permuted P-values")
parser.add_option("--true", dest="true", help="poistion of 'true' pvalue in dataset",default="No")
parser.add_option("--permuted", dest="permuted", help="poistion of 'permuted' pvalues in dataset, comma-separated list")
parser.add_option("--output", dest="out", help="output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


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

data=options.input
out=open(options.out,"w")
true=int(options.true)
permuted=map(int,options.permuted.split(","))

code={"2L":"2","2R":"3","3L":"4","3R":"5","4":"6","X":"1"}
out.write("CHR\tBP\tP\tQ\tgene\n")



## New calucaltion of permuted minimum based on ranks:

## read all permuted p-values:
permutedlist=defaultdict(list)
TRUE=[]


for l in load_data(data):
    a=l.split()
    if "NA" in l:
        continue
    TRUE.append(float(a[true]))
    for i in range(len(permuted)):
        permutedlist[i].append(float(a[permuted[i]]))

Permuted=[min(y) for y in zip(*[sorted(x) for x in permutedlist.values()])]

print("read data done")
## calculate flase discovery rate according to Jha et al 2015:

FDR=defaultdict(float)

TRUES=sorted(TRUE)
#print TRUES[:5]
#print Permuted[:5]
rlist=[]
for k in range(len(TRUES)):
    trueval=TRUES[k]
    i=bisect(Permuted,trueval)
    FDR[trueval]=i/float(k+1)
    rlist.append(i/float(k+1))

print (rlist[:10])
print("calculate FDR done")

QVAL=defaultdict(float)
qlist=[]

for i in reversed(range(len(TRUES))):
    if i==len(TRUES)-1:
        QVAL[TRUES[i]]=rlist[i]
        qlist.append(rlist[i])
        #print i,TRUES[i],rlist[i],rlist[i]
    else:
        if rlist[i]>=QVAL[TRUES[i+1]]:
            QVAL[TRUES[i]]=QVAL[TRUES[i+1]]
            qlist.append(QVAL[TRUES[i+1]])
            #print i,TRUES[i],rlist[i],QVAL[TRUES[i+1]]
        else:
            QVAL[TRUES[i]]=rlist[i]
            qlist.append(rlist[i])
            #print i,TRUES[i],rlist[i],rlist[i]

print ("qvalues calculated")

## append as last column
with open("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/permutations_geo_control/merged_permutation_pvalues_true_q.csv", "w") as out:
    for l in load_data(data):
        a=l.split()
        if "NA" in l:
            continue
        if a[0] not in code:
            continue
        R=FDR[float(a[true])]
        q=QVAL[float(a[true])]
        if q>1:
            q=1
        #out.write(l.rstrip()+"\t"+str(R)+"\t"+str(q)+"\n")
        out.write(a[0]+"\t"+a[1]+"\t"+a[true]+"\t"+str(q)+"\n")

# r('png("'+sys.argv[4]+'.png",width=1000,height=500)')
# r.assign("pvalT",robjects.vectors.FloatVector(TRUES))
# r.assign("pvalD",robjects.vectors.FloatVector(Permuted))
# r.assign("qval",robjects.vectors.FloatVector(qlist[::-1]))
#
# r('plot(pvalT,type="l",col="black",ylab="p-val")')
# r('points(pvalD,type="l",col="red")')
# r('points(qval,type="l",col="blue")')
# r('legend("topleft",legend=c("pval-TRUE","pval-DRIFT","qval"),lty=1,col=c("black","red","blue"))')
# r('dev.off()')
