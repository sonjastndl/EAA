#!/bin/sh

################################################# CREATING REQUIRED DIRECTORIES TO PROCESS AND STORE OUTPUTS ##################################################
#scriptdir="/home/sonjastndl/s3/LGA/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts"
DIRECTORY=$(pwd)

while getopts d:o:s:i:m:w:v:b:? opt 2>/dev/null
do
  case $opt in
    d) DIRECTORY=$OPTARG;;
    o) OUTPUTDIR=$OPTARG;;
    s) SAMPLES=$OPTARG;;
    i) VCFFILE=$OPTARG;;
    m) METADATA=$OPTARG;;
    w) WORLDCLIM=$OPTARG;;
    v) VARIABLE=$OPTARG;;
    b) BIOVAR=$OPTARG;;
    ?) echo "Valid parameters are: [-l, -s]" ;;
  esac
done

scriptdir="${LGAdir}/scripts"
wd="$1"
continent="$2"
samplelist="$3"
input="$4"
metadata="$5"
envdata="$6"


resultsdir="${wd}/results"
summary="${resultsdir}/summary"
data= "${wd}/data"

mkdir $wd
mkdir $resultsdir
mkdir $summray
mkdir $data

cd $wd
################################################## DOWNLOAD VCF FROM DEST.bio AND EXTRACT SAMPLENAMES  ##################################################
echo "DOWNLOADING VCF FROM DEST.bio AND EXTRACTING SAMPLENAMES"

wget --tries=inf "http://berglandlab.uvadcos.io/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz"
wget "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv"
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz


## Extract the information on the available samples
awk '{FS=","}{if (NR!=1) {print $1}}' dest_v2.samps_3May2024.csv > ${data}/samplenames.csv
cp dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz ${data}/PoolSeq2024.vcf.gz
mv  dmel-all-r6.57.gff.gz > ${data}/dmel-all-r6.57.gff.gz


module load Tools/vcftools_0.1.13

#################################################################### NAMING THE ANALYSIS ###################################################################################



awk -F, 'NR > 1 && $6 == "Europe" {print $1}' dest_v2.samps_3May2024.csv > ${data}/EuropeSamples.csv
awk -F, -v var="$continent" 'NR > 1 && $6 == var && $46 == "Pass" {print $1}' dest_v2.samps_3May2024.csv > ${data}/EuropeSamples_Pass.csv
awk -F "," '$(NF-7) !="Pass" || $(NF-9)<15 {print $1"\t"$(NF-7)"\t"$(NF-9)}' dest_v2.samps_3May2024.csv > ${data}/REMOVE.ids


echo "START FILTERING"

vcftools --gzvcf $input --keep $samplelist  --minDP 15 --stdout --recode-INFO-all --recode | grep -v "\./\." > ${resultsdir}/Subsampled_DP15.vcf.gz

pigz -dc $input | awk '$0~/^#/ || length($5)==1' | vcftools --vcf - \
        --keep $samplelist \
        --remove ${data}/REMOVE.ids \
        --stdout --recode-INFO-all \
        --recode | grep -v "\./\." | gzip  >  ${resultsdir}/Subsampled_recode_DP15.vcf.gz

echo "DONE WITH FILTERING"


Sub2="${resultsdir}/Subsampled_recode2_DP15.vcf.gz"
Sub3="${resultsdir}/Subsampled_recode3_DP15.vcf.gz"
Sub4="${resultsdir}/Subsampled_final_DP15.vcf.gz"


#################################################### RANDOMLY PICK n LINES FROM VCF ###########################################################
##
##

#echo "SUBSAMPLING WITH PYTHON"
#echo $Sub2
#
#
#
#python3 ${scriptdir}/SubsampleVCF.py \
#    --input ${Sub2} \
#    --snps all \
#    --output ${Sub3}
#
#module load Tools/bcftools-1.16 
###

###################################################    GET NEUTRAL SNPS   #######################################################################


bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $Sub3 | gzip > $Sub4
echo "GETTING INTRONIC SNPs"
get intronic SNPS

### DOWNLOAD GFF GILE ###
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz
gff_file="${data}/dmel-all-r6.57.gff.gz"
neutralSNPs="${results}/Subsampled_NeutralSNPS_80.tsv"
echo "FILTERING INTRONIC SNPs"

python3 ${scriptdir}/IntronicSNPS.py \
 --gff $gff_file \
 --vcf ${Sub4} \
 --target-length 80 \
 --output $neutralSNPs


vcftools --gzvcf $Sub4 \
   --positions $neutralSNPs \
   --recode --stdout | gzip > ${resultsdir}/Subsampled_neutral.vcf.gz


################################################### ANNOTATING SNPS IN VCF FILE ############################################################



#/usr/lib/jvm/java-11-openjdk-11.0.22.0.7-2.el8.x86_64/bin/java -jar /media/inter/ssteindl/DEST/DESTv1_Pipeline/shell/snpEff/snpEff.jar
#
module load Tools/snpEff   

annotated="${resultsdir}/Subsampled_final_DP15.ann.vcf.gz"
java  -jar $snpEff ann BDGP6.28.99 $Sub4 | gzip >> $annotated
more $annotated | gunzip | awk ' !/^#/ {split($8,a,"|"); print $1 " " $2 " " a[4]}' > ${wd}/results/annotations.txt
#awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' ${wd}/results/${arm}/Subsampled_${arm}.final_DP15.af ${wd}/results/annotations.txt > ${wd}/results/annotated_used_frqs.txt



################################################### CONVERT TO ALLELE FREQUENCIES ############################################################
#echo "CONVERTING VCF FILE TO ALLELE FREQUENCY FILE"


python3 ${scriptdir}/VCF2AF.py --input $Sub4 > ${resultsdir}/Subsampled.final_DP15.af
python3 ${scriptdir}/VCF2AF.py --input $resultsdir/Subsampled_neutral.vcf.gz > ${resultsdir}/Neutral.final.af


#### RDA


AF_file="${resultsdir}/Subsampled_fullgenome2.final_DP15.af"
metadata="${wd}/dest_v2.samps_3May2024.csv"
neutral_af="${resultsdir}/Neutral.final.af"
RDA_out="${resultsdir}/RDA"
mkdir $RDA_out

