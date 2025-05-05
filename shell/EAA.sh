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
   --recode --stdout | gzip > ${results}/Subsampled_neutral.vcf.gz

################################################### CONVERT TO ALLELE FREQUENCIES ############################################################
#echo "CONVERTING VCF FILE TO ALLELE FREQUENCY FILE"


python3 ${scriptdir}/VCF2AF.py --input $Sub4 > ${resultsdir}/Subsampled.final_DP15.af
python3 ${scriptdir}/VCF2AF.py --input $resultsdir/Subsampled_neutral.vcf.gz > ${resultsdir}/Neutral.final.af


################################################## PERFORM LINEAR REGRESSION #################################################################
echo "PERFORMING LINEAR REGRESSION"

Rscript ${scriptdir}/Plot_pvalues.R $wd $AF $envdata $arm $summary

 
################################################## PERFORM LFMM 2 (LATENT FACTOR MIXED MODEL) ###################################################
echo "PERFORMING LFMM"

LeaOut="${resultsdir}/LEA"
mkdir $LeaOut
#
variables=$(head -n 1 "$envdata" | sed 's/\r$//')
#echo $variables
##  #Use a for loop to iterate over words (assuming space-separated words)
IFS=',' read -ra header_elements <<< "$variables"

#Rscript /home/sonjastndl/s3/InstallLea.R
##test
#Rscript ${scriptdir}/LEA_RunLFMM2.R $LeaOut $AF $metadata_new "bio1" $rep

#for ((i = 1; i < 3; i++)); do
for ((i = 1; i < ${#header_elements[@]}; i++)); do
    element=${header_elements[i]} 
    echo "$element"
    #echo $LeaOut
    #echo $AF
    rep=1
    echo $rep
    # Add your processing here
    Rscript ${scriptdir}/LEA_RunLFMM2.R $LeaOut $AF $envdata $element $rep
    #Rscript ${scriptdir}/LEA_RunLFMM2.R $LeaOut $AF $metadata_new "bio1" 1
    #If needed average the Repetitions
    #Rscript Rscript ${scriptdir}/LEA_ZPcalc.R $LeaOut $nK $nR $AF $var
done

Rscript ${scriptdir}/PlotLEAPValues.r $wd $AF $metadata_new $arm $summary

Rscript ${scriptdir}/ComparePValues.R $AF ${resultsdir}/GM $LeaOut $summary

#### RDA
## Get Intronic SNPs
#
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz
#

AF_file="${resultsdir}/Subsampled_fullgenome2.final_DP15.af"
metadata="${wd}/dest_v2.samps_3May2024.csv"
neutral_af="${resultsdir}/Neutral.final.af"
RDA_out="${resultsdir}/RDA"
mkdir $RDA_out

