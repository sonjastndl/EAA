# Environmental Association Analysis
## FAIRiCUBE - Use Case 3 🪰 

##  Table of Contents

<!--ts-->

1. [About the Repository](#about-the-repository)
2. [Hypothesis and Research Questions](#hypothesis-and-research-questions)
3. [Novel Aspects](#novel-aspects)
4. [Data](#data)
   * [Genetic Data](#genetic-data)
   * [Environmental Data](#environmental-data)
5. [Project Workflow](#project-workflow)
6. [Results](#results)

<!--ts-->

---

## About the Repository

This repository is supposed to hold documentation on the Environmental Association Analysis Workflow conducted in the scope of the UC3 Project of [FAIRiCUBE](https://fairicube.nilu.no/).
The Use Case approach aims to identify how environmental factors shape genetic variation and influence evolutionary processes in european popualtions of Drosophila melanogaster. By correlating population genomics with environmental data, the study aims to uncover genetic targets affected by environmental selection pressures.
This is a standalone repository linking to other repository of the project work concerning the following topics. 

- Accessing and Processing Environmental Data with FAIRiCUBE
- Conducting EAA on the FAIRiCUBE HUB


---

##  Hypothesis and Research Questions

Many research has been conducted on the model organism Drosophila melanogaster. This previous research and collaboration efforts led to the accumulation of genetic data of D. melanogaster from these worldwide populations. The [DEST](https://dest.bio/) (Drosophila Evolution over Space and Time) dataset is a valuable collection of genomic data from natural populations across multiple locations worldwide and multiple time points. 

---
- How does environmental variation across space and time correlate with patterns of genetic diversity in Drosophila melanogaster populations?
  
- Which genomic regions or genes in Drosophila melanogaster show signatures of adaptation to specific environmental conditions?
 
- Can combinations of environmental factors predict changes in genetic structure or the presence of adaptive alleles in natural populations?
<br>

Studies on the environmental adaptation of the organism revealed that climate specific  variables like [temperature, rainfall, and wind](https://pubmed.ncbi.nlm.nih.gov/33350518/) are associated with genetic and phenotypic variations and therefore significant factors driving local adaptation in D. melanogaster populations.

In this study, the researchers used BayPass, a genome-environment association tool that accounts for population structure when identifying genetic variants linked to environmental variables. BayPass models the covariance in allele frequencies across populations, allowing for robust detection of SNPs (single nucleotide polymorphisms) and TEs (transposable elements) associated with specific environmental factors. The analysis was performed on genomic data from 36 population samples of Drosophila melanogaster across Europe and North America, alongside 59 environmental variables, including temperature, rainfall, and wind.

We are now interested in whether incorporating more populations, richer environmental data, and applying multivariate statistical methods like redundancy analysis (RDA) might yield different or more refined insights. RDA allows the simultaneous modeling of multiple environmental predictors, potentially capturing complex adaptive responses missed by univariate approaches like BayPass. This comparison could help clarify whether previous findings remain consistent under broader data and methodological scopes.







## Novel Aspects


--- 

## Data

### Environmental Data


If you want to work with customly aquired environmental data, make sure you put it in the corresponding format, the number of variables (columns) is not limited.
The file structure can be taken from [Env.csv](data/environmental/fairicube/Env_imputed.csv) which is part of this repository and looks the following. 

| Sample ID       | Variable 1      | Variable 2      |
| --------------- | --------------- | --------------- |
| AT_Kar_See_1_2014-08-17 | 44 | 255 |
| CH_Vau_Vul_1_2020-09-15| -128| 240 |
|CM_Nor_Oku_1_2004-04-15|  63|  54|
|  |  |  |

The environmental data used in this project consists of the following:
- [Rasdaman OGC Web Service Data](https://fairicube.rasdaman.com/rasdaman/ows#/services): <br>
  Part of the FAIRiCUBE infrastructure, we considered more than 150 layers including pesticide layers for our research.
  Ingestion of custom datsets during the project, a tool on how to accquire this data can be found below. 
- [WorldClim Data](https://www.worldclim.org/): <br>
  WorldClim offers 19 BioVariables derived from monthly temperature and rainfall values, generating biologically meaningful variables.
  These are often used for species distribution modeling and ecological studies.
- [Copernicus Data ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview):<br>
  ERA5 is a climate reanalysis product from the Copernicus Climate Change Service (C3S),
  providing hourly estimates of various atmospheric, land, and ocean variables.
  It's a globally available dataset, updated monthly, with a delay of approximately two months. This data can be downloaded by using the Climate Data Store (CDS) [API](https://cds.climate.copernicus.eu/how-to-api). 
- [Metdata file](EAA/data/environmental/fairicube/Env_all_meta.csv):<br>
  In this file you can see the datasets included and supporting information on source, value interpretation and Null Value encoding which is important for data cleanup strategies. 



#### Suggested Download Tool FAIRiCUBE QueryCube

- [QueryCube GitHub Repository](https://github.com/FAIRiCUBE/uc3-drosophola-genetics/tree/main/projects/QueryCube)

- [QueryCube Web Application](https://querycube.nilu.no/)

- [Code to Access and Download CDS](scripts/getCDSdata.py): <br>
  Please make sure to register at the CDS Website and agree to terms of condition of ERA5 dataset when downloading. 


### Genomic Data

This project is working with population samples from the model organism Drosophila melanogaster.

#### Download Genomic Data And Metadata

The download is part of the [EAA script](shell/EAA.sh) of the main pipeline. The data on geolocation of Drosohila melanogaster samples togehter with the according metadata, as well as the genetic data as VCF can be retrieved via DEST.bio.
Please be aware that these specific VCFs are large files and therefore download might take longer as usual. 


```bash
wget --tries=inf "http://berglandlab.uvadcos.io/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz"
wget "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv"
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz


## Extract the information on the available samples
awk '{FS=","}{if (NR!=1) {print $1}}' dest_v2.samps_3May2024.csv > ${data}/samplenames.csv
cp dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz ${data}/PoolSeq2024.vcf.gz
mv  dmel-all-r6.57.gff.gz > ${data}/dmel-all-r6.57.gff.gz

```

#### Additional Genomic Data

Later in the workflow, estimates for population structure and gene annotations will be used. We used FlyBase and snpEff to get intronic SNPs to inferr population structure and to annotate the VCF file. 

```
# Get Intronic SNPs

wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz

``` 
Move the GFF file.

```
# Create VCF of intronic SNPs only
vcf_file="${data}/dmel-all-r6.57.gff.gz"
gff_file="${wd}/data/dmel-all-r6.57.gff.gz"
neutralSNPs="${wd}/results/${arm}/Subsampled_NeutralSNPS_80.tsv"

python3 ${scriptdir}/IntronicSNPS.py \
 --gff $gff_file \
 --vcf $vcf_file \
 --target-length 80 \
 --output $neutralSNPs

vcftools --gzvcf $vcf_file \
   --positions $neutralSNPs \
   --recode --stdout | gzip > ${wd}/results/Subsampled_neutral.vcf.gz

```


Annotating SNPs with [snpEFF](http://pcingola.github.io/SnpEff/)

```
snpEff="/opt/bioinformatics/snpEff/snpEff.jar"

module load Tools/snpEff        

annotated="${wd}/results/Subsampled_ann.vcf.gz"
java -jar $snpEff ann BDGP6.28.99 $Sub4 | gzip >> $annotated
more $annotated | gunzip | awk ' !/^#/ {split($8,a,"|"); print $1 " " $2 " " a[4]}' > ${wd}/results/annotations.txt

```

--- 

## Project Workflow

Here are some of the key features of this project:

- **Aquiring Data** – See in the chapter [Data](#data) above.
  
- **Cleaning Data** – Having clean data is essential for reliable results.
  
> **Clean Up of Genetic Data included:** 
>>  -  Removing Low Quality Sites
>>  -  Removing incomplete sites (Meaning the genomic site information is not present in all populations)
>>  -  Removing monomorphic sites (The give no additional information)
>>  -  Removing extreme values (e.g. sites with extremely low/high avarage allele frequencies: 0,05 < AF < 0,95)
>>  -  Imutation of missing data (This is not performed here.
>>      However in the course of the project, a method for imputing population genomic data at allele frequency level was developed: [see here](https://github.com/FAIRiCUBE/uc3-drosophola-genetics/tree/main/projects/gap_filling).) 
>
> **Clean Up of Environmental Data involves:**
>>  -  Analysis of Completeness of Data 
>>  -  Removing non-unique variables 
>>  -  Removing highly intercorrelated variables 
>>  -  Removing variables with missing data above threshold 
>>  -  Impute missing data (We imputed categorial data by using the most common value, numerical data was imputed with inverse distance weighting in our [Rscript]().)
>
>
>
> When cleaning the data, we tested two approaches with different thresholds. You can see the entire CleanUp Strategy in our [RMD file](results/processing/cleaning/Filter_Clean_Impute_Data.Rmd) and the results in the [CleanUp Folder](results/processing/cleaning).


> | Filter 1  | Filter 1 Threshold  | Filter 2  | Filter 2 Threshold | Populations | Env | WorldClim Env | Env Total | Total Imputed|
> | --------------- | --------------- | --------------- |--------------- | --------------- | --------------- |--------------- | --------------- |  --------------- | 
> | **Sample** | **15 %** | **Env** | **15%** | **293** | **135** | **19**| **154** | **0.01794903**| 
> |  Env| 10% | Sample | 10% | 28 | 126| 19| 145 |
> | **Env** | **15%** | **Sample** | **15%** | **179** | **141**| **19**| **160**| **0.00454633**| 
> | Sample | 10% | Env | 10% | 276 | 86 | 19 | 105 |
 

![Strategy A](results/processing/cleaning/Pops_StrategyA.png)
![Strategy B](results/processing/cleaning/Pops_StrategyB.png)


- **Redundancy Analysis** – This is an R script applicable across various platforms and operating systems.
   - **Intersecting Data**
   - **Additional Data** -
   - **Variable selection** ordiR2step
   - **Variance partitioning** - 
   - **Permutations** 
   - **Threshold**
   - **GOterm analysis** 
---

## Results


### Test First PLot
![Correlation of Environemtnal Variables](results/RDA/partialRDA/plot/CorrelationEnv.png)

## Resources

Make sure you have the following installed:

- [Python](https://www.python.org/)
- [VCFTools](http://vcftools.github.io/license.html)
- [vegan R-Package](https://www.bioconductor.org/packages/release/bioc/html/LEA.html)


### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/your-username/project-name.git
   cd project-name
