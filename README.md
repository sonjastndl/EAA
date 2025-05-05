# Environemntal Association Analysis
## FAIRiCUBE - USe Case 3 🪰 

##  Table of Contents

<!--ts-->

1. [About the Repository](#about-the-repository)
2. [Hypothesis](#hypothesis)
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
This is a standalone repository linking to other repository of the project work concerning the following topics. 

- Accessing and Processing Environmental Data with FAIRiCUBE
- Conducting EAA on the FAIRiCUBE HUB


---

##  Hypothesis

## Novel Aspects

## Data

### Genetic Data

This project is working with population samples from the model organism Drosophila melanogaster.

#### Download Genomic Data

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


### Environmental Data


If you want to work with customly aquired environmental data, make sure you put it in the corresponding format, the number of variables (columns) is not limited. The file structure can be taken from [Env.csv](data/environmental/fairicube/Env_imputed.csv) which will be created in the next step of the workflow. 

| Sample ID       | Variable 1      | Variable 2      |
| --------------- | --------------- | --------------- |
| AT_Kar_See_1_2014-08-17 | 44 | 255 |
| CH_Vau_Vul_1_2020-09-15| -128| 240 |
|CM_Nor_Oku_1_2004-04-15|  63|  54|
|  |  |  |

#### Suggested Download Tool FAIRiCUBE QueryCube

[QueryCube GitHub Repository](https://github.com/FAIRiCUBE/uc3-drosophola-genetics/tree/main/projects/QueryCube)
[QueryCube Web Application](https://querycube.nilu.no/)



## Project Workflow

Here are some of the key features of this project:

- **Aquiring Data** – Simple interface and clean design.
- **Cleaning Data** – Easy to extend or integrate into other systems.
- **Intersecting Data** – Includes unit and integration tests for reliability.
- **Redundancy Analysis** – Works across different environments and operating systems.

---

## Results



### Resources

Make sure you have the following installed:

- [Python](https://www.python.org/)
- [VCFTools](http://vcftools.github.io/license.html)
- [vegan R-Package](https://www.bioconductor.org/packages/release/bioc/html/LEA.html)


### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/your-username/project-name.git
   cd project-name

