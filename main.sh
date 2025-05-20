## Main Environmental Association Analysis pipeline can be run here

wd="/media/inter/ssteindl/EAA"
continent="Europe"
input=""
metadata=""
envdata=""
AF=""
#scriptdir=""


bash  /media/inter/ssteindl/EAA/shell/EAA.sh $wd \
    $continent \
    $samplelist \
    $input \
    $metadata \
    $envdata \
    $AF \
    $scriptdir >> /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/logs/LandscapeGenomics160125_qsub.log 2>&1


## 100 permutations of the model are stored here

permuted_values=$(seq 3 101 | paste -sd, -)


python3 scripts/FDR-rank.py \
    --input "/media/inter/ssteindl/EAA_NHM/EAA/results/RDA/permutations/Pvalues_merged_with_permutations.csv" \
    --true -2 \
    --permuted -3,4,5 \
    --output "results/RDA/permutations/Pvalues_FDR2.csv"

