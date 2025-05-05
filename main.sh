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
