#!/bin/sh


##'DEFAULT VALUES IF NOT CHANGED'
DIRECTORY=$(pwd)
CONTINENT="Europe"
SAMPLES="All"
VCFFILE="Default Download"
METADATA="Default Download"
ENVDATA="Default Download"
OUTPUTDIR=$(pwd)

scriptdir=$DIRECTORY



while getopts d:o:s:i:c:? opt 2>/dev/null
do
  case $opt in
    d) DIRECTORY=$OPTARG;;
    o) OUTPUTDIR=$OPTARG;;
    s) SAMPLES=$OPTARG;;
    i) VCFFILE=$OPTARG;;
    c) CONTINENT=$OPTARG;;
    ?) echo "Valid parameters are: [-d, -s]" ;;
  esac
done

if [ -n "$DIRECTORY" ]; then
  echo "Working Directory is set to: $DIRECTORY"
fi

if [ -n "$OUTPUTDIR" ]; then
  echo "Results will be written to: $OUTPUTDIR"
  read -p " Do you want to proceed? [Y/y|N/n]" -n 1 -r
  if [[ ! $REPLY =~ [Yy]$ ]];
    then
        echo "Programm stopped"
        exit 1
    fi
fi

if [ -n "$VCFFILE" ]; then
  echo "VCFFILE is set to: $VCFFILE"
fi


if [ -n "$METADATA" ]; then
  echo "METADATA is set to: $METADATA"
fi

if [ -n "$SAMPLES" ]; then
  echo "SAMPLES is set to: $SAMPLES"
fi



scriptdir=$(pwd)/scripts
wd=$DIRECTORY
continent=$SAMPLES
#samplelist="$3"
#input="$4"
#metadata="$5"
#envdata="$6"
#