#!/bin/bash

dir="$1" # first argument
file="$2" #second argumment

NAME=`echo ${file%%_fastqc.zip}`
if [ ! -f /QC/logs_afterCutAdapt/$NAME ]; then
	ClusterStuff/FastQC/fastqc --outdir  $dir/QC_2ndTime -t 6 $file 
fi

mv Folder_name/QC_2ndTime/*.zip Folder_name/QC/logs_afterCutAdapt
mv Folder_name/QC_2ndTime/*.html Folder_name/QC/readout_afterCutAdapt
