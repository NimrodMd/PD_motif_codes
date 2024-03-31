#!/bin/bash

dir="$1" # first argument
file="$2" #second argumment

NAME=`echo ${file%%_fastqc.zip}`
if [ ! -f /QC/$NAME ]; then
	ClusterStuff/FastQC/fastqc --outdir  $dir/QC -t 6 $file 
fi

mv Folder_name/QC/*.zip Folder_name/QC/logs_Raw
mv Folder_name/QC/*.html Folder_name/QC/readout_Raw
