#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2-0:0:0
#SBATCH -o Error-%j.err

dir="/ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name"
mkdir $dir/QC_2ndTime

for f in $dir/flexbar_q/*.fastq
do
	echo $f
	sbatch -c1 --time=2-0:0:0 /ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name/4.2.fastQC_afterCutAdapt_in.sh	$dir $f
done &> FastQC_$(date +'%m.%d.%Y').log

mkdir Folder_name/QC/readout_afterCutAdapt
mkdir Folder_name/QC/logs_afterCutAdapt
