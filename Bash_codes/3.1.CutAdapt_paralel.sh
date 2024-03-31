#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2-0:0:0
#SBATCH -o Error-%j.err

dir="/ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name"
export LD_LIBRARY_PATH=/ClusterStuff/flexbar1/flexbar-3.5.0-linux:ClusterStuff/flexbar1/flexbar-3.5.0-linux:$LD_LIBRARY_PATH

mkdir $dir/flexbar_q

for f in $dir/FASTQ/*.gz
do
	echo $f
	sbatch -c1 --time=2-0:0:0 /ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name/3.2.CutAdapt_in.sh	$dir $f
done &> FastQC_$(date +'%m.%d.%Y').log

