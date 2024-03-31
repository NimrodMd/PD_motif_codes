#!/bin/bash

dir="$1" # first argument
file="$2" #second argumment

echo $file
tmp=$(echo $file| cut -d'.' -f 2)
echo $tmp	
SUBSTRING=$(echo $tmp| cut -d'/' -f 5)
echo $SUBSTRING
out=$(echo $SUBSTRING| sed 's/.fastq//') 
if [ ! -e $dir/flexbar_q/$out.fastq ]; then
	/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/flexbar1/flexbar-3.5.0-linux/flexbar -r $file -a /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/Adaptors/illumina_universal.fasta -q TAIL -qf sanger -qw 4 -min-read-length 16 -n 1 -t /ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name/flexbar_q/$SUBSTRING
else
	echo "File exists."
fi