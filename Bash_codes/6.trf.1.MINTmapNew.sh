#!/bin/bash

dir="$1"  #first argument
file="$2" #second argumment

NAME1=$(echo $file| cut -d'.' -f2)
NAME2=`echo $NAME1 | cut -d '/' -f5`
NAME=`echo Folder_name/MINTmap_output/$NAME2`
echo $NAME2
		
/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/MINTmap_files/MINTmap.pl -f $file -p $NAME -l /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/MINTmap_files/LookupTable.tRFs.MINTmap_v1.txt -s /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/MINTmap_files/tRNAspace.Spliced.Sequences.MINTmap_v1.fa -o /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/MINTmap_files/OtherAnnotations.MINTmap_v1.txt -j /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/MINTmap_files/MINTplates
