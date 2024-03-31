#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2-0:0:0
#SBATCH -o Error-%j.err

dir="/ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name"

rm $dir/QC_2ndTime

export LD_LIBRARY_PATH=/ClusterStuff/flexbar1/flexbar-3.5.0-linux:ClusterStuff/flexbar1/flexbar-3.5.0-linux:$LD_LIBRARY_PATH
export PATH=$PATH:/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/src >> ~/.bashrc
export PATH=$PATH:/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/bin >> ~/.bashrc
export PATH=$PATH:/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/essentials/bowtie-1.1.1 >> ~/.bashrc
export PATH=$PATH:/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/essentials/ViennaRNA-1.8.4/install_dir/bin >> ~/.bashrc
export PATH=$PATH:/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/essentials/randfold-2.0 >> ~/.bashrc
export PERL5LIB=$PERL5LIB:/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/lib/perl5 >> ~/.bashrc


source ~/.bashrc

mkdir $dir/MINTmap_output
mkdir $dir/miRDeep2

for f in $dir/flexbar_q/*.fastq
do
	echo $f
	sbatch -c1 --time=2-0:0:0 /ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name/6.trf.1.MINTmapNew.sh	$dir $f
	sbatch -c1 --time=2-0:0:0 /ems/elsc-labs/soreq-m/nimrod.madrer/Folder_name/6.miR.1.miRDeep_loop.sh	$dir $f
done &> TRF_$(date +'%m.%d.%Y').log
