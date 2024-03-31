#!/bin/sh

dir="$1"  #first argument
read="$2" #second argumment

 
echo $read
NAME=$(echo $read| sed 's/.fastq//')
out=`echo $NAME | cut -d '/' -f9` 
mkdir $dir/miRDeep2/$out
echo $out
cd $dir/miRDeep2/$out
        
/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/bin/mapper.pl $read -e -h -i -j -l 16 -m -p /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/Homo_sapiens.GRCh38.dna.primary_assembly_107 -s $dir/miRDeep2/$out/reads_collapsed.fa -t $dir/miRDeep2/$out/reads_collapsed_vs_genome.arf -v
         
/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/bin/quantifier.pl -p /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/hsa_hairpin.fa -m /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/hsa_mature.fa -g 0 -r $dir/miRDeep2/$out/reads_collapsed.fa -t hsa -y 16_19 -d
#-d           if parameter given pdfs will not be generated, otherwise pdfs will be generated 

/ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/mirdeep2-master/bin/miRDeep2.pl $dir/miRDeep2/$out/reads_collapsed.fa /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/Homo_sapiens.GRCh38.dna.primary_assembly_107.fa $dir/miRDeep2/$out/reads_collapsed_vs_genome.arf /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/hsa_mature.fa /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/mature.fa /ems/elsc-labs/soreq-m/nimrod.madrer/ClusterStuff/genomes/hsa_hairpin.fa -t H.sapiens 2> report.log
