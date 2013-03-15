#$ -S /bin/sh
#! /bin/bash


#/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish count -m 17 -t 32 /scgene/elephant/data/train/fq/rice_genome_1.fq -o jell_kmer_out/out_ -s 10000000
/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish  merge -o all_out.jf jell_kmer_out/out__*
