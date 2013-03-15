#$ -S /bin/sh
#! /bin/bash

# gzip -dc /scgene/elephant/rawdata/20110808LYGRICEDNA/fastq/Clean_data/110718_I270_FCC04BHABXX_L8_RIChjjRAADIAAPEI-2_1.fq.gz /scgene/elephant/rawdata/20110808LYGRICEDNA/fastq/Clean_data/110718_I270_FCC04BHABXX_L8_RIChjjRAADIAAPEI-2_2.fq.gz > /scgene/tiger/invent/guoshuguang/kmerDate/raw_all.fq

/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish count -m 27 -o /scgene/tiger/invent/guoshuguang/kmerDate/estimate_genome_size/count_27_out -C -s 900000000 -U 500 -t 32 /scgene/tiger/invent/guoshuguang/kmerDate/raw_all.fq



/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish histo /scgene/tiger/invent/guoshuguang/kmerDate/estimate_genome_size/count_27_out_0 > /scgene/tiger/invent/guoshuguang/kmerDate/estimate_genome_size/count27.histo

python /scgene/tiger/invent/guoshuguang/repo_michelia/kmer/plot_kmer.py /scgene/tiger/invent/guoshuguang/kmerDate/estimate_genome_size/count27.histo /scgene/tiger/invent/guoshuguang/kmerDate/figure/estimate_genome_27count.png

# /scgene/tiger/invent/guoshuguang/estimate_genome_size.pl-master/estimate_genome_size.pl --kmer=17 --peak=24 --fastq=/scgene/elephant/rawdata/20110808LYGRICEDNA/fastq/Clean_data/110718_I270_FCC04BHABXX_L8_RIChjjRAADIAAPEI-2_1.fq.gz /scgene/elephant/rawdata/20110808LYGRICEDNA/fastq/Clean_data/110718_I270_FCC04BHABXX_L8_RIChjjRAADIAAPEI-2_2.fq.gz

