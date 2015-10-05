#$ -S /bin/sh
#! /bin/bash

# head -n 127074388 /scgene/tiger/invent/guoshuguang/kmerDate/raw_all.fq > /scgene/tiger/invent/guoshuguang/kmerDate/1_4.fq

# /scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish count -m 17 -o /scgene/tiger/invent/guoshuguang/kmerDate/same_kmer_out/17kmer1_4 -C -s 900000000 -U 500 -t 32 /scgene/tiger/invent/guoshuguang/kmerDate/1_4.fq



/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish histo /scgene/tiger/invent/guoshuguang/kmerDate/same_kmer_out/17kmer1_4_0 > /scgene/tiger/invent/guoshuguang/kmerDate/same_kmer_out/1_4.histo

python /scgene/tiger/invent/guoshuguang/repo_michelia/kmer/plot_kmer.py /scgene/tiger/invent/guoshuguang/kmerDate/same_kmer_out/1_4.histo /scgene/tiger/invent/guoshuguang/kmerDate/same_kmer_out/1_4.png

# /scgene/tiger/invent/guoshuguang/estimate_genome_size.pl-master/estimate_genome_size.pl --kmer=17 --peak=24 --fastq=/scgene/elephant/rawdata/20110808LYGRICEDNA/fastq/Clean_data/110718_I270_FCC04BHABXX_L8_RIChjjRAADIAAPEI-2_1.fq.gz /scgene/elephant/rawdata/20110808LYGRICEDNA/fastq/Clean_data/110718_I270_FCC04BHABXX_L8_RIChjjRAADIAAPEI-2_2.fq.gz
