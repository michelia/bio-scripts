#$ -S /bin/sh
#! /bin/bash

/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish count -m 17 -t 32 /scgene/tiger/invent/guoshuguang/kmerDate/1fq.fa -o /scgene/tiger/invent/guoshuguang/kmerDate/fa_jell_out/1fqout -s 10000000

/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish  merge -o /scgene/tiger/invent/guoshuguang/kmerDate/1fqout.jf /scgene/tiger/invent/guoshuguang/kmerDate/fa_jell_out/1fqout*

/scgene/elephant/pipeline/denovo/tools/jellyfish/jellyfish-1.1.6/bin/jellyfish histo /scgene/tiger/invent/guoshuguang/kmerDate/1fqout.jf -o  /scgene/tiger/invent/guoshuguang/kmerDate/1fq_histo.txt

python  /scgene/tiger/invent/guoshuguang/repo_michelia/kmer/plot_kmer.py /scgene/tiger/invent/guoshuguang/kmerDate/1fq_histo.txt /scgene/tiger/invent/guoshuguang/kmerDate/figure/1fq.png