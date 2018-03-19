#!/bin/bash
echo ""
echo "Ahoy!\nWelcome to TopoPRED, a membrane protein predictor."
echo ""
echo "Please make sure that your test fasta file is in the folder /data with the following name format: <filename>.fa."
echo ""
echo "TopoPRED will now prepare your test file."
echo ""

#python data_prep.py

cd ../bin/

#if [ ! -f .fasta.gz ]; then
#    wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
#gunzip *.gz
#fi

#if [ ! -f *.phr ]; then
 #   makeblastdb -in uniref50.fasta -dbtype "prot" -out nr
    #makeblastdb -in uniref90 BLA BLA BLA
#fi

for file in *.fa; do
    if [ ! -f $file.pssm ]; then
#	echo "Running PSI-BLAST on $file at $(date). Please wait ..."
#	echo ""
	psiblast -query $file -evalue 10 -db /mnt/c/KB8024/bin/nr -num_iterations 2 -out ../data/$file.psiblast -out_ascii_pssm ../data/$file.pssm
#	echo ""
#	echo "PSI-BLAST on $file terminated at $(date)."
    fi
done

#Building a new DB, current time: 03/17/2018 13:13:34
#New DB name:   /mnt/c/KB8024/bin/nr
#New DB title:  uniref50.fasta
#Sequence type: Protein
