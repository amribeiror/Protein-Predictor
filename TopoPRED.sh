#!/bin/bash
echo ""
echo "Ahoy!\nWelcome to TopoPRED, a membrane protein predictor."
echo ""
echo "Please make sure that your test fasta file is in the folder /data with the following name format: <filename>.fa."
echo ""
echo "TopoPRED will now prepare your test file."
echo ""

python data_prep.py

cd ../bin/

if [ ! -f .fasta.gz ]; then
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip *.gz
fi

if [ ! -f *.phr ]; then
    makeblastdb -in uniref90.fasta -dbtype "prot" -out nr
fi

for file in *.fa; do
    if [ ! -f $file.pssm ]; then
	echo "Running PSI-BLAST on $file at $(date). Please wait ..."
	echo ""
	psiblast -query $file -evalue 0.005 -db /mnt/c/KB8024/bin/nr -num_iterations 3 -out ../data/$file.psiblast -out_ascii_pssm ../data/$file.pssm
	echo ""
	echo "PSI-BLAST on $file terminated at $(date)."
    fi
done

python SM_SVM_PSSM.py
