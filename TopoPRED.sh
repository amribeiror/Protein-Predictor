#!/bin/bash
echo ""
echo "Ahoy!\nWelcome to TopoPRED, a membrane protein predictor."

mkdir -p TopoPRED/{scripts,input,output,logs,bin}

echo ""
echo "Please make sure that your query fasta file is in the folder TopoPRED/input with the following name format: <filename>.fa."
echo ""

mv project_andre_rosa.py PROJECTS/project_one/scripts
echo "TopoPRED will now download uniref90 for a PSIBLAST run on the query proteins."
echo ""

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
	psiblast -query $file -evalue 0.005 -db nr -num_iterations 3 -out ../input/$file.psiblast -out_ascii_pssm ../input/$file.pssm
	echo ""
	echo "PSI-BLAST on $file terminated at $(date)."
    fi
done

echo "TopoPRED will now run the model SM_SVM_PSSMw15.pkl on the query proteins."

python SM_SVM_PSSMw15.py
