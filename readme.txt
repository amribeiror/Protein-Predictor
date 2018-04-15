Welcome to TopoPRED, a membrane protein topology predictor.

The main script for the predictor (TopoPRED.sh) contains the main commands.

Once run from the command line, TopoPRED.sh (1) generates the required folder structure, (2) runs a PSIBLAST to extract positional information for each amino acid residue in the query protein sequence(s) and (3) runs a linear SVM classifier using the Substitution Matrix information with a window size = 15; these are the parameters that yield more accurate results with a test dataset (confusion matrices can be found in the folder 'graphs').

For TopoPRED.sh to work, the query Fasta file **AND** the model SM_SVM_PSSMw15.pkl **AND** the python script SM_SVM_PSSMpredictor.py must be in the folder 'input'.

Additional scripts used to generate models for Random Forests (SM_RF_PSSM.py), Decision Trees (SM_DT_PSSM.py) and SVM using FM matrices (FM_SVM_PSSM.py) can be found in the folder 'scripts', but are not required for TopoPRED to work.
