def TopoPRED(w_size, x_value):
    
    print('\nTopoPRED is now running with a window size of ' + str(w_size) + ' and ' + str(x_value) + '-fold cross-validation.\nResults will be sent to folder /output.\n')

    import time
    import math
    import itertools
    import numpy as np
    import statistics as stat
    from sklearn import svm, datasets

    clockstarts=time.time()

    #####UTILS####
    ws=w_size
    x=int((ws-1)/2)
    p=[list(np.zeros(20))]
    zero=0.0
    d={1:'M', 2:'O', 3:'I'}
    c_value=1
    
    #Iterates through file F, separating protein IDs from topologies
    #Takes window size into account to create subsets of sequences:topologies pairs
    #Integrates data from PSSMatrices obtained after PSI-BLAST at evalue 0.005

    with open('../input/query.fa', 'r') as F:
        f=F.read().splitlines()

    IDs=[]
    seqs=[]
    
    for i in range(0, len(f), 3):
        IDs.append(f[i].strip('>'))
    for i in range(1, len(f), 3):
        seqs.append(f[i])

    #Defines function to normalize PSSMatrices values
    def sigmoid(x):
        return (1 / (1 + (math.exp(-x))))

    test_PSSM = []
    seq_test = []
    for i in range(0, len(IDs)):
        for j in range(0, len(IDs[i])):
            #seqs.append(ZIPPED[i][1])
            with open('../input/' + str(IDs[i]) + '.fa.pssm', 'r') as F:
                f=F.read().splitlines()
                del f[0:3]
                del f[-6:]
                SMatx=[]
                temp=[]
                split=[]
                for k in f:
                    split=k.split()
                    temp=split[2:22]
                    temp2=[]
                    for l in temp:
                        temp2.append(sigmoid(int(l)))
                    SMatx.append(temp2)
                test_PSSM.append(SMatx)

    X_test = []
    for protein in test_PSSM:
        for i in range(0, len(protein)):
            if i < x:
                SM_w=(p*(x-i)) + protein[0:ws-(x-i)]
                SM_w=[j for i in SM_w for j in i]
                X_test.append(SM_w)
            elif i >= (len(protein)-x):
                SM_p=protein[(i-x):ws-(x-i)+1]
                if len(SM_p) != ws:
                    SM_p.extend(p*(ws-len(SM_p)))
                SM_p=[j for i in SM_p for j in i]
                X_test.append(SM_p)
            else:
                SM_w=protein[(i-x):(i+x+1)]
                SM_w=[j for i in SM_w for j in i]
                X_test.append(SM_w)
        ###############################
        #CLASSIFIER_HERE
        ###############################
        CLF = joblib.load('../output/SM_SVM_PSSMw15.pkl')

        RESULTS=[]
    
        PREDICTION=CLF.predict(X_test)
        RESULTS = list(PREDICTION)
        RESULTS = [d[i] for i in RESULTS]
        RESULTS = ''.join(RESULTS)
        pred = []
        pred.append(RESULTS)
        
    with open('../output/TopoPRED_results.txt', 'w') as r:
       for i in range(0, len(IDs)) :
           r.write('>' + str(IDs[i]) + '\n')
           r.write(str(seqs[i]) + '\n')
           r.write(str(pred[i]) + '\n')

#TopoPRED(w_size=3, x_value=3)
#TopoPRED(w_size=7, x_value=3)
#TopoPRED(w_size=11, x_value=3)
TopoPRED(w_size=15, x_value=3)
#TopoPRED(w_size=19, x_value=3)
#TopoPRED(w_size=23, x_value=3)
#TopoPRED(w_size=27, x_value=3)
#TopoPRED(w_size=31, x_value=3)
