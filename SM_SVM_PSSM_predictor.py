def TopoPRED(w_size, x_value):
    
    print('\nTopoPRED is now running with a window size of ' + str(w_size) + ' and ' + str(x_value) + '-fold cross-validation.\nRBF was chosen as classifier.\nResults will be sent to folder /output.\n')

    import time
    import math
    import itertools
    import numpy as np
    import statistics as stat
    from sklearn import svm, datasets
    from sklearn.utils import shuffle
    from sklearn.externals import joblib
    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import f1_score
    from sklearn.metrics import precision_score
    from sklearn.metrics import average_precision_score
    from sklearn.metrics import recall_score
    from sklearn.metrics import roc_auc_score
    from sklearn.metrics import classification_report
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

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

    with open('../input/tm_alpha_beta_3state.3line_pssm.txt', 'r') as F:
        f=F.read().splitlines()

    IDs=[]
    seqs=[]
    
    for i in range(0, len(f), 3):
        IDs.append(f[i].strip('>'))
    for i in range(1, len(f), 3):
        seqs.append(f[i])
#    print(IDs)
#    random_IDs=[]
#    random_topols=[]

   # ZIPPED=list(zip(IDs,seqs))
#    SHUFFLED=shuffle(ZIPPED)
    
#    for i in range(0, len(SHUFFLED)):
#        random_IDs.append(SHUFFLED[i][0])
#        random_topols.append(SHUFFLED[i][1])

    #Creates cross validated datasets
#    avg=len(random_IDs)/float(x_value)
#    IDs_train=[]
#    IDs_test=[]
#    topols_train=[]
#    topols_test=[]
    
#    while zero < len(random_IDs):
#        IDs_train.append(random_IDs[0:int(zero)] + random_IDs[int(zero + avg):])
#        IDs_test.append(random_IDs[int(zero):int(zero + avg)])
#        topols_train.append(random_topols[0:int(zero)] + random_topols[int(zero + avg):])
#        topols_test.append(random_topols[int(zero):int(zero + avg)])
#        zero+=avg

    #Defines function to normalize PSSMatrices values
    def sigmoid(x):
        return (1 / (1 + (math.exp(-x))))

    #Prepares train,test datasets to sklearn
#    train_PSSM = []
#    topol_train = []
#    for i in range(0, len(IDs_train)):
#        for j in range (0, len(IDs_train[i])):
#            topol_train.append(SHUFFLED[i][1])
#            with open('../input/' + str(SHUFFLED[i][0]) + '.fa.pssm', 'r') as F:
#                f=F.read().splitlines()
#                del f[0:3]
#                del f[-6:]
#                SMatx=[]
#                temp=[]
#                split=[]
#                for k in f:
#                    split=k.split()
#                    temp=split[2:22]
#                    temp2=[]
#                    for l in temp:
#                        temp2.append(sigmoid(int(l)))
#                    SMatx.append(temp2)
#                train_PSSM.append(SMatx)

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

#    X_train = []
#    for protein in train_PSSM:
#        for i in range(0, len(protein)):
#            if i < x:
#                SM_w=(p*(x-i)) + protein[0:ws-(x-i)]
#                SM_w=[j for i in SM_w for j in i]
#                X_train.append(SM_w)
#            elif i >= (len(protein)-x):
#                SM_p=protein[(i-x):ws-(x-i)+1]
#                if len(SM_p) != ws:
#                    SM_p.extend(p*(ws-len(SM_p)))
#                SM_p=[j for i in SM_p for j in i]
#                X_train.append(SM_p)
#            else:
#                SM_w=protein[(i-x):(i+x+1)]
#                SM_w=[j for i in SM_w for j in i]
#                X_train.append(SM_w)

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
#     print(pred)
#    accuracy=[]

#    c_matrix=confusion_matrix(y_test,PREDICTION,labels=[1,2,3])
#    c_matrix_n=c_matrix.astype('float')/c_matrix.sum(axis=1)[:, np.newaxis]
#    f1_score=f1_score(y_test,PREDICTION,average='weighted')
#    p_score=precision_score(y_test,PREDICTION,average='weighted')
#    r_score=recall_score(y_test,PREDICTION,average='weighted')
#    report=classification_report(y_test,PREDICTION,target_names=classes)
#    accuracy.append(accuracy_score(y_test,PREDICTION))	
#    avg_accuracy=sum(accuracy)/len(accuracy)
    
    #Outputs to file at folder /Output
#    printfile=open('../output/pred_SM_SVM_ws_' + str(w_size) + '_x-fold-validation_' + str(x_value) + '.txt', 'w')
#    printfile.write('Average accuracy is: ' + str(avg_accuracy) + '\n')
#    printfile.write('Accuracies are: ' + str(accuracy) + '\n')
#    printfile.write('F1 score is: ' + str(f1_score) + '\n')
#    printfile.write('Precision score is: ' + str(p_score) + '\n')
#    printfile.write('Recall score is: ' + str(r_score) + '\n')
#    printfile.write(report + '\n')
        
#    from sklearn.model_selection import train_test_split
       
#    classes=[1,2,3]
#    plt.figure()
#    plt.imshow(c_matrix_n,interpolation='nearest',cmap='PuBuGn')
#    plt.title('SM_SVM_ws: ' + str(w_size) + ' fold-validation: ' + str(x_value) + ' F1: ' + str(f1_score))
#    plt.colorbar()
#    tick_marks = np.arange(len(classes))
#    plt.xticks(tick_marks,classes,rotation=0)
#    plt.yticks(tick_marks,classes)
    
#    threshold=c_matrix_n.max()/2.
#    for i, j in itertools.product(range(c_matrix_n.shape[0]), range(c_matrix_n.shape[1])):
#        plt.text(j,i,round(c_matrix_n[i,j],3),horizontalalignment="center",color="white" if c_matrix_n[i, j] > threshold else "black")
#    plt.tight_layout()
#    plt.ylabel('True label')
#    plt.xlabel('Predicted label')
#    plt.savefig('../output/pred_SM_SVM_ws_' + str(w_size) + 'fold_validation_' + str(x_value) + '.png')

#    clockends=time.time() - clockstarts
    #for i in range(0, len(IDs)) :
    with open('../output/TopoPRED_results.txt', 'w') as r:
       for i in range(0, len(IDs)) :
           r.write('>' + str(IDs[i]) + '\n')
           r.write(str(seqs[i]) + '\n')
           r.write(str(pred[i]) + '\n')

#y    print('Prediction completed in %0.2f seconds'.%(clockends))
#    joblib.dump(CLF, '../output/SM_SVM_PSSM.pkl')

#TopoPRED(w_size=3, x_value=3)
#TopoPRED(w_size=7, x_value=3)
#TopoPRED(w_size=11, x_value=3)
TopoPRED(w_size=15, x_value=3)
#TopoPRED(w_size=19, x_value=3)
#TopoPRED(w_size=23, x_value=3)
#TopoPRED(w_size=27, x_value=3)
#TopoPRED(w_size=31, x_value=3)
