in_file = '../data/tm_alpha_beta_3state.3line.txt'
out_dir = '../bin/'

#iterates through input file, parsing sequences and topologies
f=open(in_file, 'r')
f=f.readlines()

#iterates through input file, parsing sequences and topologies according to
#given window size
IDs=[]
sequences=[]
g=0
o=0
p=0
for i in range(0, len(f), 3):
    if 'G' in f[i+2]:
        g+=1
        continue
    else:
        IDs.append(f[i].split("|", 1)[0].rstrip().strip('>'))
        p+=1
#QUESTION:WHAT'S THE DIFFERENCE BETWEEN 'PASS' AND 'CONTINUE'?

for i in range(1, len(f), 3):
    if 'G' in f[i+1]:
        continue
    else:
        sequences.append(f[i].rstrip('\n'))

o=g+p   
print('Your test file contains {} protein sequences.'.format(o))
print('In total, {} membrane proteins will be used in the final dataset.'.format(p))
print('{} sequences refer to globular proteins and were removed from the final dataset.'.format(g))
print('TopoPRED will now prepare your sequences to PSI-BLAST. Please wait ...')

#RUNS PSI-BLAST ON INPUT PROTEIN SEQUENCES
for i in range(len(IDs)):
    with open((out_dir + str(IDs[i]) + '.fa'), 'w') as fa:
        fa.write(str(IDs[i]) + '\n')
        fa.write(str(sequences[i]) + '\n')

    ####GO TO BASH HERE####
