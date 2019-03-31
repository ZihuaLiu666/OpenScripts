import argparse
import datetime
import sys
import pandas as pd
from collections import defaultdict
import numpy as np
np.set_printoptions(precision=5, suppress=True)

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inmatrix', metavar='File',dest='inmatrix',help='Input matrix file',type=str,required=True)
parser.add_argument('-om','--omatrix',metavar='File',dest='omatrix',help='Output matrix file',type=argparse.FileType('w'),required=True)
parser.add_argument('-b','--breakpoint',metavar='Int',dest='breakpoint',help='break point value',type=int,default=7)
args = parser.parse_args()
###### arguments ######

def D5hmCR_finder_by_gene(ifile, breakpointvalue, omtrix):

    ifile = open(ifile)
    print('Working with TPM...')
    starttime = datetime.datetime.now()
    df = pd.read_csv(ifile, index_col=0, header=0, sep='\t')
    TotalCounts = [df[i].sum() for i in df.columns]
    column_index = [i for i in df.columns]
    outcluster = cluster_gather(column_index)
    aa = sorted(outcluster.keys())
    cc = sorted(outcluster.keys())
    cc.insert(0, 'name')
    #  Ofile writing
    omtrix.write('{}\n'.format('\t'.join(cc)))

    #  Optimization preallocation
    c = 0
    TPM = pd.DataFrame()

    while c < len(column_index):
        TPM[column_index[c]] = (df[column_index[c]]/TotalCounts[c])*1000000
        c += 1
    endtime = datetime.datetime.now()
    print('TPM done ...')
    print('time use is {} seconds'.format((endtime - starttime).seconds))
    #np.savetxt('TPM_PCA_{}_gene_{}'.format(ncp, threshold), TPM, delimiter='\t')

    print('Calculating ...')
    starttime = datetime.datetime.now()
    linenum = TPM.shape[0]
    inum = 1
    for i in TPM.index:
        inum += 1
        outputpercentage(linenum, inum)
        cross = Crossmeanstd(column_index, TPM.loc[i])
        if sum(cross.values()) <= breakpointvalue:
            s = [cross[i] for i in aa]
            s.insert(0, i)
            omtrix.write('{}\n'.format('\t'.join([str(i) for i in s])))

    endtime = datetime.datetime.now()
    print('Calculation done ...')
    print('time use is {} seconds'.format((endtime - starttime).seconds))


def cluster_gather(clusterlist):
    outcluster = defaultdict(list)
    b = set([i.split('-')[0] for i in clusterlist])
    #print(b)
    #print(len(b))
    c = [i.split('-')[0] for i in clusterlist]
    for value in b:
        outcluster[value] = [i for i, x in enumerate(c) if x == value]
    return outcluster

def Crossmeanstd(column_index, dfloc):
    #print(dfloc[0])
    outcluster = cluster_gather(column_index)
    mean = defaultdict(list)
    std_mean = defaultdict(list)
    for k in outcluster.keys():
        mean[k] = np.mean([dfloc[i] for i in outcluster[k]])
        std_mean[k] = np.std([dfloc[i] for i in outcluster[k]]) \
                      + np.mean([dfloc[i] for i in outcluster[k]])
    sorted_TPM_list = sorted(mean.items(),
                             key=lambda x: x[1],
                             )
    sort_name = [i[0] for i in sorted_TPM_list]
    #print(sort_name)
    cross = {}
    for i in sort_name:
        cross[i] = 1
    #cross = []
    c = 1
    #breakpoint = 0
    while c < len(sort_name):
        cross[sort_name[0]] = 0
        value = sort_name[c]
        std_value = sort_name[c - 1]
        if mean[value] > std_mean[std_value]:
            break
        else:
            cross[value] = 0
            #cross.append(breakpoint)
        c += 1
    #print(cross.values())
    return cross

def outputpercentage(linenum, inum):
    if inum%10 == 0:
        percent = float(inum) * 100 / float(linenum)
        sys.stdout.write("%.2f" % percent)
        sys.stdout.write("%\r")
        sys.stdout.flush()

D5hmCR_finder_by_gene(ifile = args.inmatrix,
                      breakpointvalue = args.breakpoint,
                      omtrix=args.omatrix,
                      )
