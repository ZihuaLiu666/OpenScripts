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
parser.add_argument('-o','--output',metavar='File',dest='output',help='Output file',type=argparse.FileType('w'),required=True)
parser.add_argument('-b','--breakpoint',metavar='Int',dest='breakpoint',help='break point value',type=int,default=6)
args = parser.parse_args()
###### arguments ######

def D5hmCR_finder_by_gene(ifile, breakpointvalue, ofile):

    pca_matrix_name = ifile
    ifile = open(ifile)
    print('Working with TPM...')
    starttime = datetime.datetime.now()
    df = pd.read_csv(ifile, index_col=0, header=0, sep='\t')
    TotalCounts = [df[i].sum() for i in df.columns]
    column_index = [i for i in df.columns]
    #  Ofile writing
    cc = list(column_index)
    cc.insert(0, 'name')
    ofile.write('{}\n'.format('\t'.join(cc)))

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

    td = 5
    pair_extreme = remove_low_confidence_gene_counts(TPM, td)

    print('Calculating ...')
    starttime = datetime.datetime.now()
    linenum = TPM.shape[0]
    inum = 1
    writinglist = []
    for i in TPM.index:
        inum += 1
        outputpercentage(linenum, inum)
        if i.split('_')[0] not in ['chrX', 'chrY', 'chrM'] \
                and sum(TPM.loc[i]) != 0 \
                and row_compare(TPM.loc[i], pair_extreme) == 0:
            cross = Crossmeanstd(column_index, TPM.loc[i])
            if breakpointvalue in cross:
                body = np.array(TPM.loc[i]).tolist()
                body.insert(0, i)
                ofile.write('{}\n'.format('\t'.join([str(i) for i in body])))
            else:
                continue

    endtime = datetime.datetime.now()
    print('Calculation done ...')
    print('time use is {} seconds'.format((endtime - starttime).seconds))

def remove_low_confidence_gene_counts(df, td):
    column_index = [i for i in df.columns]
    pair_extreme = []
    for col_ind in column_index:
        # pair_extreme = [range(s1,e1), range(s2,e2) ...]
        s = np.percentile(list(set(df[col_ind])), td)
        e = np.percentile(list(set(df[col_ind])), 100 - td)
        pair_extreme.append(range(int(s),int(e)))
    return pair_extreme

def row_compare(line, pair_extreme):
    linelenrange = range(0, len(line))
    # Trick!!!
    a = sum([1 for i in linelenrange if line[i] != 0 and int(line[i]) in pair_extreme[i]])
    if a < len(line) - list(line).count(0):
        return 1  # if 'outlier' is in this line...
    else:
        return 0

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
    #print(outcluster)
    mean = defaultdict(list)
    std_mean = defaultdict(list)
    for k in outcluster.keys():
        mean[k] = np.mean([dfloc[i] for i in outcluster[k]])
        std_mean[k] = np.std([dfloc[i] for i in outcluster[k]]) \
                      + np.mean([dfloc[i] for i in outcluster[k]])
    sorted_TPM_list = sorted(mean.items(),
                             key=lambda x: x[1],
                             )
    #print(sorted_TPM_list)
    #print(mean)
    #print(std_mean)
    sort_name = [i[0] for i in sorted_TPM_list]
    #cross = {}
    cross = []
    c = 1
    breakpoint = 1
    while c < len(sort_name):
        value = sort_name[c]
        std_value = sort_name[c - 1]
        #cross[value] = (mean[value], std_mean[std_value])
        if mean[value] > std_mean[std_value]:
            #cross[value] = breakpoint
            cross.append(breakpoint)
            breakpoint += 1
        else:
            #cross[value] = 0
            cross.append(breakpoint)
        c += 1
    #print(cross)
    return cross

def outputpercentage(linenum, inum):
    if inum%10 == 0:
        percent = float(inum) * 100 / float(linenum)
        sys.stdout.write("%.2f" % percent)
        sys.stdout.write("%\r")
        sys.stdout.flush()


D5hmCR_finder_by_gene(ifile = args.inmatrix,
                      breakpointvalue = args.breakpoint,
                      ofile=args.output,
                      )
