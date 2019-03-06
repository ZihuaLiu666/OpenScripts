import argparse
import pandas as pd
import math
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

import plotly
plotly.tools.set_credentials_file(username='ZaihuNiu', api_key='sQTQMEgilz1NV2Y6sgnv')

import plotly.plotly as py
import plotly.figure_factory as ff

plt.switch_backend('agg')

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-f','--fasta', metavar='File',dest='inmatrix',help='Input file',type=open,required=True)
parser.add_argument('-p','--prop',metavar='File',dest='prop_output',help='Output prop file',type=argparse.FileType('w'),required=True)
parser.add_argument('-d','--decompos',metavar='File',dest='de_output',help='Output de pca file',type=argparse.FileType('w'),required=True)
parser.add_argument('-o','--output',metavar='File',dest='output',help='output file directory',type=str,required=True)
parser.add_argument('-t','--threshold',metavar='Int',dest='threshold',help='Output TPM file',type=int,default=500)
parser.add_argument('-ncp','--n_components',metavar='Int',dest='ncp',help='n_components',type=int,default=3)
args = parser.parse_args()
###### arguments ######

def D5hmCR_finder_by_gene(ifile, ofile, threshold, prop_otp, de_otp, ncp):

    data = pd.read_table(ifile, index_col=0, header=0)
    df = pd.DataFrame(data=data)
    ###########################
    df = df.iloc[:,0:-1]
    ###########################
    #Geneid Start End Sample1 Sample2 ... SampleN
    TotalCounts = [df[i].sum() for i in df.columns]
    column_index = [i for i in df.columns]

    c = 0
    newdf = pd.DataFrame()
    bin_log_over_one = {}

    while c < len(column_index):
        newdf[column_index[c]] = (df[column_index[c]]/TotalCounts[c])*1000000
        c += 1
    for i in newdf.index:
        if len(set(newdf.loc[i])) > 1:
            midT = math.log((newdf.loc[i]).std(), 10)
            if midT > 0:
                bin_log_over_one[i] = midT

    aa = sorted(bin_log_over_one.items(), key=lambda x: x[1], reverse=True)

    co = 0
    co_list = []

    for i in aa:
        if co < threshold:
            co_list.append(i[0])
            co += 1
        else:
            break
    #print(co_list)


    ############################################### PCA part
    pca_matrix = pd.DataFrame()
    for i in newdf.index:
        if i in co_list:
            pca_matrix[i] = df.loc[i]
    #print(pca_matrix)
    pca = PCA(n_components=ncp)
    new_pca_matrix = pca.fit_transform(pca_matrix)
    pca_p = pca.explained_variance_ratio_

    # cumulative plotting
    cul_pca_p = [sum(pca_p[0:i+1]) for i in range(len(pca_p))]
    for i in cul_pca_p:
        prop_otp.write('{}\n'.format(i))
    prop_otp.close()
    plt.bar(range(len(cul_pca_p)), cul_pca_p)
    #plt.show()
    plt.savefig('{}'.format(ofile))

    ####################################### Hierarchical Clustering
    de_otp.write('{}\n'.format('\t'.join([str(i) for i in column_index])))
    new_pca_matrix_tp = np.transpose(new_pca_matrix)
    np.savetxt(de_otp, new_pca_matrix_tp, delimiter='\t')
    dendro = ff.create_dendrogram(new_pca_matrix, labels=column_index)
    dendro['layout'].update({'width': 800, 'height': 800})
    py.iplot(dendro, filename=ofile)


D5hmCR_finder_by_gene(ifile=args.inmatrix, ofile=args.output, threshold = args.threshold,
                      prop_otp=args.prop_output, de_otp=args.de_output, ncp=args.ncp)
