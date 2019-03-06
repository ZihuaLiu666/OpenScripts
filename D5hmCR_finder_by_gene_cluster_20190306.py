import argparse
import pandas as pd
import math
import numpy as np
np.set_printoptions(precision=8, suppress=True)
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram

plt.switch_backend('agg')

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inmatrix', metavar='File',dest='inmatrix',help='Input matrix file',type=open,required=True)
parser.add_argument('-t','--threshold',metavar='Int',dest='threshold',help='Output TPM file',type=int,default=500)
parser.add_argument('-ncp','--n_components',metavar='Int',dest='ncp',help='n_components',type=int,default=9)
args = parser.parse_args()
###### arguments ######

def D5hmCR_finder_by_gene(ifile, threshold, ncp):

    data = pd.read_csv(ifile, index_col=0, header=0, sep='\t')
    df = pd.DataFrame(data=data)
    df = df.iloc[:, 0:-1]
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

    sorted_newdf_list = sorted(bin_log_over_one.items(),
                               key=lambda x: x[1],
                               reverse=True,
                               )

    co = 0
    co_list = []

    for i in sorted_newdf_list:
        if co < threshold:
            co_list.append(i[0])
            co += 1
        else:break

    ############################################### PCA part
    pca_matrix = pd.DataFrame()
    for i in newdf.index:
        if i in co_list:
            pca_matrix[i] = newdf.loc[i]      #20190306
    pca = PCA(n_components=ncp)
    new_pca_matrix = pca.fit_transform(pca_matrix)
    pca_p = pca.explained_variance_ratio_

    # cumulative plotting
    cul_pca_p = [sum(pca_p[0:i+1]) for i in range(len(pca_p))]
    np.savetxt('CulPCA_{}.txt'.format(ncp), np.array(cul_pca_p), delimiter='\t')

    ####################################### Hierarchical Clustering

    plt.figure(figsize=(12, 5))
    mergings = linkage(pca_matrix, method='complete', metric='euclidean')
    dendrogram(mergings,
               labels=column_index,
               leaf_rotation=45,
               leaf_font_size=5,
               )
    plt.savefig('hierarchical_cluster_PCA_{}_gene_{}.pdf'.format(ncp, threshold))


D5hmCR_finder_by_gene(ifile=args.inmatrix,
                      threshold = args.threshold,
                      ncp=args.ncp,
                      )
