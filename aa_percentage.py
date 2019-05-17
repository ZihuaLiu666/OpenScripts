import argparse
###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-i','--ifile', metavar='File',dest='ifile',help='Input file',type=open,required=True)
parser.add_argument('-r','--region',metavar='Str',dest='region',help='region',type=str,required=True)
args = parser.parse_args()
###### arguments ######

def product_list(ifile, iregion):
    aDict={}
    for line in ifile:
        if line[0]=='>':
            key=line.strip().split()[0]
            aDict[key]=[]
        else:
            aDict.setdefault(key,[]).append(line.strip())
    s=int(iregion.split(':')[0])
    e=int(iregion.split(':')[1])+1
    for k in aDict.keys():
        value = ''.join(list(aDict[k]))
        cc=value[s:e]

    single_cc = {}
    for i in cc:
        if i not in single_cc.keys():
            single_cc[i] = 0
    for j in cc:
        single_cc[j] += 1

    lzh = sorted(single_cc.items(),
                             key=lambda x: x[1],
                             reverse=True
                             )
    i = 0
    while i < 6:
        print('{}: {}'.format(lzh[i][0], round(lzh[i][1]/e, 3)))
        i += 1

product_list(ifile=args.ifile,
             iregion=args.region)
