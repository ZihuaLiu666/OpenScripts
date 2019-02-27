import argparse
import pysam
import re

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-b','--bam', metavar='File',dest='bam',help='Input file',required=True)
parser.add_argument('-f','--fasta', metavar='File',dest='fasta_file',help='Input file',type=open,required=True)
parser.add_argument('-o','--Output',metavar='File',dest='output',help='Output file',required=True)
args = parser.parse_args()
###### arguments ######

# IF YOU HAVE ANY QUESTIONS ABOUT THE SCRIPTS
# PLEASE FEEL FREE TO CONTACT WITH THE AUTHOR
# E-MAIL: zihua.liu666@gmail.com
# HOPE YOU CAN ENJOY THIS TOY!


def product_list(fasta_file):
    aDict={}
    for line in fasta_file:
        if line[0]=='>':
            key=line.strip().split()[0][1::]
            aDict[key]=[]
        else:
            aDict.setdefault(key,[]).append(line.strip())

    long_string_Dict = {}
    for k in aDict.keys():
        value = ''.join(list(aDict[k]))
        long_string_Dict[k] = value
    print('Compile reference fasta file completed!')
    return long_string_Dict

################################################################################################

def order_inser(cigar):
    order_event = [int(i) for i in re.findall("\d+", cigar)]
    all_M = order_event[0:100:2]
    all_event = order_event[1:100:2]
    # print(order_I, all_M, all_I)
    i = 0
    order_insertion = []
    while i < len(all_event):
        order_insertion.append((sum(all_M[0:i + 1]), all_event[i]))
        i += 1
    return order_insertion
        # print(order_insertion)


def deduction_D_I(cigar):
    # cigar = '2M1D2M1I2M'
    if cigar.find('I') > cigar.find('D'):
        order_I = [int(i) for i in re.findall("\d+", cigar.split('I')[0])]
        order_D = [int(i) for i in re.findall("\d+", cigar.split('D')[1])]
        newI = str(sum(order_I[0:-1])) + 'M' + str(order_I[-1]) + 'I' + cigar.split('I')[1]
        newD = cigar.split('D')[0] + 'D' + str(sum(order_D)) + 'M'
        #print(newD, newI)
        return newI, newD
    else:
        order_D = [int(i) for i in re.findall("\d+", cigar.split('D')[0])]
        order_I = [int(i) for i in re.findall("\d+", cigar.split('I')[1])]
        newD = str(sum(order_D[0:-1])) + 'M' + str(order_D[-1]) + 'D' + cigar.split('D')[1]
        newI = cigar.split('I')[0] + 'I' + str(sum(order_I)) + 'M'
        #print(newD, newI)
        return newI, newD


def add_N(ori_sequence, pos_info):
    new_sequence = ''
    pos_info.insert(0, (0, 0))
    i = 0
    while i < len(pos_info) - 1:
        new_sequence += ori_sequence[pos_info[i][0]:pos_info[i + 1][0]] + pos_info[i + 1][1] * 'N'
        i += 1
    new_sequence += ori_sequence[pos_info[-1][0]::]
    # print(new_sequence)
    return new_sequence


def trans_query_reference_by_cigar(query, reference, cigar):
    # I or D?
    query_I = [i.start() - 1 for i in re.finditer('I', cigar)]
    query_D = [i.start() - 1 for i in re.finditer('D', cigar)]

    order_insertion = order_inser(cigar)

    if query_D == []:
        new_reference = add_N(reference, order_insertion)
        #print(query, new_reference)
        return query, new_reference
        # print(new_sequence)
    elif query_I == []:
        new_query = add_N(query, order_insertion)
        #print(new_query, reference)
        return new_query, reference
    else:
        newIcigar, newDcigar = deduction_D_I(cigar)
        newI_oi = order_inser(newIcigar)
        newD_oi = order_inser(newDcigar)
        new_query = add_N(query, newD_oi)
        new_reference = add_N(reference, newI_oi)
        #print(new_query, new_reference)
        return new_query, new_reference


################################################################################################

def readbam(bam_file, fasta_file, wfile):

    genomeDict = product_list(fasta_file)
    #print(genomeDict['chr10'][3244897:3244897+153])

    bamfile = pysam.AlignmentFile(bam_file,"rb")
    ofile = pysam.AlignmentFile(wfile, "w", template=bamfile)

    for r in bamfile:
        query = r.seq
        reference = genomeDict[r.reference_name][r.pos:r.pos+len(query)+2]
        reference = reference.upper()
        #print(reference)
        #print(r.cigarstring)
        # THIS SCRIPT HAS NO SENSITIVITY FOR DELETION OR INSERTION CONDITIONS
        if 'I' in r.cigarstring or 'D' in r.cigarstring:
            new_query, new_reference = trans_query_reference_by_cigar(query, reference, r.cigarstring)
            #if r.cigarstring == '81M4I65M':
                #print(new_query)
                #print(new_reference)
            #if discard_TH_and_DA(new_query, new_reference) == 1:
            if discard_TH_and_DA(new_query, new_reference) != 1:
                ofile.write(r)
            #continue
        else:
            #if discard_TH_and_DA(query, reference) == 1:
            if discard_TH_and_DA(query, reference) != 1:
                ofile.write(r)
    print('Congratulations! You can move to the next steps now!')

def discard_TH_and_DA(query, reference):
    # H represents A/C/T NOT G
    query_T = [i.start() for i in re.finditer('T', query)]
    reference_C = [j.start() for j in re.finditer('C', reference)]
    # D represents A/G/T NOT C
    query_A = [i.start() for i in re.finditer('A', query)]

    # Position overlapping
    plus_T = set(query_T).intersection(set(reference_C))
    base_behind_plus_T = [query[i+1] for i in plus_T if i+1 < len(query) and len(plus_T) >= 3]
    #print(base_behind_plus_T)
    plus_A = set(query_A).intersection(set(reference_C))
    base_behind_plus_A = [query[i - 1] for i in plus_A if i > 0 and len(plus_A) >= 3]

    judge_T = len(base_behind_plus_T) - base_behind_plus_T.count('G')
    judge_A = len(base_behind_plus_A) - base_behind_plus_A.count('C')

    if judge_T >=3 or judge_A >=3: return 1
    else: return 0

readbam(bam_file = args.bam, fasta_file=args.fasta_file, wfile=args.output)
