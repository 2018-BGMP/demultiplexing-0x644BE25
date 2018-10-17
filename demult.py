import gzip
import time
from contextlib import ExitStack
import argparse
start_time = time.time()

def get_args():
    """ get command line input """
    parser=argparse.ArgumentParser(description="Generate and plot assembly stats")
    parser.add_argument("-i1", "--index1", help="index read 1 fastq.gz file", required=True, type=str)
    parser.add_argument("-i2", '--index2', help="index read 2 fastq.gz file", required=True, type=str)
    parser.add_argument("-r1", "--read1", help="bio read 1 fastq.gz file", required=True, type=str)
    parser.add_argument("-r2", '--read2', help="bio read 2 fastq.gz file", required=True, type=str)
    parser.add_argument("-ix", "--indices", help="file of indices", required=True, type=str)
    parser.add_argument("-qc", "--quality", help="quality cutoff for indices",default=15, type=int)
    return(parser.parse_args())

def revComp(seq):
    """ return reverse complement of sequence """
    res =""
    for nt in seq:
        res = rcDict.get(nt) + res
    return(res)

def buildIndexDict(file):
    """ given tsv file of index names + indices, return dict of pairs """
    f = open(file, 'r')
    toName = {}
    for l in f.readlines():
        curr = l.split("\t")
        toName[curr[1].strip("\n")]=curr[0].strip("\n")
    return(toName)

def buildToFile(myDict):
    res = {}
    i=0
    for val in myDict.values():
        print(val)
        res[val]=i
        i+=1
    res['BAD']=i
    return(res)

def meanQual(seq):
    length = len(seq)
    score = 0
    for i in range (length):
        score += ord(seq[i])-33
    return score/length

def getIndex(seq1,seq2):
    """ figure out how to route read pair"""
    global nCount, hopCount, badIndexCount, goodCount, ixDict
    # check for N in read
    if 'N' in seq1 or 'N' in seq2:
        nCount +=1
        return('BAD')
    
    # lookup first index
    ix1 = ixDict.get(seq1)
    
    # if good, get lookup rev comp for second index
    if ix1 != None:
        ix2 = ixDict.get(revComp(seq2))
        
    # if not good, use rev comp for first, normal for second
    else:
        ix1 = ixDict.get(revComp(seq1))
        ix2 = ixDict.get(seq2)
        
    # if one had no lookup = bad index
    if ix1 == None or ix2 == None:
        badIndexCount +=1
        return('BAD')
    
    # if they don't match = hopping
    if ix1 != ix2:
        hopCount +=1
        return('BAD')
    
    # if we make it here, the pair is good
    else:
        goodCount +=1 
        return(ix1)

def writeRead(index,ix1,ix2):
    # get our global vars
    global r1, r2, r1List, r2List, toFile
    
    # get index of file 
    i = toFile.get(index)

    # write read 1
    file = r1List[i]
    header = r1.readline().strip()+':'+ix1+'\n'
    file.write(header)
    file.write(r1.readline())
    file.write(r1.readline())
    file.write(r1.readline())
       
    # write read 2
    file = r2List[i]
    header = r2.readline().strip()+':'+ix2+'\n'
    file.write(header)
    file.write(r2.readline())
    file.write(r2.readline())
    file.write(r2.readline())
    

    pass

def openEm():
    global cutoff, lowQualCount
    index1 = gzip.open(args.index1, 'rt')
    index2 = gzip.open(args.index2, 'rt')
    for val in ixDict.values():
        r1List
    # make our output files
    for line in index1:
        # get the rest of our index1  entry
        header1 = line.strip()
        #print(header1)
        ix1 = index1.readline().strip()
        #print(ix1)
        sep1 = index1.readline().strip()
        #print(sep1)
        qual1 = index1.readline().strip()
        #print(qual1)
        
        # get index2 entry
        header2 = index2.readline().strip()
        ix2 = index2.readline().strip()
        sep2 = index2.readline().strip()
        qual2 = index2.readline().strip()
        
        # check for low quality
        if min(meanQual(qual1),meanQual(qual2)) < cutoff:
            lowQualCount +=1
            indexName = ('BAD')
        else:
            # find correct index file
            indexName = getIndex(ix1,ix2)
        
        # write out the read
        #print(indexName)
        writeRead(indexName,ix1,ix2)
    pass
    
    resFile = open('demult_result.txt', 'w+')
    resFile.write('index pairs containing N: '+str(nCount)+'\n')
    resFile.write('hopped pairs: '+str(hopCount)+'\n')
    resFile.write('low quality pairs: '+str(lowQualCount)+'\n')
    resFile.write('invalid indexes: '+str(badIndexCount)+'\n')
    resFile.write('good index pairs: '+str(goodCount)+'\n')
    total = nCount+hopCount+lowQualCount+badIndexCount+goodCount
    resFile.write('TOTAL: '+str(total))
    
def closeFiles(dictionary):
    for key, val in dictionary.items():
        val.close()
    resFile.close()

# get our passed arguments
args=get_args()

# make some global vars
nCount, hopCount, lowQualCount, badIndexCount, goodCount = 0,0,0,0,0
rcDict = dict([('A','T'),('T','A'),('C','G'),('G','C'),('N','N')])

r1 = gzip.open(args.read1, 'rt')
r2 = gzip.open(args.read2, 'rt')
ixDict = buildIndexDict(args.indices)
toFile = buildToFile(ixDict)
cutoff = int(args.quality)

with ExitStack() as stack:
    r1List = [stack.enter_context(open(fname+'_R1.fastq', 'w+')) for fname in toFile.keys()]
    r2List = [stack.enter_context(open(fname+'_R2.fastq', 'w+')) for fname in toFile.keys()]
    openEm()

print("--- %s seconds ---" % (time.time() - start_time))
