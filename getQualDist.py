import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import gzip

def get_args():
    """ get user input for file, read length, and plot title"""
    parser=argparse.ArgumentParser(description="Generate and plot FASTQ quality stats")
    parser.add_argument("-f", "--file", help="path of FASTQ file", required=True, type=str)
    parser.add_argument("-l", '--readlen', help="length of reads", required=True, type=int)
    parser.add_argument("-r", "--reads", help="number of reads in the file", required=True, type=int)
    parser.add_argument("-t", "--title", help="title for plot", default="untitled plot", type=str)
    return(parser.parse_args())

def convert_phred(letter):
    """Converts a single character into a phred-33 score"""
    return (int(ord(letter)-33)) #I love you, ord

def init_ndarray(x,y,t):
    """ initialize desired size/shape/datatype of ndarray and return """
    return(np.zeros(shape=(x,y), dtype=t))

def get_scores(file,reads,readlen):
    """ given file, # of reads, and lenght of reads return 2darray of scores """
    i=0
    with gzip.open(file, 'rt') as f:
        content=f.readlines()
        scores=init_ndarray(readlen, reads, int)
        for line in content:
            if i%4==3:
                for j in range(readlen):
                    temp=i//4
                    scores[j,temp]=convert_phred(line[j])
            i+=1
    #print("ndarray built!")
    return(scores)    

def get_stats(scores,readlen):
    """ given ndarray of scores + length of reads, return ndarrays of mean, variance, stdev by position"""
    means=np.zeros(shape=readlen,dtype=float)
    var=np.zeros(shape=readlen,dtype=float)
    stdev=np.zeros(shape=readlen,dtype=float)
    for i in range(readlen):
        pos=scores[i]
        means[i]=pos.mean()
        var[i]=pos.var()
        stdev[i]=pos.std()
    #print("stats computed...")
    return(means,var,stdev)

def plot_dist(m,st,l,t):
    """ given arrays of mean (m) and stdev (st), read length l) and a title (t), plot it """
    x=range(l)
    y=range(40)
    fig=plt.figure(figsize=(16,8))
    ax1=fig.add_subplot(111)
    ax1.plot(m, label="mean score")
    ax1.plot(st,label="stdev")
    plt.legend(loc="center right");
    plt.title(t)
    plt.ylabel("Phred score")
    plt.xlabel("# Base Pair")
    plt.show()
    fig.savefig(t.replace(" ","_")+".pdf")
    return(True)
    
def do_all(file,reads,readlen,title):
    """ populate score array, calculate stats, and plot """
    scores=get_scores(file,reads,readlen)
    means,var,stdev=get_stats(scores,readlen)
    plot_dist(means,stdev,readlen,title)
    return(True)

a=get_args()
do_all(a.file,a.reads,a.readlen,a.title)
