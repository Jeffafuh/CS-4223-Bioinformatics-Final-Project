from MEMEFunctions import *
from testingFunctions import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pytictoc import TicToc
#default params are nSeqs=20, seqL=20, k=7, nMuts=1

def testmotifL():
    params = np.arange(5, 16, 1)
    trials = params.shape[0]
    runtimes = np.zeros((trials,))
    timer = TicToc()

    for i in range(trials):
        seqL =20
        nSeqs = 20
        motifL = params[i]
        nMuts = 1
        (seqs, consensus, locs) = generateSeqsWithMotif(nSeqs, seqL, motifL, nMuts)

        timer.tic()
        runBasicMEME(seqs, motifL, 0.7, 0.001)
        runtimes[i] = timer.tocvalue()

    fig = plt.figure()
    plt.plot(params, runtimes, marker='o')
    plt.ylabel('Runtime (seconds)')
    plt.xlabel('motifL')
    plt.title('MEME Runtime by Motif Length')
    plt.xticks(params)
    fig.tight_layout()
    fig.savefig('./images/motifL.jpeg', dpi=500)

    print('Done!')

def testseqL():
    params = np.arange(10, 31, 1)
    trials = params.shape[0]
    runtimes = np.zeros((trials,))
    timer = TicToc()

    for i in range(trials):
        seqL = params[i]
        nSeqs = 20
        motifL = 7
        nMuts = 1
        (seqs, consensus, locs) = generateSeqsWithMotif(nSeqs, seqL, motifL, nMuts)

        timer.tic()
        runBasicMEME(seqs, motifL, 0.7, 0.001)
        runtimes[i] = timer.tocvalue()

    fig = plt.figure()
    plt.plot(params, runtimes, marker='o')
    plt.ylabel('Runtime (seconds)')
    plt.xlabel('seqL')
    plt.title('MEME Runtime by Sequence Length')
    plt.xticks(params)
    fig.tight_layout()
    fig.savefig('./images/seqL.jpeg', dpi=500)

    print('Done!')

def testnSeqs():
    params = np.arange(10, 31, 1)
    trials = params.shape[0]
    runtimes = np.zeros((trials,))
    timer = TicToc()

    for i in range(trials):
        seqL = 20
        nSeqs = params[i]
        motifL = 7
        nMuts = 1
        (seqs, consensus, locs) = generateSeqsWithMotif(nSeqs, seqL, motifL, nMuts)

        timer.tic()
        runBasicMEME(seqs, motifL, 0.7, 0.001)
        runtimes[i] = timer.tocvalue()

    fig = plt.figure()
    plt.plot(params, runtimes, marker='o')
    plt.ylabel('Runtime (seconds)')
    plt.xlabel('nSeqs')
    plt.title('MEME Runtime by Number of Sequences')
    plt.xticks(params)
    fig.tight_layout()
    fig.savefig('./images/nSeqs.jpeg', dpi=500)

    print('Done!')

testmotifL()