import numpy as np
import pandas as pd
import matplotlib as plt

#Given a set of kmers, generate the consensus position-weight matrix
def getNewWeights(kmers):
    pwm = np.zeros((4, kmers.shape[1]))

    for x in range(4):
        pwm[x, :] = np.count_nonzero(kmers == x, axis=0) / kmers.shape[0]

    return pwm

#Init weights to a random sampling of the sequences
def initWeightsSample(sequences, motifL):
    locs = np.random.randint(0, sequences.shape[1] - motifL+1, (sequences.shape[0],))
    randomSample = getKmers(sequences, locs, motifL)
    weights = getNewWeights(randomSample)
    weights = np.concatenate((np.zeros((4,1))+0.25,weights), axis=1)

    return weights

#Init the weights to a uniform probability
def initWeightsUniform(motifL):
    return np.zeros((4, motifL))+0.25

#Given a bunch of sequences and locations, extract the kmers of length k from those locations and return them
def getKmers(seqs, locs, motifL):
    kmers = np.zeros((seqs.shape[0], motifL))
    for x in range(locs.shape[0]):
        kmers[x, :] = seqs[x, locs[x]:locs[x]+motifL]
    
    return kmers

#Given a bunch of sequences, find the match of the weight matrix on each seqeunce
def findBestMatches(sequences, pwm, baseLine, pM):
    probs = np.zeros(sequences.shape[0])-1
    locs = np.zeros(sequences.shape[0], dtype=int)
    for x in range(sequences.shape[0]):
        (locs[x], probs[x]) = findBestMatch(sequences[x, :], pwm, baseLine, pM)
    
    return (locs, probs)

#Given a position-weight matrix and (long) sequence, find the segment (best match) that has the highest probability of occuring
def findBestMatch(sequence, pwm, baseLine, pM):
    prob = -1
    loc = 0
    for x in range(len(sequence)-pwm.shape[1]+1):
        tempProb = getProbability(sequence[x:x+pwm.shape[1]], pwm, baseLine, pM)
        if tempProb > prob:
            prob = tempProb
            loc = x

    return (loc, prob)

#Given a kmer and position-weight matrix, calculate the probability of that kmer appearing given the matrix
def getProbability(kmer, pwm, baseLine, pM):
    pWGivenM = 1
    for x in range(pwm.shape[1]):
        pWGivenM *= pwm[kmer[x]][x]
    pWGivenB = 1

    for x in range(len(kmer)):
        pWGivenB *= baseLine[kmer[x]]

    #pM=0.1
    pB=1-pM
    pW = pWGivenM*pM + pWGivenB*pB
    pMGivenW = pWGivenM*pM / pW

    return pMGivenW