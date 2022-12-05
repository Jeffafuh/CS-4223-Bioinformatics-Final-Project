import numpy as np
from random import sample

""" TODO: generating the unique subsequences isn't the bottle neck, taking a step with them is
          so maybe make an alternate memeInitPWM that only considers a random subset of the subsequences
          maybe like only 1/4 or 1/5 of them?
"""

def getZLocs(z):
    return np.argmax(z, axis=1)

def runBasicMEME(seqs, motifL, motifStartingProb, threshold, initSamplingSize):
    pwm, probCur = memeInitPWM(seqs, motifL, motifStartingProb, initSamplingSize)
    z = initZMatrix(seqs, pwm)
    e_step(seqs, pwm, z)
    m_step(seqs, pwm, z)
    probNew = calcLogLikelyhood(seqs, pwm, z)
    while abs(probCur-probNew) > threshold:
        probCur = probNew
        e_step(seqs, pwm, z)
        m_step(seqs, pwm, z)
        probNew = calcLogLikelyhood(seqs, pwm, z)
    
    return (pwm, z, probNew)

def memeInitPWM(seqs, motifL, motifProb, initSamplingSize):
    subseqs = getAllUniqueSubSeqs(seqs, motifL)
    subseqs = sample(subseqs, int(len(subseqs) * initSamplingSize))
    bestPWM = None
    bestProb = -np.inf
    for subseq in subseqs:
        pwm = initPWMForSeq(np.array(subseq), motifProb)
        z = initZMatrix(seqs, pwm)
        e_step(seqs, pwm, z)
        m_step(seqs, pwm, z)
        likelyhood = calcLogLikelyhood(seqs, pwm, z)
        if likelyhood > bestProb:
            bestProb = likelyhood
            bestPWM = pwm

    return (bestPWM, bestProb)

#Given a set of sequences, generate the set of all unique substrings of size length found in the sequences
def getAllUniqueSubSeqs(seqs, length):
    subseqs = set()

    for i in range(seqs.shape[1]-length+1):
        sliced = seqs[:, i:i+length]
        for j in range(seqs.shape[0]):
            subseqs.add(tuple(sliced[j, :]))

    return subseqs

#Given a set of sequences, a pwm, and a z matrix, calculate to log probability ratio that this pwm is the true pwm
def calcLogLikelyhood(seqs, pwm, z):
    prob = 0
    for i in range(seqs.shape[0]):
        posProb = 0
        for j in range(z.shape[1]):
            posProb += getProbInSeqLoc(seqs[i, :], pwm, j)
        prob += np.log(posProb/seqs.shape[0])

    return prob

#Initialize a pwm for a given input sequence where each symbol in the sequence has a probability of x ocurring
def initPWMForSeq(seq, x):
    pwm = np.zeros((4, seq.shape[0]+1)) + (1-x)/3
    pwm[:, 0] = 0.25
    for i in range(seq.shape[0]):
        pwm[seq[i], i+1] = x
    
    return pwm

#Re-estimate the pwm matrix given a new Z matrix
def m_step(seqs, pwm, z):
    #motif probabilites (all but last column)
    for i in range(1, pwm.shape[1]-1):
        for j in range(4):
            pwm[j,i] = (z[seqs[:, i-1:i-pwm.shape[1]+1] == j].sum())

    #motif probabilities for last column
    for j in range(4):
        pwm[j, -1] = (z[seqs[:, pwm.shape[1]-1-1:] == j].sum())

    #background probabilities
    for j in range(4):
        pwm[j, 0] = (np.count_nonzero(seqs == j) - np.sum(pwm[j, 1:]))
    
    #The 'random' +4 to the denominator and +1 to each numerator is for smoothing out of probs
    pwm += 1
    pwm[:, 1:] /= seqs.shape[0]+4 #equivalent to doing z.sum()+4
    pwm[:, 0] /= np.sum(pwm[:, 0])

#Initialize the Z matrix to the appropriate shape
def initZMatrix(seqs, pwm):
    return np.zeros((seqs.shape[0], seqs.shape[1]-pwm.shape[1]+2))

#Re-estimate the Z matrix given a new weight matrix
# z = np.zeros((seqs.shape[0], seqs.shape[1]-pwm.shape[1]+2))
def e_step(seqs, pwm, z):
    for i in range(seqs.shape[0]):
        getProbSequence(seqs[i, :], pwm, z[i, :])

#Given a sequence and a weight matrix, update the Z matrix in-place to the appropriate values/probabilities
def getProbSequence(seq, pwm, subZ):
    for i in range(subZ.shape[0]):
        subZ[i] = getProbInSeqLoc(seq, pwm, i)
    #normalize row
    subZ /= subZ.sum()

#Given a sequence, weight matrix, and a hypothesized location where the motif is, calculate the probability of this occuring
def getProbInSeqLoc(seq, pwm, loc):
    prob = 1

    #get probabilities before the hypothesized start location in the background
    for i in range(loc):
        prob *= pwm[seq[i]][0]

    #get probabilites for the motif
    for i in range(pwm.shape[1]-1):
        prob *= pwm[seq[i+loc]][i+1]

    #get background prob for rest of the seq
    for i in range(loc+pwm.shape[1]-1, seq.shape[0]):
        prob *= pwm[seq[i]][0]

    return prob

#Given a bunch of sequences and locations, extract the kmers of length k from those locations and return them
def getKmers(seqs, locs, motifL):
    kmers = np.zeros((seqs.shape[0], motifL))
    for x in range(locs.shape[0]):
        kmers[x, :] = seqs[x, locs[x]:locs[x]+motifL]
    
    return kmers