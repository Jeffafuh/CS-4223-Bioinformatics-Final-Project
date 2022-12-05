import numpy as np
from random import sample

def runBasicMEME(seqs, motifL, motifStartingProb, threshold, initSamplingSize=1):
    """
    Runs the MEME algorithm on a set of sequences.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    motifL : int
        Length of the planted motif.
    motifStartingProb : float
        Probability of the motif occuring. (used to init the PWM).
    threshold : float
        The cutoff point of when to stop the algorithm.
    initSamplingSize: float, optional
        If provided, only searches through a random sampling of the unique subsequences proportional to the input value.

    Returns
    -------
    pwm : final PWM calculated
    z : final Z matrix calculated
    probNew: probability that the final PWM is the true PWM of the sequences
    """
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

def memeInitPWM(seqs, motifL, motifProb, initSamplingSize=1):
    """
    Search all the unique subsequences found throughout the sequences and select the subsequence
    that has the highest probability of occuring after on E/M step. This chosen subsequence is then
    used to initialize a PWM.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    motifL : int
        Length of the planted motif.
    motifProb : float
        Probability of the motif occuring. (used to init the PWM).
    initSamplingSize: float, optional
        If provided, only searches through a random sampling of the unique subsequences proportional to the input value.

    Returns
    -------
    bestPWM : initialized PWM to the best sequence found
    bestProb : probability that the best sequence found occurs in the sequences
    """
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

def getAllUniqueSubSeqs(seqs, length):
    """
    Enumerates all unique substrings found within the sequences of the specified length.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    length : int
        Length of the substrings to be enumerated.

    Returns
    ----------
    subseqs : the set of the unique substrings stored as tuples
    """
    subseqs = set()

    for i in range(seqs.shape[1]-length+1):
        sliced = seqs[:, i:i+length]
        for j in range(seqs.shape[0]):
            subseqs.add(tuple(sliced[j, :]))

    return subseqs

def calcLogLikelyhood(seqs, pwm, z):
    """
    Calculates the log probability ratio that a PWM is the true PWM of the sequences.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    pwm : ndarray
        Position-Weighted Matrix to be considered.
    z : ndarray
        Z matrix to use for calculations.

    Returns
    ----------
    prob : float of the calculated probability
    """
    prob = 0
    for i in range(seqs.shape[0]):
        posProb = 0
        for j in range(z.shape[1]):
            posProb += getProbInSeqLoc(seqs[i, :], pwm, j)
        prob += np.log(posProb/seqs.shape[0])

    return prob

def initPWMForSeq(seq, x):
    """
    Initializes PWM to a given sequence.

    For i > 0, PWM[c, i] = x if c occurs at seq[i-1] and PWM[c, i] = 1-x otherwise.

    The background probabilities (i=0) are initialized to a uniform distribution.

    Parameters
    ----------
    seq : ndarray
        Sequence to initialize the PWM to.
    x : float
        Probability of the motif occuring.

    Returns
    ----------
    pwm : ndarray representing the PWM
    """
    pwm = np.zeros((4, seq.shape[0]+1)) + (1-x)/3
    pwm[:, 0] = 0.25
    for i in range(seq.shape[0]):
        pwm[seq[i], i+1] = x
    
    return pwm

def m_step(seqs, pwm, z):
    """
    Updates the PWM (in-place) given a Z matrix.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    pwm : ndarray
        Position-Weighted Matrix to be updated.
    z : ndarray
        Z matrix to use for calculations.
    """

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

def initZMatrix(seqs, pwm):
    """
    Initializes an empty Z matrix to the correct shape based of the PWM and sequences' shape.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    pwm : ndarray
        Position-Weighted Matrix.

    Returns
    ----------
    All-zero ndarray representing the Z matrix
    """
    return np.zeros((seqs.shape[0], seqs.shape[1]-pwm.shape[1]+2))

def e_step(seqs, pwm, z):
    """
    Updates the Z matrix (in-place) given a PWM.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    pwm : ndarray
        Position-Weighted Matrix to use for calculations.
    z : ndarray
        Z matrix to be updated.
    """
    for i in range(seqs.shape[0]):
        getProbSequence(seqs[i, :], pwm, z[i, :])

def getProbSequence(seq, pwm, subZ):
    """
    Updates the one row of the Z matrix (in-place) given a PWM.

    Parameters
    ----------
    seq : ndarray
        Sequence to be considered for the Z row.
    pwm : ndarray
        Position-Weighted Matrix to use for calculations.
    subZ : ndarray
        Row of the Z matrix to be updated.
    """
    for i in range(subZ.shape[0]):
        subZ[i] = getProbInSeqLoc(seq, pwm, i)
    #normalize row
    subZ /= subZ.sum()

def getProbInSeqLoc(seq, pwm, loc):
    """
    Calculates the probability that the motif starts at the specified location in a sequence given a PWM.

    Parameters
    ----------
    seq : ndarray
        Sequence to be considered.
    pwm : ndarray
        Position-Weighted Matrix to use for calculations.
    loc : int
        Starting location of the motif.

    Returns
    -------
    prob : the probability as a single float.
    """
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

def getKmers(seqs, locs, motifL):
    """
    Extract the kmers of length motifL from the set of sequences at the specified locations.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    locs : ndarray
        Locations to extract the kmers from.
    motifL : int
        Length of the motif to be extracted.

    Returns
    -------
    kmers : ndarray of the kmers extracted from the sequences.
    """
    kmers = np.zeros((seqs.shape[0], motifL))
    for x in range(locs.shape[0]):
        kmers[x, :] = seqs[x, locs[x]:locs[x]+motifL]
    
    return kmers

def getZLocs(z):
    """
    Gets the position with the highest probability of the motif occuring for each sequence.

    Parameters
    ----------
    z : ndarray
        Z matrix containing the probabilities.

    Returns
    ----------
    An ndarray containing the location of the motif for each sequence.
    """
    return np.argmax(z, axis=1)