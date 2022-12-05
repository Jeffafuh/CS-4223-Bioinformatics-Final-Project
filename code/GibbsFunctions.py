import numpy as np

def runGibbs(seqs, motifL, minIters, maxIters):
    """
    Runs the Gibbs Sampling algorithm on a set of sequences.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    motifL : int
        Length of the planted motif.
    minIters : int
        Minimum number of iterations to run the algorithm for.
    maxIters : int
        Maximum number of iterations to run the algorithm for.

    Returns
    -------
    pwm : final PWM calculated.
    locsNew : the final locations of the predicted motifs.
    i: number of iterations the algorithm performed.
    """
    pwm = initPWMMatrix(motifL)
    #Do one iteration, then enter while-loop
    locsCur = initRandomSamplePos(seqs, motifL)
    locsNew = locsCur.copy()
    idxToRemove = np.random.randint(seqs.shape[0])
    predict_step(seqs, pwm, locsCur, idxToRemove)
    sampling_step(seqs, pwm, locsNew, idxToRemove)
    
    i = 1
    while (np.any(locsCur != locsNew) or i < minIters) and i < maxIters:
        locsCur = locsNew.copy()
        idxToRemove = np.random.randint(seqs.shape[0])
        predict_step(seqs, pwm, locsCur, idxToRemove)
        sampling_step(seqs, pwm, locsNew, idxToRemove)
        i += 1

    #print('Total iters = ', i)
    return (pwm, locsNew, i)


def initPWMMatrix(motifL):
    """
    Initializes an empty PWM to the correct shape based on the motif length.

    Parameters
    ----------
    motifL : int
        Length of the motif to be considered.

    Returns
    ----------
    All-zero ndarray representing the PWM.
    """
    return np.zeros((4, motifL+1))

def initRandomSamplePos(sequences, motifL):
    """
    Initializes the predicted motif locations to random locations.

    Parameters
    ----------
    sequences : ndarray
        Array representing the set of sequences.
    motifL : int
        Length of the motif to be considered.

    Returns
    ----------
    locs: ndarray of the randomized locations for each sequence.
    """
    locs = np.random.randint(0, sequences.shape[1] - motifL+1, (sequences.shape[0],))

    return locs

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

def getBackground(seqs, locs, motifL):
    """
    Extract the characters AROUND the kmers of length motifL from the set of sequences at the specified locations.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    locs : ndarray
        Locations to of the kmers to cut out.
    motifL : int
        Length of the kmers to be extracted.

    Returns
    -------
    background : ndarray the sequences without the specified kmers.
    """
    background = np.zeros((seqs.shape[0], seqs.shape[1] - motifL))
    for i in range(seqs.shape[0]):
        background[i, :] = np.delete(seqs[i, :], np.arange(locs[i], locs[i]+motifL))
    return background
    
def predict_step(seqs, pwm, locs, idxToRemove):
    """
    Given a sequence to remove, remove the sequence from seqs then recalculate the pwm matrix.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    pwm : ndarray
        Position-Weighted Matrix to be updated.
    locs : ndarray
        Locations of all the predicted motif locations.
    idxToRemove : int
        Index of the sequence to remove.
    """
    motifL = pwm.shape[1] - 1
    seqs = np.delete(seqs, idxToRemove, axis=0)
    locs = np.delete(locs, idxToRemove)

    kmers = getKmers(seqs, locs, motifL)
    background = getBackground(seqs, locs, motifL)

    #loop through the characters and count them up
    for i in range(4):
        #kmers
        pwm[i, 1:] = np.count_nonzero(kmers == i, axis=0)
        #background
        pwm[i, 0] = np.count_nonzero(background == i)

    #For laplace smoothing to make sure no one prob is ever too close to 0, add 1 and divide by +4
    pwm += 1
    pwm[:, 1:] /= seqs.shape[0] + 4
    pwm[:, 0] /= seqs.shape[0] * (seqs.shape[1]-motifL) + 4

def getDistributionForSeq(seq, pwm):
    """
    Calculate the distribution array for a sequence given a PWM.
    This distribution is the likelyhood of the PWM occuring at each possible loc in the sequence.

    Parameters
    ----------
    seq : ndarray
        The sequence to calculate the distribution for.
    pwm : ndarray
        Position-Weighted Matrix to use for calculations.

    Returns
    ----------
    distribution : ndarray of the motif distribution for the sequence.
    """
    distribution = np.zeros((seq.shape[0] - pwm.shape[1] + 2,))

    for i in range(distribution.shape[0]):
        distribution[i] = np.exp(calcLikelyhoodOfSeq(seq, pwm, i))
    
    #normalize
    distribution /= np.sum(distribution)

    return distribution

def calcLikelyhoodOfSeq(seq, pwm, startLoc):
    """
    Calculates the probability that the motif starts at the specified location in a sequence given a PWM.

    Parameters
    ----------
    seq : ndarray
        Sequence to be considered.
    pwm : ndarray
        Position-Weighted Matrix to use for calculations.
    startLoc : int
        Starting location of the motif.

    Returns
    -------
    prob : the probability as a single float.
    """
    prob = 0

    for i in range(1, pwm.shape[1]):
        prob += np.log(pwm[seq[startLoc + i - 1], i] / pwm[seq[startLoc + i - 1], 0])

    return prob

def sampling_step(seqs, pwm, locs, locIdx):
    """
    For the removed sequence, calculate the distribution array
    and choose a new location (in-place) based on the new distribution array.

    Parameters
    ----------
    seqs : ndarray
        Array representing the set of sequences.
    pwm : ndarray
        Position-Weighted Matrix to use for calculations.
    locs : ndarray
        Locations of all the predicted motif locations.
    locIdx : int
        Location of the removed sequence to update.
    """
    distribution = getDistributionForSeq(seqs[locIdx, :], pwm)
    locs[locIdx] = np.random.choice(distribution.shape[0], 1, p=distribution)

def getConsensus(kmers):
    """
    Gets the consensus string/array from a set of kmers.

    Parameters
    ----------
    kmers : ndarray
        Array representing the set of kmers.

    Returns
    -------
    consensus : ndarray of ints representing the common consensus.
    """
    #count up the char occurances for each pos
    counts = np.zeros((4, kmers.shape[1]))
    for i in range(4):
        counts[i, :] = np.count_nonzero(kmers == i, axis=0)

    #find which char has the highest count for each pos
    consensus = np.argmax(counts, axis=0)
    
    return consensus