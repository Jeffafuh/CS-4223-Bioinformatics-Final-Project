import numpy as np

"""
    Returns the indices of the maximum values along an axis.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int, optional
        By default, the index is into the flattened array, otherwise
        along the specified axis.
    out : array, optional
        If provided, the result will be inserted into this array. It should
        be of the appropriate shape and dtype.

    Returns
    -------
    index_array : ndarray of ints
"""
def runGibbs(seqs, motifL, minIters, maxIters):
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


#Initialize the PSSM matrix to the appropriate shape
def initPWMMatrix(motifL):
    return np.zeros((4, motifL+1))

#Init weights to a random sampling of the sequences
def initRandomSamplePos(sequences, motifL):
    locs = np.random.randint(0, sequences.shape[1] - motifL+1, (sequences.shape[0],))

    return locs

def getKmers(seqs, locs, motifL):
    kmers = np.zeros((seqs.shape[0], motifL))
    for x in range(locs.shape[0]):
        kmers[x, :] = seqs[x, locs[x]:locs[x]+motifL]
    
    return kmers

def getBackground(seqs, locs, motifL):
    background = np.zeros((seqs.shape[0], seqs.shape[1] - motifL))
    for i in range(seqs.shape[0]):
        background[i, :] = np.delete(seqs[i, :], np.arange(locs[i], locs[i]+motifL))
    return background
    

#Given a sequence to remove, remove the sequence from seqs then recalculate the pwm matrix
def predict_step(seqs, pwm, locs, idxToRemove):
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

    return idxToRemove

#Given one sequence, a location in that sequence, and a pwm, calculate the log likelyhood that this pwm occurs at the startLoc
def calcLikelyhoodOfSeq(seq, pwm, startLoc):
    prob = 0

    for i in range(1, pwm.shape[1]):
        prob += np.log(pwm[seq[startLoc + i - 1], i] / pwm[seq[startLoc + i - 1], 0])

    return prob

#Given one sequence, calculate the likelyhood of the pwm occuring at each possible loc in the sequence
def getDistributionForSeq(seq, pwm):
    distribution = np.zeros((seq.shape[0] - pwm.shape[1] + 2,))

    for i in range(distribution.shape[0]):
        distribution[i] = np.exp(calcLikelyhoodOfSeq(seq, pwm, i))
    
    #normalize
    distribution /= np.sum(distribution)

    return distribution

#For the removed sequence, calculate a distribution matrix and choose a new location based on the new distribution returned
def sampling_step(seqs, pwm, locs, locIdx):
    distribution = getDistributionForSeq(seqs[locIdx, :], pwm)
    locs[locIdx] = np.random.choice(distribution.shape[0], 1, p=distribution)

#Given a set of kmers, find the consensus among them
def getConsensus(kmers):
    #count up the char occurances for each pos
    counts = np.zeros((4, kmers.shape[1]))
    for i in range(4):
        counts[i, :] = np.count_nonzero(kmers == i, axis=0)

    #find which char has the highest count for each pos
    consensus = np.argmax(counts, axis=0)
    
    return consensus