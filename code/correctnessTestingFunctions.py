from MEMEFunctions import *
from GibbsFunctions import *
from testingFunctions import *

def test_e_m_steps(nSeqs, seqL, motifL):
    (seqs, consensus, locs) = generateSeqsWithMotif(nSeqs, seqL, motifL, 0)

    print(seqs)
    weights = memeInitPWM(seqs, motifL, 0.7)
    z = np.zeros((seqs.shape[0], seqs.shape[1]-weights.shape[1]+2))
    print(weights)
    e_step(seqs, weights, z)
    print(z)
    print(z.sum(axis=1))
    m_step(seqs, weights, z)
    print(np.sum(weights, axis=0))
    print(weights)

#a=0, c=1, g=2, t=3

def testMEMEToy():
    seqs = np.array([[2, 3, 1, 0, 2, 2],
                 [2, 0, 2, 0, 2, 3],
                 [0, 1, 2, 2, 0, 2],
                 [1, 1, 0, 2, 3, 1]])

    pwm, z, prob = runBasicMEME(seqs, 3, 0.7, 0.001)
    print(getKmers(seqs, getZLocs(z), 3))
    #print(consensus)
    print(getZLocs(z))
    #print(locs)
    #print(getConsensus(pwm))
    #print(consensus)
    print(prob)

def testGibbs(nSeqs, seqL, motifL, nMuts):
    """ nSeqs = 20
    seqL = 500
    motifL = 15
    nMuts = 6 """
    (seqs, consensus, locs) = generateSeqsWithMotif(nSeqs, seqL, motifL, nMuts)

    pssm, newLocs = runGibbs(seqs, motifL, 5, 300)
    newConsensus = getConsensus(getKmers(seqs, newLocs, motifL))
    print('Found Motif: ', newConsensus)
    print('Found Locs', newLocs)

    print('Planted Motif', consensus)
    print('Planted Locs', locs)

    print('loc hamming dist: ', locHammingDistance(newLocs, locs))
    print('consensus hamming dist: ', kmerHammingDistance(newConsensus, consensus))