import numpy as np
from MEMEFunctions import *
from GibbsFunctions import *
from pytictoc import TicToc

# this function generates nSeqs of length seqL, and insert a kmer of
# length motifL into each sequence. Each kmer has exactly nMuts mutations from
# a common censensus sequence. 

# Returns the seqs as an nSeqs x seqL array of integers between 1 and 4,
# the consensus as an 1 x motifL array, and the locs where the kmers were
# inserted as a nSeqs x 1 array.
def generateSeqsWithMotif(nSeqs, seqL, motifL, nMuts):
    seqs = np.random.randint(0, 4, (nSeqs, seqL))
    consensus = np.random.randint(0, 4, (motifL,))
    locs = np.random.randint(0, seqL - motifL+1, (nSeqs,))
    charSet = np.arange(4)

    for i in range(nSeqs):
        pattern = consensus.copy()
        mutLocs = np.random.choice(motifL, size=(nMuts,), replace=False)
        for j in mutLocs:
            pool = np.setdiff1d(charSet, pattern[j])
            pattern[j] = np.random.choice(pool, 1)
        
        seqs[i, locs[i]:locs[i]+motifL] = pattern

    return (seqs, consensus, locs)

#Given two kmers, return the hamming distance between them
def kmerHammingDistance(kmer1, kmer2):
    return np.count_nonzero(kmer1 - kmer2)

#Given two arrays of locations, measure the difference in locations between them
def locHammingDistance(loc1, loc2):
    return np.sum(np.abs(loc1 - loc2))

def generateInputFile(fileOut, params, numTestsPerParam=1, randomState=None):
    """
    Generates an input file given a list of parameters to generate the sequences from

    Parameters
    ----------
    fileOut : string
        Name of the file to be written to (original file will first be truncated)
    params : list
        The list of parameters to generate the sequences from.
        The order of the parameters in the list are as follows:
        [nSeqs, seqL, motifL, nMuts]
        Whichever parameter is going to be changing will be in the form of a list in the respective index.
    numTestsPerParam : int, optional
        Determines how many sequences to generate per unique parameter combination
    randomState : float, optional
        If provided, sets the np.random.seed() to the provided value

    File Outout
    ----------
    Number_of_different_parameters_to_test numTests
    nSeqs seqL motifL nMuts

    <Flattened array of the sequences generated for test 1>

    <Consensus for the sequences generated for test 1>

    <Locations of the planted motifs for test 1>

    ...

    <Flattened array of the sequences generated for test numTests>

    <Consensus for the sequences generated for test NumTests>

    <Locations of the planted motifs for test NumTests>

    nSeqs seqL motifL nMuts numTests

    ...
    """
    if randomState != None:
        np.random.seed(randomState)

    iterIdx = 0
    for idx, param in enumerate(params):
        if isinstance(param, list):
            iterIdx = idx
            break
    
    def writeListToFile(file, arr):
        for i in arr:
            file.write(str(i)+' ')
        file.write('\n')
    
    def writeNDarrayToFile(file, arr, dim):
        arr = arr.tolist()
        if dim == 1:            
            writeListToFile(file, arr)
        else: #dim == 2
            for r in arr:
                writeListToFile(file, r)
    
    writeParams = params.copy()
    with open(fileOut, 'w') as file:
        file.write(str(len(params[iterIdx]))+' '+str(numTestsPerParam)+'\n')
        for param in params[iterIdx]:
            writeParams[iterIdx] = param
            writeListToFile(file, writeParams)
            for _ in range(numTestsPerParam):
                (seqs, consensus, locs) = generateSeqsWithMotif(writeParams[0], writeParams[1], writeParams[2], writeParams[3])
                writeNDarrayToFile(file, seqs, 2)
                writeNDarrayToFile(file, consensus, 1)
                writeNDarrayToFile(file, locs, 1)

def testInputFileOnMEME(fileIn, fileOut, randomState=None, verbose=True):
    if randomState != None:
        np.random.seed(randomState)

    timer = TicToc()
    with open(fileIn, 'r') as file:
        numParams, numTests = file.readline().strip().split()
        numParams, numTests = int(numParams), int(numTests)
        runtimes = np.zeros((numParams, numTests))
        
        for i in range(numParams):
            nSeqs, seqL, motifL, _= file.readline().strip().split()
            nSeqs, seqL, motifL = int(nSeqs), int(seqL), int(motifL)

            for j in range(numTests):
                seqs = []
                for _ in range(nSeqs):
                    row = file.readline().strip().split()
                    seqs.append(list(map(int, row)))
                seqs = np.array(seqs)

                file.readline()
                file.readline()

                timer.tic()
                runBasicMEME(seqs, motifL, 0.7, 0.001, 0.1)
                runtimes[i, j] = timer.tocvalue()

                if verbose:
                    print('Finished with test', j, 'for param', i)
            if verbose:
                print('Done with param', i)       

        runtimes = runtimes.mean(axis=1)
        np.savetxt(fileOut, runtimes)

def testInputFileOnGibbs(fileIn, fileOut, randomState=None, verbose=True):
    if randomState != None:
        np.random.seed(randomState)

    timer = TicToc()
    with open(fileIn, 'r') as file:
        numParams, numTests = file.readline().strip().split()
        numParams, numTests = int(numParams), int(numTests)
        runtimes = np.zeros((2, numParams, numTests))
        
        for i in range(numParams):
            nSeqs, seqL, motifL, _= file.readline().strip().split()
            nSeqs, seqL, motifL = int(nSeqs), int(seqL), int(motifL)

            for j in range(numTests):
                seqs = []
                for _ in range(nSeqs):
                    row = file.readline().strip().split()
                    seqs.append(list(map(int, row)))
                seqs = np.array(seqs)

                file.readline()
                file.readline()

                timer.tic()
                _, _, iters = runGibbs(seqs, motifL, 5, 300)
                runtimes[0, i, j] = timer.tocvalue()
                runtimes[1, i, j] = iters

                if verbose:
                    print('Finished with test', j, 'for param', i)
            if verbose:
                print('Done with param', i)   

        runtimes = runtimes.mean(axis=2)
        np.savetxt(fileOut, runtimes)