import numpy as np
from MEMEFunctions import *
from GibbsFunctions import *
from pytictoc import TicToc

def kmerHammingDistance(kmer1, kmer2):
    """
    Calculates the hamming distance between 2 kmers.

    Parameters
    ----------
    kmer1 : ndarray
        First kmer to be considered.
    kmer2 : ndarray
        Second kmer to be considered.

    Returns
    ----------
    The hamming distance between kmer1 and kmer2.
    """
    return np.count_nonzero(kmer1 - kmer2)

#Given two arrays of locations, measure the difference in locations between them
def locHammingDistance(loc1, loc2):
    """
    Calculates the total point-wise difference between location ndarrays.

    Parameters
    ----------
    loc1 : ndarray
        First array to be considered.
    loc2 : ndarray
        Second array to be considered.

    Returns
    ----------
    The total point-wise difference between loc1 and loc2.
    """
    return np.sum(np.abs(loc1 - loc2))

def generateSeqsWithMotif(nSeqs, seqL, motifL, nMuts):
    """
    This function generates nSeqs of length seqL, and inserts a kmer of
    length motifL into each sequence. Each kmer has exactly nMuts mutations from
    a common censensus sequence. 

    Returns the sequences, motif, and motif locations planted within the sequences.

    Parameters
    ----------
    nSeqs : int
        Number of sequences to be generated.
    seqL : int
        Length of each sequence.
    motifL : int
        Length of the planted motif.
    nMuts : int
        Number of mutations for the motif.

    Returns
    -------
    seqs : ndarray of shape nSeqs x seqL with integers between 0 and 3, inclusive
    consensus : ndarray containing the original planted motif
    locs : ndarray containing the location of the planted motif for each sequence
    """
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

def generateInputFile(fileOut, params, numTestsPerParam=1, randomState=None):
    """
    Generates an input file given a list of parameters to generate the sequences from.

    Parameters
    ----------
    fileOut : string
        Name of the file to be written to (original file will first be truncated).
    params : list
        The list of parameters to generate the sequences from.
        The order of the parameters in the list are as follows:
        [nSeqs, seqL, motifL, nMuts]
        Whichever parameter is going to be changing will be in the form of a list in the respective index.
    numTestsPerParam : int, optional
        Determines how many sequences to generate per unique parameter combination.
    randomState : float, optional
        If provided, sets the np.random.seed() to the provided value.

    File Format
    ----------
    Number_of_different_parameters_to_test numTests
    nSeqs seqL motifL nMuts

    <Array of the sequences generated for test 1>

    <Consensus for the sequences generated for test 1>

    <Locations of the planted motifs for test 1>

    ...

    <Array of the sequences generated for test numTests>

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
    """
    Generates an output file of the runtimes for MEME
    given an input file to read the sequence sets from.
    The input file is in the format specified in generateInputFile().

    Parameters
    ----------
    fileIn : string
        Name of the file to read the sequences from.
    fileOut : string
        Name of the file to be written to (original file will first be truncated).
    randomState : float, optional
        If provided, sets the np.random.seed() to the provided value.
    verbose : boolean, optional
        Determines whether or not to print out messages to stdout after the completion of each test.

    File Output
    ----------
    The runtimes for the algorithm are written as an ndarray to a .txt file using np.savetxt().
    """
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
    """
    Generates an output file of the runtimes for Gibbs Sampling
    given an input file to read the sequence sets from.
    The input file is in the format specified in generateInputFile().

    Parameters
    ----------
    fileIn : string
        Name of the file to read the sequences from.
    fileOut : string
        Name of the file to be written to (original file will first be truncated).
    randomState : float, optional
        If provided, sets the np.random.seed() to the provided value.
    verbose : boolean, optional
        Determines whether or not to print out messages to stdout after the completion of each test.

    File Output
    ----------
    The runtimes for the algorithm are written as an ndarray to a .txt file using np.savetxt().
    """
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