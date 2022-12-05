from testingFunctions import *

if __name__ == "__main__":
    """
    Makes a bunch of calls to generateInputFile() with the appropriate file names and parameter settings
    for nSeqs, seqL, motifL, and nMuts individually.
    Random state is set to 0 for reproducability.
    """

    #Default values
    #[nSeqs, seqL, motifL, nMuts]
    defaults = [20, 40, 7, 1]
    #nSeqs
    params = defaults.copy()
    params[0] = list(range(10, 51, 2))
    generateInputFile('data/nSeqs_Inputs.txt', params, 5, 0)

    #seqL
    params = defaults.copy()
    params[1] = list(range(10, 51, 2))
    generateInputFile('data/seqL_Inputs.txt', params, 5, 0)

    #motifL
    params = defaults.copy()
    params[2] = list(range(5, 21, 1))
    generateInputFile('data/motifL_Inputs.txt', params, 5, 0)

    #nMuts
    params = defaults.copy()
    params[3] = list(range(0, 7, 1))
    generateInputFile('data/nMuts_Inputs.txt', params, 5, 0)