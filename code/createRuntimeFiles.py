from testingFunctions import *

if __name__ == "__main__":
    """
    Makes a bunch of calls to testInputFileOnX() with the appropriate file names.
    Random state is set to 0 for reproducability.
    """

    #nSeqs
    testInputFileOnMEME('data/nSeqs_Inputs.txt', 'data/nSeqs_MEME_Runtime.txt', 0)
    testInputFileOnGibbs('data/nSeqs_Inputs.txt', 'data/nSeqs_Gibb_Runtime.txt', 0)
    print('nSeqs Done!')

    #seqL
    testInputFileOnMEME('data/seqL_Inputs.txt', 'data/seqL_MEME_Runtime.txt', 0)
    testInputFileOnGibbs('data/seqL_Inputs.txt', 'data/seqL_Gibb_Runtime.txt', 0)
    print('seqL Done!')

    #motifL
    testInputFileOnMEME('data/motifL_Inputs.txt', 'data/motifL_MEME_Runtime.txt', 0)
    testInputFileOnGibbs('data/motifL_Inputs.txt', 'data/motifL_Gibb_Runtime.txt', 0)
    print('motifL Done!')

    #nMuts
    testInputFileOnMEME('data/nMuts_Inputs.txt', 'data/nMuts_MEME_Runtime.txt', 0)
    testInputFileOnGibbs('data/nMuts_Inputs.txt', 'data/nMuts_Gibb_Runtime.txt', 0)
    print('nMuts Done!')