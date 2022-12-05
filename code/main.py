from MEMEFunctions import *
from GibbsFunctions import *
from testingFunctions import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pytictoc import TicToc
from random import sample

defaults = [20, 40, 7, 1]
params = defaults.copy()
params[3] = list(range(0, 8, 1))
generateInputFile('data/nMuts_Inputs.txt', params, 5, 0)

testInputFileOnMEME('data/nMuts_Inputs.txt', 'data/nMuts_MEME_Runtime.txt', 0)
testInputFileOnGibbs('data/nMuts_Inputs.txt', 'data/nMuts_Gibb_Runtime.txt', 0)
print('nMuts Done!')
