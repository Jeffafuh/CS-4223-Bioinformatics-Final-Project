# CS-4223-Bioinformatics-Final-Project

This GitHub repo is part of the final project assigned for the class CS-4223 (Bioinformatics) at UTSA during the fall 2022 semester taught by Dr. Jianhua Ruan.

The written report for this project can be found in the repo as "CS_4223_Project_Milestone_2.pdf"

## Version Details
All of the code for this project was written in Python 3.8.10 using VSCode. Additional libraries include:
- numpy 1.21.2 (used for the algorithm code)
- matplotlib 3.4.3 (used for plotting)
- pytictoc 1.5.2 (used to measure runtime)

##Installation
All code needed to just run MEME and Gibbs Sampling can be found in the files MEMEFunctions.py and GibbsFunctions.py, respectively. Additionally, all of the code used to test the runtimes of the algorithms can be found in their respective files. Proper documentation for all of the functions are included in the code as well. Once the appropriate file(s) are downloaded, simply place them in the appropriate project directory and include them in the python files as needed.

##Runtime testing
In order to create reproducable results, I stored all of the sequences used to test the algorithms in plain-text files should you wish to use them for your own tests located under `code/data`. Documentation and the functions needed to generate new input files or use those input files to test the runtime can be found in testingFunctions.py, createInputFiles.py, and createRuntimeFiles.py, respectively. Additionally, the format of these files can also be found in testingFunctions.py.
