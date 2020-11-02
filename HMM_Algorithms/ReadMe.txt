BC-L3: 

This program corresponds to the implementation of Viterbi and Forward-Backward
algorithms for the given HMM model on Group III.a). These results are 
to be used on questions III.b), III.c) and III.d).

# How to run:
To run this program you only need to run the python file (HMM_Algorithms.py)
using python3. The program will ask you for an input sequence. The input
sequence must be composed of symbols that belong to the alphabet (ACTG, 
symbols that correspond to the nitrogen-containing nucleobases of DNA).

The program will output the optimal sequence of hidden states (pi*, using
the viterbi algorithm), a probability P(S) (using the forward algorithm) and
the posterior probability matrix (using the forward-backward algorithm).

# Run example

(Using the input S = ATCG)

$ python3 HMM_Algorithms.py
Please input the DNA sequence: ATCG

The most likely sequence of states for the given input is 1122.
This sequence has a probability of 0.003519333333333334.
The posterior probabilities matrix is:
 [[0.         0.         0.         0.         0.        ]
 [0.         0.40329608 0.45633643 0.         0.30488729]
 [0.         0.10170013 0.20002605 0.79363989 0.69511271]
 [0.         0.49500379 0.34363753 0.20636011 0.        ]]



# TODO: Insert pretty HMM Model graph on github
# TODO: Explain in readMe the -inf values.
# TODO: Explain in readMe that this algorithm is bruteforce coded
# to a HMM Model, serves more as an example and not the main algorithm
# (can't be applied to diferent HMM). Explain also why we use the log formula

