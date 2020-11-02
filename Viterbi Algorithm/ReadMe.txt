BC-L3: 

This program corresponds to the implementation of Viterbi and Forward 
algorithms for the given HMM model on Group III.a). These results are 
to be used on questions III.b) and III.c).

# How to run:
To run this program you only need to run the python file (viterby_forward.py)
using python3. The program will ask you for an input sequence. The input
sequence must be composed of symbols that belong to the alphabet (ACTG, 
symbols that correspond to the nitrogen-containing nucleobases of DNA).

The program will output the optimal sequence of hidden states (pi*, using
the viterbi algorithm) and a probability P(S) (using the forward algorithm).

# Run example

Using the input S = CATGCGGGTTATAAC:

$ python3 viterby_forward.py
Please input the DNA sequence: CATGCGGGTTATAAC

The most likely sequence of states for the given input is 211222221111112.
This sequence has a probability of 9.38645997061689e-10.

