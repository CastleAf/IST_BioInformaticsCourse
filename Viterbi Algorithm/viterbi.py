import sys
import math
import random
import numpy as np

# TODO: Insert pretty HMM Model graph on github
# TODO: Explain in readMe the -inf values.
# TODO: Explain in readMe that this algorithm is bruteforce coded
# to a HMM Model, serves more as an example and not the main algorithm
# (can't be applied to diferent HMM)

# FIXME: Viterbi does not generate the optimal path, it generates the most likely
# path for a certain input. Change vars names

# Function that does a traceback on the viterbi matrix to find the optimal path 
def doTraceback(matrix, prevCellsIndex):

    optimalPath = []
    columnNumber = len(matrix[0]) - 1

    # Firstly, we'll pick the state with the maximum value in the last column 
    lastStates = matrix[1:, columnNumber]
    bestValue = max(lastStates[0], lastStates[1], lastStates[2])
    bestIndex = getBestValIndex(bestValue, lastStates[0], lastStates[1], lastStates[2])
    optimalPath.append(bestIndex)

    myState = bestIndex

    # Tracebacks following prevCellsIndex to get the optimal path
    while columnNumber > 1:

        newState = int(prevCellsIndex[myState, columnNumber])
        optimalPath.append(newState)
        myState = newState
        columnNumber -= 1

    # Path is backwards, .reverse() method reverses it.
    optimalPath.reverse()
    return optimalPath


# Function that returns the max value's index (+ 1) out of 3 values. 
# In case some of them are equal, returns one of them randomly
def getBestValIndex(value, firstVal, secondVal, thirdVal):

    # Pick a random state in case more than one path is optimal
    if (value == firstVal and value == secondVal and value == thirdVal):
        return random.choice([1, 2, 3])

    if (value == firstVal and value == secondVal):
        return random.choice([1, 2])
    
    if (value == firstVal and value == thirdVal):
        return random.choice([1, 3])
    
    if (value == secondVal and value == thirdVal):
        return random.choice([2, 3])

    # Picks only the optimal path in case there aren't more than one
    if (value == firstVal):
        return 1

    if (value == secondVal):
        return 2

    if (value == thirdVal):
        return 3
    
    # Shouldn't get here
    print("Error: Unable to get the best value's Index. Exiting.")
    raise SystemExit


# Function to avoid doing logarithm of 0 values. log(0) can happen 
# when a certain state does not emit a symbol or when we 
# test a state transition that does not happen (probability of 0).
def doLog(value):
    
    # In case the value to log is 0, returns an infinite
    # negative value, so that that option is never picked
    if (value == 0):
        return -float("inf")
    
    return math.log(value)


# Function to return the index of symbols in the alphabet
def getSymb(sequence, i):

    symb = sequence[i]

    if (symb == 'A'):
        return 0
    
    if (symb == 'C'):
        return 1
    
    if (symb == 'G'):
        return 2
    
    if (symb == 'T'):
        return 3

    # Shouldn't get here
    print("Error: Symbol Index not identified in sequence. Exiting.")
    raise SystemExit


# Function that will create a matrix with vk(i) values.
# Returns that matrix and prevIndex, a similar matrix
# where each cell represents the state from 
# which that cell was generated
def calcViterbiMatrix(seq, transProb, emissProb):
    
    j = 1
    seqSize = len(seq)

    # 1. Initialization
    myMatrix = np.zeros((4, seqSize + 1))
    myMatrix[0, 0] = 1
    prevIndex = np.zeros((4, seqSize + 1))

    # 2. Recursion
    while j < (seqSize + 1):

        symbIndex = getSymb(seq, j - 1)
    
        # On the first column, all values will generate from state 0
        if (j == 1):

            myMatrix[1, 1] = doLog(emissProb[1, symbIndex]) + 1 + doLog(transProb[0, 1])
            myMatrix[2, 1] = doLog(emissProb[2, symbIndex]) + 1 + doLog(transProb[0, 2])
            myMatrix[3, 1] = doLog(emissProb[3, symbIndex]) + 1 + doLog(transProb[0, 3])
            j += 1
            continue
        
        for k in range(1, 4):
            
            firstStateVal = doLog(emissProb[k, symbIndex]) + myMatrix[1, j - 1] + doLog(transProb[1, k])
            secondStateVal = doLog(emissProb[k, symbIndex]) + myMatrix[2, j - 1] + doLog(transProb[2, k])         
            thirdStateVal = doLog(emissProb[k, symbIndex]) + myMatrix[3, j - 1] + doLog(transProb[3, k])
            
            bestVal = max(firstStateVal, secondStateVal, thirdStateVal)

            # Fill prevIndex.
            # prevIndex[k, j] = 1: The similar cell in myMatrix was generated from state 1
            # prevIndex[k, j] = 2: The similar cell in myMatrix was generated from state 2
            # prevIndex[k, j] = 3: The similar cell in myMatrix was generated from state 3
            prevIndex[k, j] = getBestValIndex(bestVal, firstStateVal, secondStateVal, thirdStateVal)
            myMatrix[k, j] = bestVal

        j += 1

    print(myMatrix)
    print(prevIndex)
    return myMatrix, prevIndex


# Main Program:
def main():

    # Receive user input. Upper method makes the string uppercase.
    sequence = input("Please input the DNA sequence: ").upper()
    alphabet = ['A', 'C', 'G', 'T']

    # Checks if input is valid
    for char in sequence:
        if char not in alphabet:
            print("Error: Input must be a DNA Sequence.")
            return
    
    # Transition probabilities from state [line] to state [column]
    transProb = np.array([[0, 1/3, 1/3, 1/3],[0, 0.6, 0.4, 0],[0, 0.25, 0.5, 0.25],[0, 0.25, 0.25, 0.5]])

    # Lines -> State Numbers, Columns -> Emissions (0 = 'A', 1 = 'C', 2 = 'G', 3 = 'T')
    emissProb = np.array([[0, 0, 0, 0],[0.4, 0, 0.3, 0.3],[0.1, 0.4, 0.4, 0.1],[0.4, 0.3, 0, 0.3]])

    # solvedMatrix   -> Matrix with vl(i) values calculated
    # prevCellsIndex -> Indexes from what line were the values were generated
    solvedMatrix, prevCellsIndex = calcViterbiMatrix(sequence, transProb, emissProb)

    # Calculate optimal path doing a traceback on the viterbi matrix
    optimalPath = doTraceback(solvedMatrix, prevCellsIndex)

    print(optimalPath)



# Runs main program if this module isn't being imported
if __name__ == '__main__':
    sys.exit(main())