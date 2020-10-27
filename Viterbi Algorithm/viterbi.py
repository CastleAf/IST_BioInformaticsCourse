import sys
import math
import random
import numpy as np

# TODO: Meter grafico bonito das transicoes como imagem no github

def doTraceback(matrix, prevCellsIndex):

    path = []
    columnNumber = len(matrix[0]) - 1

    # Firstly, we'll pick the state with the maximum value in the last column 
    lastStates = matrix[1:, columnNumber]
    bestValue = max(lastStates[0], lastStates[1], lastStates[2])
    bestIndex = getBestValIndex(bestValue, lastStates[0], lastStates[1], lastStates[2])
    path.append(bestIndex)

    myState = bestIndex

    # Tracebacks following prevCellsIndex to get the optimal path
    while columnNumber > 1:

        newState = int(prevCellsIndex[myState, columnNumber])
        print(newState)
        path.append(newState)
        myState = newState
        columnNumber -= 1

    # Path is backwards. This method reverses it.
    path.reverse()

    return path


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
    print("Error: Unable to get the Best Value Index. Exiting.")
    raise SystemExit

# Function to avoid doing log of 0 values
def doLog(value):
    
    # ln(0) is impossible, in case we have a 0 emission or
    # transition probability value, returns an infinite
    # negative value, so that option is never picked
    if (value == 0):
        return -float("inf")
    
    return math.log(value)

def getSymb(sequence, i):

    symb = sequence[i]

    if symb == 'A':
        return 0
    
    if symb == 'C':
        return 1
    
    if symb == 'G':
        return 2
    
    if symb == 'T':
        return 3

    # Shouldn't get here
    print("Error: Symbol Index not identified in sequence. Exiting.")
    raise SystemExit


def calcViterbiMatrix(seq, transProb, emissProb):
    
    j = 1
    seqSize = len(seq)

    # 1. Initialization
    matrix = np.zeros((4, seqSize + 1))
    matrix[0, 0] = 1
    prevIndex = np.zeros((4, seqSize + 1))

    # 2. Recursion

    while j < (seqSize + 1):

        symbIndex = getSymb(seq, j - 1)
    
        # On the first column, all values will generate from state 0
        if (j == 1):

            matrix[1, 1] = doLog(emissProb[1, symbIndex]) + 1 + doLog(transProb[0, 1])
            matrix[2, 1] = doLog(emissProb[2, symbIndex]) + 1 + doLog(transProb[0, 2])
            matrix[3, 1] = doLog(emissProb[3, symbIndex]) + 1 + doLog(transProb[0, 3])
            prevIndex[1, 1] = 0
            prevIndex[2, 1] = 0
            prevIndex[3 ,1] = 0
            j += 1
            continue
        
        for k in range(1, 4):
            
            firstStateVal = doLog(emissProb[k, symbIndex]) + matrix[1, j - 1] + doLog(transProb[1, k])
            secondStateVal = doLog(emissProb[k, symbIndex]) + matrix[2, j - 1] + doLog(transProb[2, k])         
            thirdStateVal = doLog(emissProb[k, symbIndex]) + matrix[3, j - 1] + doLog(transProb[3, k])
            
            bestVal = max(firstStateVal, secondStateVal, thirdStateVal)
            prevIndex[k, j] = getBestValIndex(bestVal, firstStateVal, secondStateVal, thirdStateVal)

            matrix[k, j] = bestVal

        j += 1

    print(matrix)
    print(prevIndex)
    return matrix, prevIndex


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

    # solvedMatrix -> Matrix with vl(i) values calculated
    # prevCellsIndex -> Indexes from what line were the values were generated
    solvedMatrix, prevCellsIndex = calcViterbiMatrix(sequence, transProb, emissProb)

    optimalPath = doTraceback(solvedMatrix, prevCellsIndex)

    print(optimalPath)



# Runs main program if this module isn't being imported
if __name__ == '__main__':
    sys.exit(main())