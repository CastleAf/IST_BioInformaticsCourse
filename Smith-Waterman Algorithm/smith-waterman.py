# Bioinformatics / Computational Biology Course. 1st Semester 2020/2021.
# Usage: python3 smith-waterman.py

import sys
import itertools
import numpy as np
import time


# Function to load the scoring model matrix
def createMatrix(file):
    txtContents = open(file).read()
    return [item.split() for item in txtContents.split('\n')[:-1]]


# Get Blossom50 values
def findScore(upperLetter, leftLetter, scoreModel):

    if (upperLetter in scoreModel[0]):
        columnIndex = scoreModel[0].index(upperLetter)

    else:
        if (upperLetter == leftLetter):
            return 1
        else: return -5

    if (leftLetter in scoreModel[0]):
        rowIndex = scoreModel[0].index(leftLetter)
    else: return -5

    return int(scoreModel[rowIndex][columnIndex])


# Creates the alignment matrix with scores
def createScoringMatrix(a, b, scoreModel, gapPenalty):

    numRows = len(b) + 1
    numColumns = len(a) + 1
    bestAlignVal = 0
    bestAlignsIndex = []
    
    # Initialize matrix with int 0 values
    myMatrix = np.zeros((numRows, numColumns), dtype = int)

    # prevCells[i][j] will track from which direction was the score at i,j generated
    prevCells = np.zeros((numRows, numColumns), dtype = int)

    for i, j in itertools.product(range(1, numRows), range(1, numColumns)):

        scoreModelVal = myMatrix[i - 1, j - 1] + findScore(a[j - 1], b[i - 1], scoreModel)
        upperGapVal = myMatrix[i - 1, j] + gapPenalty
        lowerGapVal = myMatrix[i, j - 1] + gapPenalty

        bestVal = max(0, scoreModelVal, upperGapVal, lowerGapVal)

        # Fill prevCells matrix
        if (bestVal == scoreModelVal):
            prevCells[i, j] = 1 # Cell score was generated from diagonal value
        if (bestVal == upperGapVal):
            prevCells[i, j] = 2 # Cell score was generated from an upper gap    
        if (bestVal == lowerGapVal):
            prevCells[i, j] = 3 # Cell score was generated from a lower gap

        if (bestVal == scoreModelVal and bestVal == upperGapVal):
            prevCells[i, j] = 4 # Cell score can be generated from an upper gap or a diagonal value
        if (bestVal == scoreModelVal and bestVal == lowerGapVal):
            prevCells[i, j] = 5 # Cell score can be generated from a lower gap or a diagonal value 

        if (bestVal == 0):
            prevCells[i, j] = 0 # Cell score was negative, cell will have a 0 val
    
        # Assign the score to the current cell
        myMatrix[i, j] = bestVal

        # Keep track of the best aligns
        if ((bestVal != 0) and (bestVal == bestAlignVal)):
            bestAlignsIndex.append([i, j])
        if (bestVal > bestAlignVal):
            bestAlignVal = bestVal
            bestAlignsIndex.clear() # Clears older best aligns indexes
            bestAlignsIndex.append([i, j])

    return bestAlignVal, bestAlignsIndex, myMatrix, prevCells


# Function to traceback the matrix and get an align
def calcPaths(nextStep, a, b, row, column, prevCellsMat):

    top = []
    bot = []
    
    # Flag to mark if there are secundary paths to handle
    hasSecundaryPath = False

    # Assign names to number values
    DIAG, UP, BOT, DIAG_UP, DIAG_BOT = range(1, 6)

    while(nextStep):

        if (nextStep == DIAG): # Diagonal path
            top.append(a[column - 1])
            bot.append(b[row - 1])
            row -= 1
            column -= 1

        if (nextStep == UP): # Upper gap path
            top.append('-')
            bot.append(b[row - 1])
            row -= 1

        if (nextStep == BOT): # Lower gap path
            top.append(a[column - 1])
            bot.append('-')
            column -= 1

        if (nextStep == DIAG_UP): # Diagonal and upper gap paths

            # Will do the upper gap path next
            hasSecundaryPath = True
            prevCellsMat[row, column] = UP

            # Diagonal path
            top.append(a[column - 1])
            bot.append(b[row - 1])
            row -= 1
            column -= 1

        if (nextStep == DIAG_BOT): # Diagonal and lower gap paths
            
            # Will do the lower gap path next
            hasSecundaryPath = True
            prevCellsMat[row, column] = BOT

            top.append(a[column - 1])
            bot.append(b[row - 1])
            row -= 1
            column -= 1
        
        nextStep = prevCellsMat[row, column]
    
    return top, bot, hasSecundaryPath
    

# Function to calculate all the best aligns
def calcAligns(bestAlignsIndex, prevCellsMatrix, seqOne, seqTwo):

    # Lists to save the upper and lower alignments
    topAlignments = []
    botAlignments = []

    # Duplicates prevCellsMatrix to reset the matrix to the original
    # values later, so we can find all path combinations
    originalPrevMatrix = prevCellsMatrix.copy()
    doReset = 0

    for indexes in bestAlignsIndex:

        row = indexes[0]
        column = indexes[1]
        next_step = prevCellsMatrix[row, column]

        top, bot, hasSecundaryPath = calcPaths(next_step, seqOne, seqTwo, row, column, prevCellsMatrix)

        # Append found alignments
        topAlignments.append(top)
        botAlignments.append(bot)
            
        # If at the second iteration, resets the prevCellsMatrix so it's possible to arrange all alignments combinations
        if (doReset):
            prevCellsMatrix = originalPrevMatrix
        
        while (hasSecundaryPath):

            top, bot, hasSecundaryPath = calcPaths(next_step, seqOne, seqTwo, row, column, prevCellsMatrix)
            
            # Condition to reject duplicates
            if (not (top in topAlignments and bot in botAlignments)):
                topAlignments.append(top)
                botAlignments.append(bot)

        # Marks first iteration completed
        doReset = 1

    return topAlignments, botAlignments


# Main Program:
def main():

    # Receive user input. Upper method makes the string uppercase.
    seqOne = str(input("First Amino Acid Sequence: ")).upper()
    seqTwo = str(input("Second Amino Acid Sequence: ")).upper()

    # Ideally a negative value
    gapPenalty = int(input("Gap Penalty Value: "))

    # Define Score Model, in this case we use the BLOSUM 50 score model.
    myScoreModel = createMatrix('blosum50.txt')

    # Saves timestamp to calculate the program's execution time
    start_time = time.time()

    # Creating the alignment matrix with scores, saves the matrix,
    # the previous cells matrix, the best aligns index and their score
    bestAlignScore, alignsIndex, alignmentMatrix, prevCellsMatrix = createScoringMatrix(seqOne, seqTwo, myScoreModel, gapPenalty)

    print()
    print("Here's the alignment score matrix: ")
    print()
    print(alignmentMatrix)

    topAlignments, botAlignments = calcAligns(alignsIndex, prevCellsMatrix, seqOne, seqTwo)

    # Printing the results
    nOptimalAligns = len(topAlignments)

    if (nOptimalAligns == 1):
        print()
        print("Found one optimal alignment with score " + str(bestAlignScore) + ":")
        print(''.join(reversed(topAlignments[0])))
        print(''.join(reversed(botAlignments[0])))

    if (nOptimalAligns > 1):
        print()
        print("Found " + str(len(topAlignments)) + " optimal alignments with score " + str(bestAlignScore) + ":")
        print()

        for n in range(len(topAlignments)):
            print("Alignment number " + str(n + 1) + ":")
            print(''.join(reversed(topAlignments[n])))
            print(''.join(reversed(botAlignments[n])))


    print()
    end_time = time.time()
    print("This algorithm took " + str(end_time - start_time) + " seconds to execute.")


# Runs main program if this module isn't being imported
if __name__ == '__main__':
    sys.exit(main())
