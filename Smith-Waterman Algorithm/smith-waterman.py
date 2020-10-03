import sys
import itertools
import numpy as np

# TODO: Verificar que funciona com seq2 > seq1 e com dados com varios local alignments.
# TODO: Tratar de resultados de tempo

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
    # 1 = diagonal
    # 2 = up
    # 3 = left
    # 0 = placeholder
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
            prevCells[i, j] = 3 # Cell score was generated from a left gap
    
        # Assign the score to the current cell
        myMatrix[i, j] = bestVal


        if ((bestVal != 0) and (bestVal == bestAlignVal)):
            bestAlignsIndex.append([i, j])

        if (bestVal > bestAlignVal):
            bestAlignVal = bestVal
            bestAlignsIndex.clear() # Clears older best aligns indexes
            bestAlignsIndex.append([i, j])

    return bestAlignVal, bestAlignsIndex, myMatrix, prevCells


def calcAligns(bestAlignsIndex, alignmentMatrix, prevCellsMatrix, seqOne, seqTwo):

    # Lists to save the upper and lower alignments
    topAlignments = []
    botAlignments = []

    for indexes in bestAlignsIndex:

        top = []
        bot = []

        row = indexes[0]
        column = indexes[1]
        next_step = prevCellsMatrix[row, column]

        # Loop will stop next_step == 0
        while (next_step):
            
            if (next_step == 1):
                top.append(seqOne[column - 1])
                bot.append(seqTwo[row - 1])
                row -= 1
                column -= 1

            if (next_step == 2):
                top.append('-')
                bot.append(seqTwo[row - 1])
                row -= 1

            if (next_step == 3):
                top.append(seqOne[column - 1])
                bot.append('-')
                column -= 1

            next_step = prevCellsMatrix[row, column]

        topAlignments.append(top)
        botAlignments.append(bot)


    return topAlignments, botAlignments


# Main Function
def main():

    # Receive user input. Upper method makes the string uppercase.
    seqOne = str(input("First Amino Acid Sequence: ")).upper()
    seqTwo = str(input("Second Amino Acid Sequence: ")).upper()
    gapPenalty = int(input("Gap Penalty Value: "))

    # Define Score Model, in this case we use the BLOSUM 50 score model
    myScoreModel = createMatrix('blosum50.txt')

    # Creating the alignment matrix with scores, saves the matrix, the previous cells matrix, the best aligns index and their score
    bestAlignScore, alignsIndex, alignmentMatrix, prevCellsMatrix = createScoringMatrix(seqOne, seqTwo, myScoreModel, gapPenalty)

    print()
    print("Here's the alignment score matrix: ")
    print()
    print(alignmentMatrix)

    topAlignments, botAlignments = calcAligns(alignsIndex, alignmentMatrix, prevCellsMatrix, seqOne, seqTwo)

    nOptimalAligns = len(topAlignments)

    if (nOptimalAligns == 1):
        print()
        print("Found one optimal alignment with score " + str(bestAlignScore) + ":")
        print(''.join(reversed(topAlignments[0])))
        print(''.join(reversed(botAlignments[0])))

    if (nOptimalAligns > 1):
        print()
        print("Found " + str(len(topAlignments)) + " optimal alignments with score " + str(bestAlignScore) + ":")

        for n in range(len(topAlignments)):
            print("Alignment number " + str(n + 1) + ":")
            print(''.join(reversed(topAlignments[n])))
            print(''.join(reversed(botAlignments[n])))


# Runs main program if this module isn't being imported
if __name__ == '__main__':
    sys.exit(main())
