import sys
import itertools
import numpy as np

# TODO: Verificar que funciona com seq2 > seq1



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
    
    # Initialize matrix with int 0 values
    myMatrix = np.zeros((numRows, numColumns), dtype = int)

    for i, j in itertools.product(range(1, numRows), range(1, numColumns)):

        scoreModelVal = myMatrix[i - 1, j - 1] + findScore(a[j - 1], b[i - 1], scoreModel)
        upperGapVal = myMatrix[i - 1, j] - gapPenalty
        lowerGapVal = myMatrix[i, j - 1] - gapPenalty

        myMatrix[i, j] = max(0, scoreModelVal, upperGapVal, lowerGapVal)

    return myMatrix


# Main Function
def main():

    # Receive user input. Upper method makes the string uppercase.
    seqOne = str(input("First Amino Acid Sequence: ")).upper()
    seqTwo = str(input("Second Amino Acid Sequence: ")).upper()
    gapPenalty = int(input("Gap Penalty Value: "))

    # Define Score Model, in this case we use the BLOSUM 50 score model
    myScoreModel = createMatrix('blosum50.txt')

    
    alignmentMatrix = createScoringMatrix(seqOne, seqTwo, myScoreModel, gapPenalty)

    print(alignmentMatrix)


# Runs main program if this module isn't being imported
if __name__ == '__main__':
    sys.exit(main())
