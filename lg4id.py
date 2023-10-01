import re
import numpy as np

def getTestString(filename="testG4sequences.txt"):
    #print(filename)
    lineList = []
    with open(filename, 'r') as file:
        for line in file:
            lineList.append(line.strip())
    #re.sub(r"\r?\n", "", outString)
    return "".join(lineList)


def getMatches(inString):
    stringLength = len(inString)

    patternLength = 3
    correctedLength = stringLength - patternLength + 1
    boolArray = np.zeros( (correctedLength , 2) )
    
    stringIterator = iter(range(correctedLength))
    for i in stringIterator:

        tempString = inString[i:i + patternLength]
        if (tempString.upper() == "CCC"):
            boolArray[i, 0] = 1

            for _ in range(patternLength - 1):
                try:
                    next(stringIterator)

                except StopIteration:
                    continue
        elif (tempString.upper() == "GGG"):
            boolArray[i, 1] = 1

            for _ in range(patternLength - 1):
                try:
                    next(stringIterator)
                    #print(i)
                except StopIteration:
                    continue

    return boolArray

def processDictEntry(lG4List, currentStart, currentEnd, inString):
    if(len(lG4List) > 0):
        #get previous entry
        lastEntry = lG4List[-1]
        #go back to zero index
        lastStart = lastEntry['binStart'] - 1
        lastEnd = lastEntry['binEnd']
        
        #need to merge
        if(currentStart <= lastEnd):
            #remove the last entry
            lG4List.pop()
            currentStart = lastStart
            #current End stays the same
    
    sequenceString = inString[currentStart:currentEnd ]
    tempDict = {'binStart': currentStart + 1, 'binEnd': currentEnd, 'sequence' : sequenceString }
    lG4List.append(tempDict)

def driver(inString, windowSize=500, mustExceed=40):
    stringLength = len(inString)

    lG4List = []
    inLG4 = False
    currentStart = 0
    currentEnd = 0
    runningTotalGGG = 0
    runningTotalCCC = 0
    
    matchArray = getMatches(inString)
    
    #calculatate sums for first bin to initialize
    for i in range(windowSize):
        tempCCC = matchArray[i, 0]
        tempGGG = matchArray[i, 1]
        
        runningTotalCCC += tempCCC
        runningTotalGGG += tempGGG
        
    maxRunningTotal = max(runningTotalCCC, runningTotalGGG)
    if (maxRunningTotal > mustExceed):
        #this means the first bin was an LG4
        inLG4 = True
        currentStart = 0
        currentEnd = windowSize
        
        
    #now go through next bin through end, adding next entry to running total and subtracting last entry of last bin to total
    for i in range(windowSize, len(matchArray)):

        lastBinStart = i - windowSize
        
        tempCCC = matchArray[i, 0]
        tempGGG = matchArray[i, 1]
 
        tempCCClastBin = matchArray[lastBinStart, 0]
        tempGGGlastBin = matchArray[lastBinStart, 1]
 
        runningTotalCCC += tempCCC
        runningTotalGGG += tempGGG
        
        runningTotalCCC -= tempCCClastBin
        runningTotalGGG -= tempGGGlastBin

        maxRunningTotal = max(runningTotalCCC, runningTotalGGG)
        #print(maxCount)
        if (maxRunningTotal > mustExceed):
            if ( inLG4 == False ):
                currentStart = i - windowSize + 1
                currentEnd = i + 1
                inLG4 = True
            else:
                currentEnd = i + 1
        else:
            if ( inLG4 == True ):
                inLG4 = False
                #sequence string applies to previous bin (because current failed), so - 1
                
                #check that not overlapping with previous entry
                #can modify lG4List as side effect, adding or removing entries
                processDictEntry(lG4List, currentStart, currentEnd, inString)
                #sequenceString = inString[currentStart:currentEnd ]
                #tempDict = {'binStart': currentStart + 1, 'binEnd': currentEnd, 'sequence' : sequenceString }
                #lG4List.append(tempDict)
    #on last iteration, check if in LG4 to write outString
    if ( inLG4 == True ):
        inLG4 = False
        #sequence string applies to previous bin (because current failed), so - 1
        processDictEntry(lG4List, currentStart, currentEnd, inString)
    return lG4List

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
        prog="LG4ID",
        description="Long-G4 (Guanine Quadruplex) Identifier")
    parser.add_argument('-i', '--infile', required=True, help="Input Fasta File")
    parser.add_argument('-o', '--outfile', default="out.csv", required=False, help="Output File: default 'out.csv'")
    parser.add_argument('-w', '--window', default=1500, help="Window Size (default 1500)")
    parser.add_argument('-m', '--min', default=120, help="Must exceed/minimum hits per window (default 120)")
    
    
    args = parser.parse_args()
    
    
    testString = getTestString(args.infile)

    driverOutList = driver(testString, windowSize=1500, mustExceed=120)
    
    import csv
    with open(args.outfile, 'w', newline='') as csvfile:
        fieldnames = ['binStart', 'binEnd', 'sequence']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in driverOutList:
            print(entry)
            writer.writerow(entry)
            
