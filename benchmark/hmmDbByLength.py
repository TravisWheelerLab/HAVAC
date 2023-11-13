import sys
import os

hmmLengthsFileSrc = "hmmLenghts.txt"


def writeUpToLength(rfamSourceFile, outputFile, desiredLength):
  currentLength = 0
  for line in rfamSourceFile.readlines():
    splitOnWhitespace = line.split()
    if splitOnWhitespace[0] == "HMMER3/f" and currentLength>= desiredLength:
      print("found models up to desired length {}. actual length of db is {}".format(desiredLength, currentLength))

      with open(hmmLengthsFileSrc, "a") as lengthsFile:
        lengthsFile.write("{}, ".format(currentLength))
      return
    elif splitOnWhitespace[0] == "LENG":
      length = splitOnWhitespace[1]
      currentLength += int(length)
    
    outputFile.write(line);

  print("could not generate up to desired length {}. reached Length {}".format(desiredLength, currentLength))
  with open(hmmLengthsFileSrc, "a") as lengthsFile:
    lengthsFile.write("{}, ".format(currentLength))


def main():

  if len(sys.argv) != 2:
    print("Error: requires 1 arg: file src to the full RFAM hmm file.")
    exit(-1)

  fullRfamHmmFileSrc = sys.argv[1]


  modelDbLengths = [1000,5000,10000,20000,30000,40000,50000,
                    60000,70000,80000,90000,100000,150000,]


  #write the beginning of the lengths array file
  with open(hmmLengthsFileSrc, "w") as lengthsFile:
    lengthsFile.write("lengths = [")

  for dbLength in modelDbLengths:
    outputFileSrc = "rfam_{}.hmm".format(dbLength)
    with open(fullRfamHmmFileSrc, "r") as rfamDbFile:
      with open(outputFileSrc, "w") as outputFile:
        writeUpToLength(rfamDbFile, outputFile, dbLength);
  
  with open(hmmLengthsFileSrc, "a") as lengthsFile:
    lengthsFile.write("]")

if __name__ == "__main__":
  main()