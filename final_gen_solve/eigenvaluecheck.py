import EigensetReader
import numpy as np
import math

def verify(file_str):
    e = EigensetReader.Eigenset()
    e.read(file_str)
    num_eigenvectors = e.numberEigenvectors
    print("Read " + repr(num_eigenvectors) + " eigenvectors")
    
    for i in range (0, num_eigenvectors):
        currentPair = e.eigenpairs[i]
        print("Bx = " + repr(currentPair.Bx) + ", Bz = " + repr(currentPair.Bz) + ", J = " + repr(currentPair.J) +  ", Eigenvalue = " + repr(currentPair.Eigenvalue))
        for j in range(0, e.eigenvectorSize):  
            print(repr(currentPair.Eigenvector[j]))
        print()


def main():
    filename = "results.eigenset"
    
    x = verify(filename)
    

    


main()


