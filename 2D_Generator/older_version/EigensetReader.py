#Python reader
#Andy Grace, Josh Szachar 2021

import sys 
import struct
import copy

class IsingEigenset:
    J = None
    Bx = None
    Bz = None
    Eigenvalue = None
    Eigenvector = []

    def __init__(self, J, Bx, Bz, Eigenvalue):
        self.J = J
        self.Bx = Bx 
        self.Bz = Bz
        self.Eigenvalue = Eigenvalue

    def addEigenvectorValue(self, value):
        self.Eigenvector.append(copy.deepcopy(value))

    def clearEigenvector(self):
        self.Eigenvector.clear()

class Eigenset:
    version = 1
    isBigEndian = False
    eigenvectorSize = None
    numberEigenvectors = None
    eigenpairs = []

    def isSysBigEndian(self):
        if sys.byteorder == "little":
            self.isBigEndian = False
            return False
        else:
            self.isBigEndian = true
            return True

    def convertInt(self, integer, fromEndian, toEndian):
        if fromEndian != toEndian:
            return struct.unpack("<I", struct.pack(">I", integer))[0]
        return integer

    def convertDouble(self, double, fromEndian, toEndian):
        if fromEndian != toEndian:
            return struct.unpack("<d", struct.pack(">d", double))[0]
        return double

    def read(self, filename):
        bigEndianMode = False
        format = 0
        isBigEndianInt = 0

        file = open(filename, "rb")

        format = int.from_bytes(file.read(4), "little")
        isBigEndianInt = int.from_bytes(file.read(4), "little")
        bigEndianMode = isBigEndianInt == 1
        isSysB = self.isSysBigEndian()
        self.eigenvectorSize =int.from_bytes(file.read(4), "big" if bigEndianMode else "little")
        self.numberEigenvectors = int.from_bytes(file.read(4), "big" if bigEndianMode else "little")
        #print("Eigenvector size: " + repr(self.eigenvectorSize));

        for x in range(0, self.numberEigenvectors):
            tempJ = self.convertDouble( (struct.unpack('d', file.read(8))[0]), bigEndianMode, isSysB)
            tempBz = self.convertDouble( (struct.unpack('d', file.read(8))[0]), bigEndianMode, isSysB)
            tempBx = self.convertDouble( (struct.unpack('d', file.read(8))[0]), bigEndianMode, isSysB)
            tempEV = self.convertDouble( (struct.unpack('d', file.read(8))[0]), bigEndianMode, isSysB)

            tempEigset = IsingEigenset(tempJ, tempBz, tempBx, tempEV)
            vector = []
            for j in range (0, self.eigenvectorSize):
                vector.append(self.convertDouble(struct.unpack('d', file.read(8))[0], bigEndianMode, isSysB))
            tempEigset.Eigenvector = copy.deepcopy(vector)
            self.eigenpairs.append(copy.deepcopy(tempEigset))

        
