#Python reader
#Andy Grace, Josh Szachar 2021

class IsingEigenset:
    J
    Bx
    Bz
    Eigenvalue
    Eigenvector = []

    def __init__(self, J, Bx, Bz, Eigenvalue)
        self.J = J
        self.Bx = Bx 
        self.Bz = Bz
        self.Eigenvalue = Eigenvalue

    def addEigenvectorValue(value)
        Eigenvector.append(value)

class Eigenset:
    version = 1
    isBigEndian = false
    eigenvectorSize
    numberEigenvectors
    eigenpairs = []

    def isSysBigEndian(self):
        if sys.byteorder == "little":
            self.isBigEndian = false
            return false;
        else:
            self.isBigEndian = true
            return true;

    def convertInt(integer, fromEndian, toEndian):
        if fromEndian != toEndian:
            return struct.unpack("<I", struct.pack(">I", integer))[0]
        return integer

    def convertDouble(double, fromEndian, toEndian):
        if fromEndian != toEndian:
            return struct.unpack("<d", struct.pack(">d", double))[0]
        return double

    def read(self, filename):
        bigEndianMode = false
        format = 0
        isBigEndianInt = 0

        file = open(filename, "rb")

        format = int.fromBytes(file.read(4), "little")
        isBigEndianInt = int.fromBytes(file.read(4), "little")
        bigEndianMode = isBigEndianInt == 1
        isSysB = isSysBigEndian()
        self.numberEigenvectors = convertInt( int.fromBytes(file.read(4), "big" if bigEndianMode else "little"))
        self.eigenvectorSize = convertInt( int.fromBytes(file.read(4), "big" if bigEndianMode else "little"))

        for x in range(0, numberEigenvectors)
            tempJ = convertDouble(struct.unpack('d', file.read(8))[0], bigEndianMode, isSysB)
            tempBz = convertDouble(struct.unpack('d', file.read(8))[0], bigEndianMode, isSysB)
            tempBx = convertDouble(struct.unpack('d', file.read(8))[0], bigEndianMode, isSysB)
            tempEV = convertDouble(struct.unpack('d', file.read(8))[0], bigEndianMode, isSysB)

            tempEigset = IsingEigenset(tempJ, tempBz, tempBx, tempEV)

            for j in range (0, eigenvectorSize)
                tempEigset.addEigenvectorValue(convertDouble(struct.unpack('d', file.read(8))[0], bigEndianMode, isSysB))
            
            self.eigenpairs.append(tempEigset)
