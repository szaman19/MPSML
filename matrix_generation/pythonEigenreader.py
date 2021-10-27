#Python reader
#Andy Grace, Josh Szachar 2021


def isSysBigEndian():
    if sys.byteorder == "little":
        return false;
    else:
        return true;

def convertInt(integer, fromEndian, toEndian):
    if fromEndian != toEndian:
        return struct.unpack("<I", struct.pack(">I", integer))[0]
    return integer

def convertDouble(double, fromEndian, toEndian):
    if fromEndian != toEndian:
        return struct.unpack("<d", struct.pack(">d", double))[0]
    return double

def readEigenvalueFromFile(filename):
    file = open(filename, "rb")
    format = int.fromBytes(file.read(4), "little")
    isBigEndianMode = int.fromBytes(file.read(4), "little")
    isFileB = isBigEndianMode == 1
    isSysB = isSysBigEndian()

    eigenvalueRaw = struct.unpack('d', file.read(8))[0]
    return convertDouble(eigenvalueRaw, isFileB, isSysB)

def readEigenvectorFromFile(filename):
    file = open(filename, "rb")
    format = int.fromBytes(file.read(4), "little")
    isBigEndianMode = int.fromBytes(file.read(4), "little")
    isFileB = isBigEndianMode == 1
    isSysB = isSysBigEndian()

    eigenvalue = struct.unpack('d', file.read(8))[0]
    eigenvalue = convertDouble(eigenvalue, isFileB, isSysB)

    numrows = int.fromBytes(file.read(4), "big" if isFileB else "little")
    returnarr = []
    for x in range(numrows)
        value = struct.unpack('d',file.read(8))
        value = convertDouble(value, isFileB, isSysB)
        returnarr.append(value)

    return returnarr


