import sys
import copy
import struct
import scipy.sparse
import scipy.sparse
from scipy.sparse import dok_matrix


def read_next_eigenvector(file):
    curr_size = int.from_bytes(file.read(4), "little")
    curr_eigval = struct.unpack('d',file.read(8))
    curr_j_val = struct.unpack('d',file.read(8))[0]
    curr_bx_val = struct.unpack('d',file.read(8))[0]
    curr_bz_val = struct.unpack('d',file.read(8))[0]
    print(repr(curr_j_val) + " " + repr(curr_bx_val) + " " + repr(curr_bz_val))
    curr_eigenvector = []
    eig_size = (curr_size - 36) / 8
    for i in range(0, int(eig_size)):
        curr_eigenvector.append(struct.unpack('<d',file.read(8))[0])

    return curr_eigval, curr_j_val, curr_bx_val, curr_bz_val, curr_eigenvector

    # outputFile.write((char * ) &size, sizeof(int));
    #         outputFile.write((char * ) &eigenvalue, sizeof(double));
    #         outputFile.write((char * ) &j_val, sizeof(double));
    #         outputFile.write((char * ) &bx_val, sizeof(double));
    #         outputFile.write((char * ) &bz_val, sizeof(double));
