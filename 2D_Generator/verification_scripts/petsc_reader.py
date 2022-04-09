import sys
import copy
import struct
import scipy.sparse
import scipy.sparse
from scipy.sparse import dok_matrix




def convert_double(double):
    return struct.unpack("<d", struct.pack(">d", double))[0]

def form_scipy_from_petsc(filename):
    file = open(filename, "rb")
    format = int.from_bytes(file.read(4), "big")
    num_rows = int.from_bytes(file.read(4), "big")
    num_cols = int.from_bytes(file.read(4), "big")
    num_nonzero = int.from_bytes(file.read(4), "big")
    print("got matrix info: " + repr(num_rows) +" rows, " + repr(num_cols) +" cols")
    num_nonzeros_in_each_row = []
    for i in range(0,num_rows):
        num_nonzeros_in_each_row.append(int.from_bytes(file.read(4), "big"))
    
    indicies_of_nonzeros = []

    for i in range(0, num_rows):
        nonzero_indicies_in_this_row = []
        for j in range(0, num_nonzeros_in_each_row[i]):
            nonzero_indicies_in_this_row.append(int.from_bytes(file.read(4), "big"))
        indicies_of_nonzeros.append(copy.deepcopy(nonzero_indicies_in_this_row))

    #Construct the matrix

    output_matrix = dok_matrix((num_rows, num_cols))
    for i in range(0, num_rows):
        for j in indicies_of_nonzeros[i]:
            value = struct.unpack('>d',file.read(8))
            output_matrix[i, j] = value
    
    return output_matrix

        