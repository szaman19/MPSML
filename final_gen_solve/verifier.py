import EigensetReader
import numpy as np
import math

def seq_to_magnetization(arr_seq, num_qubits):
    mag_vec = []
    for elem in arr_seq:
        magnetization = 0
        for char in elem:
            temp = (int(char)*-2)+1
            magnetization += temp 
        mag_vec.append(magnetization)
    mag_vec = np.array(mag_vec)
    mag_vec = mag_vec / (num_qubits)
    return mag_vec   

def seq_gen(num_q):
    if num_q == 2:
        return ['00','01', '10','11']
    else:
        temp = []
        
        smaller_vals = seq_gen(num_q-1)
        for i in ['0','1']:
            for each in smaller_vals:
                temp.append(i+each)
        return temp 


def verify(file_str):
    e = EigensetReader.Eigenset()
    e.read(file_str)
    num_qbit = e.eigenvectorSize
    num_qbit_2 = math.sqrt(num_qbit)
    
    num_evec = e.numberEigenvectors
    print(num_evec)
    print(num_qbit)
    wavefuncs = np.zeros(( num_evec, num_qbit))
    for i in range(0, num_evec):
        for j in range (0, num_qbit ):
            wavefuncs[i,j] = e.eigenpairs[i].Eigenvector[j] 
    
    magvec = seq_to_magnetization(seq_gen(num_qbit_2), num_qbit_2)
    print(magvec)
    arr = np.matmul(wavefuncs**2, magvec)

    for i in range (0, num_evec):
        print( repr(e.eigenpairs[i].Bx) + ", " + repr(e.eigenpairs[i].Bz) + ", " + repr(arr[i]))

def verify_surface_plot(file_str):
    e = EigensetReader.Eigenset()
    e.read(file_str)
    num_qbit = e.eigenvectorSize
    num_qbit_2 = math.sqrt(num_qbit)
    
    num_evec = e.numberEigenvectors
    Bxs = []
    Bzs = []
    for i in range(0, num_evec):
        if(e.eigenpairs[i].Bx not in Bxs): Bxs.append(e.eigenpairs[i].Bx)
        if(e.eigenpairs[i].Bz not in Bzs): Bzs.append(e.eigenpairs[i].Bz)

    Bxs.sort()
    Bzs.sort()

    wavefuncs = np.zeros(( num_evec, num_qbit))
    for i in range(0, num_evec):
        for j in range (0, num_qbit ):
            wavefuncs[i,j] = e.eigenpairs[i].Eigenvector[j] 

    magvec = seq_to_magnetization(seq_gen(num_qbit_2), num_qbit_2)
    arr = np.matmul(wavefuncs**2, magvec)

    returnarr = [[0 for i in range(0,len(Bxs))] for j in range(0, len(Bzs))]

    for i in range (0, num_evec):
        index1 = Bzs.index(e.eigenpairs[i].Bz)
        index2 = Bxs.index(e.eigenpairs[i].Bx)
        returnarr[index1][index2] = arr[i]

    #print the top indices
    print(",", end = '')
    for x in range(0, len(Bxs)):
        print(repr(Bxs[x]) + ",", end = '')
    print()

    for x in range(0, len(Bzs)):
        print(repr(Bzs[x]) + ",", end='')
        for y in range(0, len(Bxs)):
            print(repr(returnarr[x][y]) + ",", end='')
        print()


def main():
    filename = "results.eigenset"
    
    x = verify_surface_plot(filename)
    print(x)

    


main()