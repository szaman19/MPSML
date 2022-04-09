import numpy as np
import math
from time import perf_counter
from scipy.linalg import eigh
import scipy.sparse as sps
import scipy.sparse.linalg
from tensor import tensor
import new_eigenset_reader

pauli_z = tensor(2, {(0, 0): 1, (1, 1): -1})
pauli_x = tensor(2, {(0, 1): 1, (1, 0): 1})
id = tensor(2, {(0, 0): 1, (1, 1): 1})


id_np = np.array([[1,0],[0,1]])
pauli_z_np = np.array([[1,0],[0,-1]])
pauli_x_np = np.array([[0,1],[1,0]])

def generateSelfInteraction(num_qubits, _site, field=pauli_z):
    """Generate the matrix for a single site in x or z transverse field
       num_qubits (int): Number of qubits in the system
       _sits (int): Site for which the matrix is being generated
       field (tensor): Expects either a pauli_x or pauli_z tensor 
    """
    _tensor = id
    if _site == 0:
        _tensor = field
    for i in range(1, num_qubits):
        if (i == _site):
            _tensor = _tensor.tensor_prod(field)
        else:
            _tensor = _tensor.tensor_prod(id)
    return _tensor


def generateNeighborInteractions(interaction_vector):
    """
      Generate the interaction between two sites
      interaction_vector (list): A list of length (num_qubits) with 1 on index i if the 
                                 site is a constituent of the interaction. For example, 
                                 for interaction between site 1 and 2 in a 4 qubit system, 
                                 the vector would be,[1, 1, 0, 0] and for stie 1 and 3, [1, 0, 1, 0]
    """
    pauli_z = tensor(2, {(0, 0): 1, (1, 1): -1})
    id = tensor(2, {(0, 0): 1, (1, 1): 1})

    _tensor = pauli_z if interaction_vector[0] == 1 else id
    for site_i in interaction_vector[1:]:
        field = pauli_z if site_i == 1 else id
        _tensor = _tensor.tensor_prod(field)
    return _tensor


def gen_bz(num_qubits):
    bz = generateSelfInteraction(num_qubits, 0)
    for i in range(1, num_qubits):
        _x = generateSelfInteraction(num_qubits, i)
        bz = tensor.join_tensors(bz, _x)
    return bz


def gen_bx(num_qubits):
    bx = generateSelfInteraction(num_qubits, 0, pauli_x)
    for i in range(1, num_qubits):
        _x = generateSelfInteraction(num_qubits, i, pauli_x)
        bx = tensor.join_tensors(bx, _x)
    return bx

def gen_j(num_qubits):
    unique_interactions = {}
    for i in range(num_qubits):
        lattice_size = round(math.sqrt(num_qubits))
        site_1 = i
        row = site_1 // lattice_size
        col = site_1 % lattice_size
        x_hop = lattice_size*row + (col + 1) % lattice_size
        y_hop = lattice_size*((row + 1) % lattice_size) + col
  
      
        for j in [x_hop, y_hop]:
            interaction_vector = [0]*num_qubits
            interaction_vector[site_1] = 1
            interaction_vector[j] = 1
            def _hash(x): return sum(
                [(2 * interaction_vector[i])**i for i in range(num_qubits)])
            unique_interactions[_hash(interaction_vector)] = interaction_vector
    unique_interactions = list(unique_interactions.values())

    j = generateNeighborInteractions(unique_interactions[0])
    for i in unique_interactions[1:]:
        _x = generateNeighborInteractions(i)
        j = tensor.join_tensors(j, _x)


    return j


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
        return ['00', '01', '10', '11']
    else:
        temp = []

        smaller_vals = seq_gen(num_q-1)
        for i in ['0', '1']:
            for each in smaller_vals:
                temp.append(i+each)
        return temp


def gen_mag_mat(num_qbit, initBx, initBz, stopBx, stopBz, numBx, numBz):
    #j = gen_j(num_qbit)
    j = gen_j_numpy_only(num_qbit)
    bz = gen_bz(num_qbit)
    bx = gen_bx(num_qbit)
    magnetization_vector = seq_to_magnetization(seq_gen(num_qbit), num_qbit)
    magnetizations = np.zeros((numBz, numBx))
    row = 0
    for z_val in np.linspace(initBz, stopBz, numBz):
        col = 0
        for x_val in np.linspace(initBx, stopBx, numBx):
            #test_mat = (j.numpy())- z_val * bz.numpy() - x_val * bx.numpy()
            test_mat = -1.0 * j - z_val * bz.numpy() - x_val * bx.numpy()
            # pretty_print(test_mat, True)
            eigenvalue, eigenvector = eigh(test_mat, eigvals=(0, 0))
            #assert np.allclose(test_mat @ eigenvector, eigenvalue*eigenvector)
            m = np.dot(magnetization_vector, eigenvector**2)
            magnetizations[row][col] = m
            col += 1
        row += 1
    return magnetizations


def check_eigenvalues(num_qbit, initBx, initBz, stopBx, stopBz, numBx, numBz, file_str):
    #j = gen_j(num_qbit)
    print("Generating J", flush=True)
    j = gen_j(num_qbit)
    print("Done. Generating Bz", flush=True)
    bz = gen_bz(num_qbit)
    print("Done. Generating Bx", flush=True)
    bx = gen_bx(num_qbit)

    
    petsc_eigenvalues = []
    gen_eigenvalues = []

    file = open(file_str, "rb")

    while True:
        eof_test = file.read(1)
        if not eof_test:
            break
        file.seek(-1,1)
        petsc_eigval, j_s, bx_s, bz_s, petsc_eigvec = new_eigenset_reader.read_next_eigenvector(file)
        petsc_eigenvalues.append(petsc_eigval)
        test_mat = (-1.0 * j.get_as_scipy_sparse()) - (bz_s * bz.get_as_scipy_sparse()) - (bx_s * bx.get_as_scipy_sparse())
        eigenvalue, eigenvector = sps.linalg.eigsh(test_mat)
        gen_eigenvalues.append(eigenvalue[0])


    error = 0
    for i in range(0, len(petsc_eigenvalues)):
        error += ((gen_eigenvalues[i] - petsc_eigenvalues[i])**2)
   
    np.savetxt("gen_eig.txt", gen_eigenvalues, fmt='%.3f')
    np.savetxt("petsc_eig.txt", petsc_eigenvalues, fmt='%.3f')
    return error / len(petsc_eigenvalues)




def main():


    initBx = 0
    stopBx = 2
    numBx = 20

    initBz = 0
    stopBz = 1
    numBz = 20

    num_qbit = 16

    filename = "results.eigenset"

    #magmat = gen_mag_mat(num_qbit, initBx, initBz, stopBx, stopBz, numBx, numBz)
    
    
    #x = verify(filename, magmat)

    #print(x)

    err = check_eigenvalues(num_qbit, initBx, initBz, stopBx, stopBz, numBx, numBz, filename)
    print(err)

#error log: 2/3/2022, mean squared eigenvalue error of 1260.74?
#error log: 2/10/22, .27289725016561583 mean eigenvalue error


main()