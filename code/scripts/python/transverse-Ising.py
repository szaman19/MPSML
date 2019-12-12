import numpy as np
from qutip import *
import sys
import os


def transverse(N, site):
    op = []
    
    for i in range(N):
        op.append(qeye(2))
        
    op[site] = sigmax()
    
    return tensor(op)

def longitudinal(N, site):
    op = []
    
    for i in range(N):
        op.append(qeye(2))
        
    op[site] = sigmaz()
    
    return tensor(op)

def coupling(N, site):
    if (site == N-1):
        return longitudinal(N,site) * longitudinal(N, 0)
    else:
        return longitudinal(N,site) * longitudinal(N, site+1)

def H(N, Bx, Bz, J):
    term1 = coupling(N, 0) 
    term2 = transverse(N, 0)
    term3 = longitudinal(N, 0)
       
    for site in range(1, N):
        term1 += coupling(N,site)
        term2 += transverse(N,site)
        term3 += longitudinal(N,site)
    
    return (J/N) * term1 - (Bx/N) * term2 - (Bz/N) * term3

def mag(N):
    tmp = longitudinal(N, 0)
    for site in range(1, N):
        tmp += longitudinal(N,site)
    return np.diag(tmp.full().real) / N

def BxBz_meshgrid(Bx_min, Bx_max, Bz_min, Bz_max, nBx, nBz):
    
    mat = np.zeros((nBx*nBz, 2), np.float64)

    x = np.linspace(Bx_min, Bx_max, nBx)
    z = np.linspace(Bz_min, Bz_max, nBz)

    c = 0
    for Bx in x:
        for Bz in z:
            mat[c, :] = Bx,Bz
            c += 1
            
    np.random.shuffle(mat)
    
    return mat


# In[ ]:


def groundstate_wavefunction_from_qutip(qutip_groundstate_object):
    return qutip_groundstate_object[1].full().real.flatten()

def groundstate_eigenenergy_from_qutip(qutip_groundstate_object):
    return qutip_groundstate_object[0].real


# In[ ]:


def values_array_from_qutip_hamiltonian(qutip_hamiltonian):
    tmp = qutip_hamiltonian.full().real
    dim = qutip_hamiltonian.shape[0]
    iu1 = np.triu_indices(dim)
    return tmp[iu1][np.nonzero(tmp[iu1])]
    
def locations_array_from_qutip_hamiltonian(qutip_hamiltonian):
    tmp = qutip_hamiltonian.full().real
    dim = qutip_hamiltonian.shape[0]
    iu1 = np.triu_indices(dim)
    return np.nonzero(tmp[iu1])[0]


# In[ ]:


def flattened_hams_gsvecs_and_enes_on_meshgrid(N, J, meshgrid):
    instances = np.size(meshgrid,0)
    dim = 2**N
    nnz = len(values_array_from_qutip_hamiltonian(H(N, 1, 0.01, -1)))
    
    vecs = np.empty((instances, dim), np.float64)
    enes = np.empty((instances, 1), np.float64)
    matval = np.empty((instances, nnz), np.float64)
    matloc = np.empty((instances, nnz), np.int32)
    
    for i in range(instances):  
        ham = H(N, meshgrid[i,0], meshgrid[i,1], J)
        gs = ham.groundstate()
        vecs[i,:] = groundstate_wavefunction_from_qutip(gs)
        enes[i,:] = groundstate_eigenenergy_from_qutip(gs)
        matval[i,:] = values_array_from_qutip_hamiltonian(ham)
        matloc[i,:] = locations_array_from_qutip_hamiltonian(ham)
        
    return matval.flatten(),matloc.flatten(),vecs.flatten(),enes.flatten()


def write_data(N, J, meshgrid, dir_prefix):
    # instances = np.size(meshgrid,0)
    # dim = 2**N
    # nnz = len(values_array_from_qutip_hamiltonian(H(N, 1, 0.01, -1)))
    
    # header = np.array([instances, nnz, dim], np.int32)
    # fields = meshgrid[:,0].flatten()
    
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    
    szdiag_path = f"{dir_prefix}szdiag.bin"
    # header_path = f"{dir_prefix}header.bin"
    # fields_path = f"{dir_prefix}fields.bin"
    # matval_path = f"{dir_prefix}matval.bin"
    # matloc_path = f"{dir_prefix}matloc.bin"
    # wavefx_path = f"{dir_prefix}wavefx.bin"
    # energy_path = f"{dir_prefix}energy.bin"

    # matval,matloc,wavefx,energy = flattened_hams_gsvecs_and_enes_on_meshgrid(N, J, meshgrid)

    with open(szdiag_path, 'wb') as file:
        mag(N).tofile(file)
    
    # with open(header_path, 'wb') as file:
        # header.tofile(file)
    
    # with open(fields_path, 'wb') as file:
        # fields.tofile(file)

    # with open(matval_path, 'wb') as file:
        # matval.tofile(file)
    
    # with open(matloc_path, 'wb') as file:
        # matloc.tofile(file)
        
    # with open(wavefx_path, 'wb') as file:
        # wavefx.tofile(file)
        
    # with open(energy_path, 'wb') as file:
        # energy.tofile(file)

J = -1

root = "/home/csingh5/Documents/Guided-Machine-Learning/code/chris/data/input/Ising/"
qubits = [sys.argv[1]]
#dsizes = [40,60,80,100]
dsizes = [int(sys.argv[2])]
sampts = ["single-phase/", "dual-phase/"]
splits = ["training/", "validation/", "testing/"]

Bz_min = 0.01; Bz_max = 0.01; nBz = 1; 
Bx_min = 0.01
    
for qubit in qubits:
    for dsize in dsizes:
        for sampt in sampts:
            for split in splits:
                dir_prefix = root + str(qubit) + "-qubits/" + str(dsize) + "k-instances/" + sampt + split
                
                # figure out max transverse field 
                if sampt == "dual-phase/" or split == "testing/":
                    Bx_max = 2
                else:
                    Bx_max = 0.5
                
                # figure out number of instances 
                if split == "training/":
                    nBx = dsize * 700
                elif split == "validation/":
                    nBx = dsize * 200
                else:
                    nBx = dsize * 100
                    
                grid = BxBz_meshgrid(Bx_min, Bx_max, Bz_min, Bz_max, nBx, nBz)
                write_data(int(qubit), J, grid, dir_prefix)
                print("finished:", dir_prefix, nBx)

