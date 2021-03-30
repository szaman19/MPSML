#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from qutip import *



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


#plt.spy(H(6,1,1,1).full())


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


# In[5]:


def unitary_from_qutip(qutip_eigenstates_object):
    dim = len(qutip_eigenstates_object[0])
    
    U = np.zeros( (dim,dim), np.float64 )
    
    for i in range(dim):
        U[i,:] = np.real(qutip_eigenstates_object[1][i].full()).flatten()
        
    return U

def eigenvalues_from_qutip(qutip_eigenstates_object):

    return np.real(qutip_eigenstates_object[0])


# In[6]:


def flattened_hamiltonian_on_meshgrid(N, J, meshgrid):
    instances = np.size(meshgrid,0)
    dim = 2**N
    
    mat = np.zeros( (instances, dim**2) , np.float64)
    
    for i in range(instances):
        mat[i,:] = np.real(  H(N, meshgrid[i,0], meshgrid[i,1], J).full()  ).flatten()
    
    return mat.flatten()


# In[7]:


def flattened_eigenstates_on_meshgrid(N, J, meshgrid): 
    instances = np.size(meshgrid,0)
    dim = 2**N
    
    vecs = np.zeros((instances, dim**2), np.float64)
    enes = np.zeros((instances, dim**2), np.float64)
    
    for i in range(instances):
        tmp = H(N, meshgrid[i,0], meshgrid[i,1], J).eigenstates()
        vecs[i,:] = unitary_from_qutip(tmp).flatten()
        enes[i,:] = np.diag(eigenvalues_from_qutip(tmp)).flatten()
        
    return vecs.flatten(),enes.flatten()


# In[8]:


def single_eigenstate_from_flattened_meshgrid(meshgrid, unitary, eigenvalues, dim, state):
    instances = np.size(meshgrid, 0)
    
    eigenstate = np.zeros((instances, dim), np.float64)
    eigenenerg = np.zeros((instances, dim), np.float64)
    
    for i in range(instances):
        lb = state * dim + i * dim**2
        ub = state * dim + i * dim**2 + dim
        eigenstate[i,:] = unitary[lb:ub]
        eigenenerg[i,:] = eigenvalues[lb:ub]
        
    return eigenstate.flatten(), eigenenerg.flatten()



def write_bin_data(N, J, Bx_min, Bx_max, Bz_min, Bz_max, nBx, nBz, dir_prefix, test_or_train):

    meshgrid = BxBz_meshgrid(Bx_min, Bx_max, Bz_min, Bz_max, nBx, nBz)
    
    instances = np.size(meshgrid,0)
    dim = 2**N

    field_input_gs_train_header = np.array([instances,      2, dim]   , np.intc)
    hamil_input_gs_train_header = np.array([instances, dim**2, dim]   , np.intc)   
    field_input_un_train_header = np.array([instances,      2, dim**2], np.intc)
    hamil_input_un_train_header = np.array([instances, dim**2, dim**2], np.intc)
    
    field_input_gs_train_path = f"{dir_prefix}{N}-site-field-gs-{test_or_train}.bin"
    hamil_input_gs_train_path = f"{dir_prefix}{N}-site-hamil-gs-{test_or_train}.bin"
    field_input_un_train_path = f"{dir_prefix}{N}-site-field-un-{test_or_train}.bin"
    hamil_input_un_train_path = f"{dir_prefix}{N}-site-hamil-un-{test_or_train}.bin"
    
    gs_prior_header = np.array([instances, dim, 1]  , np.intc)
    un_prior_header = np.array([instances, dim, dim], np.intc)
    
    gs_prior_path = f"{dir_prefix}{N}-site-gs-prior-{test_or_train}.bin"
    un_prior_path = f"{dir_prefix}{N}-site-un-prior-{test_or_train}.bin"

    hamiltonian = flattened_hamiltonian_on_meshgrid(N, J, meshgrid)
    
    unitary,eigenvalues = flattened_eigenstates_on_meshgrid(N, J, meshgrid)
    
    gss,gse = single_eigenstate_from_flattened_meshgrid(meshgrid, unitary, eigenvalues, dim, 0)
    
    with open(field_input_gs_train_path, 'wb') as file:
        field_input_gs_train_header.tofile(file)
        meshgrid.tofile(file)
        gss.tofile(file)
        
    with open(hamil_input_gs_train_path, 'wb') as file:
        hamil_input_gs_train_header.tofile(file)
        hamiltonian.tofile(file)
        gss.tofile(file)
        
    with open(field_input_un_train_path, 'wb') as file:
        field_input_un_train_header.tofile(file)
        meshgrid.tofile(file)
        unitary.tofile(file)
        
    with open(hamil_input_un_train_path, 'wb') as file:
        hamil_input_un_train_header.tofile(file)
        hamiltonian.tofile(file)
        unitary.tofile(file)        
        
    with open(gs_prior_path, 'wb') as file:
        gs_prior_header.tofile(file)
        hamiltonian.tofile(file)
        gse.tofile(file)
        
    with open(un_prior_path, 'wb') as file:
        un_prior_header.tofile(file)
        hamiltonian.tofile(file)
        eigenvalues.tofile(file)





