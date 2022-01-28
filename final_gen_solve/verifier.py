import EigensetReader
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

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

def showFig(mat, xmin, xmax, xnum, ymin, ymax, ynum):
    X = np.linspace(xmin, xmax, xnum)
    Y = np.linspace(ymin, ymax, ynum)
    X, Y = np.meshgrid(X,Y)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    surf = ax.plot_surface(X, Y, mat, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    
    ax.set_zlim(0, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    # ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig("plot.png", format="png")


def verify(file_str, print_surface_plot):
    e = EigensetReader.Eigenset()
    e.read(file_str)
    num_eigenvectors = e.numberEigenvectors
    eigenvectorSize = e.eigenvectorSize
    #Find the number of cubits
    #2^(number of qbits)^2
    
    num_qbits_raw = math.sqrt(math.log2(eigenvectorSize))

    num_qbits = round(num_qbits_raw)
    print("Got " + repr(num_qbits) + " qbits")
    magvec = seq_to_magnetization(seq_gen(num_qbits**2), num_qbits**2)

    #Check what BX and BZ values we have in the eigenset. This will help us output the magnetization values later.
    BxValues = []
    BzValues = []

    for eigenpair in e.eigenpairs:
        if(eigenpair.Bx not in BxValues ):
            BxValues.append(eigenpair.Bx)
        if(eigenpair.Bz not in BzValues):
            BzValues.append(eigenpair.Bz)
    
    BxValues.sort()
    BzValues.sort()

    magnetization_matrix = np.zeros((len(BzValues),len(BxValues)))

    for eigenpair in e.eigenpairs:
        eigenvector_np = np.array(eigenpair.Eigenvector)
        magnetization = np.dot(magvec, eigenvector_np**2)
        index_i = BzValues.index(eigenpair.Bz)
        index_j = BxValues.index(eigenpair.Bx)
        magnetization_matrix[index_i][index_j] = magnetization

    if(print_surface_plot == True):
        #Print Bz values along the top

        for BxVal in BxValues:
            print("," + repr(BxVal), end='')
        print()

        #Print the meat and cheese
        for i in range(0, len(BzValues)):
            print(repr(BzValues[i]), end='')
            for j in range(0, len(BxValues)):
                print("," + repr(magnetization_matrix[i][j]), end='')
            
            print()
    
    showFig(magnetization_matrix,  BzValues[0], BzValues[len(BzValues)-1], len(BzValues), BxValues[0], BxValues[len(BxValues)-1], len(BxValues) )

    return magnetization_matrix


def main():
    filename = "results.eigenset"
    
    x = verify(filename, True)
    #print(x)

    


main()