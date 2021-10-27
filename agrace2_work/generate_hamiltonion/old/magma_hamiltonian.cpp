
#include <iostream>
#include <tgmath.h>
#include "DynamicMatrix.cpp"
/*
Hamiltonian Generator and Solver
Andrew T Grace, Binghamton University.
This code probably has some license attached to it, so contact my superiors before you copy it.
*/

/*
This code uses the MAGMA solver. Column-major layout for general matrices

*/

DynamicMatrix sigma_hat_x_i(2, 2);

DynamicMatrix sigma_hat_z_i(2, 2);
DynamicMatrix I4(2, 2);
DynamicMatrix generateSigma(char spinDirection, int i, int N)
{

    //Finds sigma (spin) (i)
    //Basically I tensored with itself N-1 times, then tensored with corresponding sigma
    //then THAT all tensored with I tensored with itself N times.

    //Generate new dynamic matrices
    DynamicMatrix output(1, 1);
    output.set(0, 0, 1.0);

    for (int x = 1; x < N+1; x++)
    {

        //std::cout << "x = " << x  <<std::endl << output << std::endl <<std::endl;
        if (x == i)
        {
            if (spinDirection == 'x')
            {
                
                output = output.tensor(sigma_hat_x_i);
            }
            else
            {
                
                output = output.tensor(sigma_hat_z_i);
            }
        }
        else
        {
            
            output = output.tensor(I4);
        }
        
    }

    return output;
}

DynamicMatrix generateDualInteraction(int i, int j, int N){
    DynamicMatrix output(1, 1);
    output.set(0, 0, 1.0);
    
    for (int x = 1; x < N+1; x++)
    {
        //std::cout << "x = " << x  <<std::endl << output << std::endl <<std::endl;
        if (x == i)
        {
            output = output.tensor(sigma_hat_z_i);
        }
        else if(x == j){
            output = output.tensor(sigma_hat_z_i);
        }
        else
        {
            output = output.tensor(I4);
        }
    }
    return output;


}



DynamicMatrix generateHamiltonian(int latice_size, double J, double Bx, double Bz)
{
    /*
    Hamiltonian Definition

    H(sigma) = J * { Sum[i] (Sum[j = i+x, i+y] (sigma[z,i]*sigma[z,j]))} - {Bx * Sum[i](sigma[x,i])} - {Bz * Sum[i](sigma[z, i])}  

    Notes: First part of equation: J * (sigmoid of all adjacent sigmoid values)

    spin[x,i] =     | 0 1 |
                    | 1 0 |

    spin[y,i] =     | 0 -i |
                    | i  0 |

    spin[z,i] =     | 1 0 |
                    | 0 1 |

    Calculating sigma:
    sigma[u,i] = I (x) I (x) ... (x) spin(u, i) (x) ... (x) IN

    where I = 2x2 identity matrix
    */

    /*
    Assuming 9 qubits

    0 - 1 - 2
    |   |   |
    3 - 4 - 5
    |   |   |
    6 - 7 - 8

    0-1  q + 1
    0-3  q + sqrt(9)
    1-2  q + 1
    1-4  1 + sqrt(9)
    2-5  2 + sqrt(9)
    3-4  q + 1
    3-6  q + sqrt(9)
    4-5  q + 1
    4-7  q + sqrt(3)
    5-8  q + sqrt(3)
    6-7  q + 1
    7-8  q + 1


    0 - 1
    |   |
    2 - 3

    0-1
    0-2
    1-3
    2-3
    

    */
    int N = latice_size * latice_size;
    DynamicMatrix JTerms(N * N, N * N);

    for (int q = 1; q < N+1; q++)
    {
        //Find uncounted adjacent terms

        if (q + 1 < (N+1) && (q + 2) % latice_size != 0)
        {            
            JTerms = JTerms + generateDualInteraction(q, q+1, N);
        }
        if (q + latice_size < (N+1))
        {
            JTerms = JTerms + generateDualInteraction(q, q + latice_size, N);
        }

    }

    DynamicMatrix BzTerms(1, 1);
    BzTerms.set(0, 0, 1.0);
    for (int q = 1; q < N+1; q++)
    {
        if (q == 1)
        {
            BzTerms = generateSigma('z', q, N);
        }
        else
        {
            BzTerms = BzTerms + generateSigma('z', q, N);
        }
    }

    DynamicMatrix BxTerms(1, 1);
    BxTerms.set(0, 0, 1.0);
    for (int q = 1; q < N+1; q++)
    {
        if (q == 1)
        {
            BxTerms = generateSigma('x', q, N);
        }
        else
        {
            BxTerms = (BxTerms + generateSigma('x', q, N));
        }
    }

    JTerms = JTerms * J;
    BxTerms = BxTerms * Bx;
    BzTerms = BzTerms * Bz;

    DynamicMatrix output = JTerms + (BxTerms * -1) + (BzTerms * -1);
    return output;
}

int main()
{
    sigma_hat_x_i.set(0, 1, 1.0);
    sigma_hat_x_i.set(1, 0, 1.0);

    sigma_hat_z_i.set(0, 0, 1.0);
    sigma_hat_z_i.set(1, 1, -1.0);

    I4.set(0, 0, 1.0);
    I4.set(1, 1, 1.0);

    std::cout << "J Term:" << std::endl
               << generateHamiltonian(2, 1, 0, 0).printLatex() << std::endl;
    std::cout << "Bx Term:" << std::endl
              << generateHamiltonian(2, 0, -1, 0).printLatex() << std::endl;
    std::cout << "Bz Term:" << std::endl
              << generateHamiltonian(2, 0, 0, -1).printLatex() << std::endl;

    std::cout << "Complete, all scalars = 1:" << std::endl << generateHamiltonian(2, 1,1,1).printLatex() <<std::endl;

    std::cout << "sigma x2 " << std::endl << generateSigma('x', 2, 4);

    
}