
#include<iostream>
#include<tgmath.h>
#include"DynamicMatrix.cpp"
/*
Hamiltonian Generator and Solver
Andrew T Grace, Binghamton University.
This code probably has some license attached to it, so contact my superiors before you copy it.
*/

/*
This code uses the MAGMA solver. Column-major layout for general matrices

*/

DynamicMatrix sigma_hat_x_i(2,2);

DynamicMatrix sigma_hat_z_i(2,2);
DynamicMatrix generateSigma(char spinDirection, int i, int N){

    //Finds sigma (spin) (i)
    //Basically I tensored with itself N-1 times, then tensored with corresponding sigma
    //then THAT all tensored with I tensored with itself N times.

    //Generate new dynamic matrices
    DynamicMatrix output(1,1);
    output.set(0,0,1.0);
    
    for(int x = 0; x < N; x++){
        if(x == i-1 ){
            if(spinDirection == 'x'){
                output = output.tensor(sigma_hat_x_i);
            }
            else{
                output = output.tensor(sigma_hat_z_i);
            }
        }
        else{
            output = output.tensor(sigma_hat_z_i);
        }

    }

    return output;

}


DynamicMatrix generateHamiltonian(int latice_size, double J, double Bx, double Bz){
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

    DynamicMatrix JTerms(latice_size * latice_size,latice_size * latice_size);
    
    for(int q = 0; q < latice_size; q++){
        //Find uncounted adjacent terms
        DynamicMatrix start = generateSigma('z',q,latice_size);
        if(q + 1 < latice_size){
            DynamicMatrix leftAdjacent = generateSigma('z',q + 1,latice_size);
            JTerms = JTerms + (start * leftAdjacent);
        }
        if(q + log2(latice_size) < latice_size){
            DynamicMatrix underAdjacent = generateSigma('z',q + log2(latice_size),latice_size);
            JTerms = JTerms + (start * underAdjacent);
        }
    }

    DynamicMatrix BzTerms(1,1);
    BzTerms.set(0,0,1.0);
    for(int q = 0; q < latice_size; q++){
        if(q == 0){
            BzTerms = generateSigma('z', q, latice_size);
        }
        else{
            BzTerms = BzTerms + generateSigma('z', q, latice_size);
        }
    }
    
    DynamicMatrix BxTerms(1,1);
    BxTerms.set(0,0,1.0);
    for(int q = 0; q < latice_size; q++){
        if(q == 0){
            BxTerms = generateSigma('x', q, latice_size);
            std::cout << "x1:\n" << BxTerms << std::endl;
        }
        else{
            BxTerms = ( BxTerms + generateSigma('x', q, latice_size));
            std::cout << "xn:\n" << generateSigma('x', q, latice_size)<< std::endl;
        }
        
    }

    
    JTerms = JTerms * J;
    BxTerms = BxTerms * Bx;
    BzTerms = BzTerms * Bz;

    
    DynamicMatrix output = JTerms + (BxTerms * -1) + (BzTerms * -1);

    return output;
    
}



int main(){
    sigma_hat_x_i.set(0,1,1.0);
    sigma_hat_x_i.set(1,0,1.0);

    sigma_hat_z_i.set(0,0,1.0);
    sigma_hat_z_i.set(1,1,1.0);

    DynamicMatrix x = generateHamiltonian(2, 1,1,1);
    std::cout << x;



}