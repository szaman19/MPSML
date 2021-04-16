
#include<iostream>
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

void generateHamiltonian(int latice_size){
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

    JTerms(1,1);
    BzTerms.set(0,0,1.0);
    for(int q = 0; q < N; q++){

        
    }

    DynamicMatrix BzTerms(1,1);
    BzTerms.set(0,0,1.0);
    for(int q = 0; q < N; q++){
        if(q == 0){
            BzTerms = generateSigma('z', q, latice_size);
        }
        else{
            BzTerms += generateSigma('z', q, latice_size);
        }
    }

    DynamicMatrix BxTerms(1,1);
    BxTerms.set(0,0,1.0);
    for(int q = 0; q < N; q++){
        if(q == 0){
            BxTerms = generateSigma('x', q, latice_size);
        }
        else{
            BxTerms += generateSigma('x', q, latice_size);
        }
    }



}


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

int main(){
    sigma_hat_x_i.set(0,1,1.0);
sigma_hat_x_i.set(1,0,1.0);

    sigma_hat_z_i.set(0,0,1.0);
    sigma_hat_z_i.set(1,1,1.0);

    DynamicMatrix x = generateSigma('x', 2, 4);
    std::cout << x;



}