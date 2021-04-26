#include<iostream>
#include<armadillo>

arma::Mat<double> generateSigma(char spinDirection, int i, int N){

    //Finds sigma (spin) (i)
    //Basically I tensored with itself N-1 times, then tensored with corresponding sigma
    //then THAT all tensored with I tensored with itself N times.

    arma::Mat<double> sigmaX("0 1; 1 0");
    arma::Mat<double> sigmaZ("1 0; 0 -1");
    arma::Mat<double> I("1 0; 0 1");
    //Generate new dynamic matrices
    arma::Mat<double> output("1");
    
    for(int x = 0; x < N; x++){

        //std::cout << "x = " << x  <<std::endl << output << std::endl <<std::endl;
        if(x == i-1 ){
            if(spinDirection == 'x'){
                output = arma::kron(output, sigmaX);

            }
            else{
                output = arma::kron(output, sigmaZ);
            }
        }
        else{
            output = arma::kron(output, I);
        }

    }

    return output;

}


int main(){
    int latice_size = 2;
    int N = latice_size * latice_size;
    arma::Mat<double> JTerms(N*N, N*N);
    //std::cout << JTerms << std::endl;
    for(int q = 0; q < N; q++){
        //Find uncounted adjacent terms
        arma::Mat<double> start = generateSigma('z',q,N);
        
        if(q + 1 < N && (q+1) % latice_size != 0){
            arma::Mat<double> leftAdjacent = generateSigma('z',q + 1,N);
            //std::cout << "LA " << q << std::endl << leftAdjacent << std::endl << std::endl;
            JTerms +=  (start * leftAdjacent);

            
        }
        if(q + latice_size < N ){
            arma::Mat<double> underAdjacent = generateSigma('z',q + latice_size,N);
            //std::cout << "UA " << q << std::endl << underAdjacent << std::endl << std::endl;
            JTerms +=  (start * underAdjacent);
            
        }
        std::cout << JTerms << std::endl << std::endl;
    }
    std::cout << JTerms << std::endl;

    /*
    std::cout << generateSigma('x', 0, 3) <<std::endl <<std::endl;
    std::cout << generateSigma('x', 1, 3) <<std::endl <<std::endl;
    std::cout << generateSigma('x', 2, 3) <<std::endl <<std::endl;


    */

}

