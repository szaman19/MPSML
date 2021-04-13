#include<iostream>
#include<armadillo>

int main(){
    arma::Mat<double> I(2,2,arma::fill::eye);
    arma::Mat<double> sigmaX("0 1; 1 0");
    arma::Mat<double> sigmaZ("0 1; 0 1");

    arma::mat term1 = arma::kron(I, sigmaX);
    arma::mat term2 = arma::kron(I, I);
    arma::mat resultant = arma::kron(term1, term2);
    std::cout << resultant << std::endl;


}