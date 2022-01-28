#pragma ONCE
#include<iostream>
#include<vector>
#include "EigensetLib.cpp"


int main(){
    std::string filename = "results.eigenset";

    Eigenset e;

    e.read(filename);
    for(int i = 0; i < e.eigenpairs.size(); i++){

        std::cout << e.eigenpairs[i].J << 
        " " << e.eigenpairs[i].Bx <<
        " " << e.eigenpairs[i].Bz <<
        " " << e.eigenpairs[i].Eigenvalue <<
        std::endl;
    }


}