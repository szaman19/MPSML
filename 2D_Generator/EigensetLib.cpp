#pragma ONCE
#include<iostream>
#include<fstream>
#include<istream>
#include<ostream>
#include<vector>

/*  
Eigenset File Format
    - File Format version - integer
    - Is Big Endian - integer
    - Number of Eigenvalues - int
    - Number of Eigenvector values per eigenvalue - int
    - Array of values: (1 + 1 + 1 + 1 + Number of Eigenvector values)* Number of Eigenvalues doubles
        - Format: J, Bx, Bz, Eigenvalue, <vectors>
*/

class IsingEigenpair{
    public:
    double J;
    double Bx;
    double Bz;
    double Eigenvalue;
    std::vector<double> Eigenvector;

    IsingEigenpair(double J, double Bx, double Bz, double Eigenvalue, std::vector<double> &Eigenvector){
        this->J = J;
        this->Bx = Bx;
        this->Bz = Bz;
        this->Eigenvalue = Eigenvalue;
        this->Eigenvector = Eigenvector;
    }



};

class Eigenset{
    public:
    
    int version = 1;
    int isBigEndian;
    int eigenvectorSize;
    int numberEigenvectors;

    std::vector<IsingEigenpair> eigenpairs;

    void write(std::string fileName){
        std::cout << "eigenvector size: " << eigenvectorSize << "\n";
        int format = 1;
        if(!hasCheckedEndian) detectEndianess();
        int isBigEndianInt = (isBigEndian) ? 1 : 0;

        swapEndianness((char*) &format, sizeof(int), isBigEndian);
        swapEndianness((char*) &isBigEndianInt, sizeof(int), isBigEndian);


        std::ofstream outputFile(fileName, std::ios::out | std::ios::binary );
        outputFile.write((char*) &format, sizeof(int));
        outputFile.write((char*) &isBigEndianInt, sizeof(int));
        outputFile.write((char*) &eigenvectorSize, sizeof(int));
        outputFile.write((char*) &numberEigenvectors, sizeof(int));

        for(int i = 0; i < numberEigenvectors; i++){
            IsingEigenpair pair = eigenpairs[i];
            outputFile.write(reinterpret_cast<char*>(&pair.J), sizeof(double));
            outputFile.write(reinterpret_cast<char*>(&pair.Bx), sizeof(double));
            outputFile.write(reinterpret_cast<char*>(&pair.Bz), sizeof(double));
            outputFile.write(reinterpret_cast<char*>(&pair.Eigenvalue), sizeof(double));
            for(int j = 0; j < pair.Eigenvector.size(); j++)
                outputFile.write(reinterpret_cast<char*>(&pair.Eigenvector[j]), sizeof(double));
 
        }

        outputFile.close();
    }

    void addEigenpair(IsingEigenpair x){
        numberEigenvectors++;
        if(eigenvectorSize == 0) eigenvectorSize = x.Eigenvector.size();
        eigenpairs.push_back(x);
         
    }

    void clear(){
        numberEigenvectors = 0;
        eigenpairs.clear();
    }

    void read(std::string fileName){
        std::ifstream inputFile(fileName, std::ios::out | std::ios::binary);
        bool bigEndianMode = false;
        int format;
        int isBigEndianInt;
        //eigenvalue
        int numRows;

        inputFile.read((char*) &format, sizeof(int));
        inputFile.read((char*) &isBigEndianInt, sizeof(int));
        inputFile.read((char*) &eigenvectorSize, sizeof(int));
        inputFile.read((char*) &numberEigenvectors, sizeof(int));
        swapEndianness((char*) &format, sizeof(int), false);
        swapEndianness((char*) &isBigEndianInt, sizeof(int), false);
        swapEndianness((char*) &numberEigenvectors, sizeof(int), false);
        swapEndianness((char*) &eigenvectorSize, sizeof(int), false);


        // std::cout << format << std::endl;
        // std::cout << isBigEndianInt << std::endl;
        // std::cout << numberEigenvectors << std::endl;
        // std::cout << eigenvectorSize << std::endl;

        if(isBigEndianInt == 1) bigEndianMode = true;

        for(int i = 0; i < numberEigenvectors; i++){
            double tempJ, tempBx, tempBz, tempEV;
            inputFile.read((char*) &tempJ, sizeof(double));
            swapEndianness((char*) &tempJ, sizeof(double), bigEndianMode);

            inputFile.read((char*) &tempBx, sizeof(double));
            swapEndianness((char*) &tempBx, sizeof(double), bigEndianMode);

            inputFile.read((char*) &tempBz, sizeof(double));
            swapEndianness((char*) &tempBz, sizeof(double), bigEndianMode);

            inputFile.read((char*) &tempEV, sizeof(double));
            swapEndianness((char*) &tempEV, sizeof(double), bigEndianMode);

            std::vector<double> tempVecVals;
            for(int j = 0; j < eigenvectorSize; j++){
                double t;
                inputFile.read((char*) &t, sizeof(double));
                swapEndianness((char*) &t, sizeof(double), bigEndianMode);
                tempVecVals.push_back(t);
                //std::cout << t << std::endl;
            }
            eigenpairs.push_back(IsingEigenpair(tempJ, tempBx, tempBz, tempEV, tempVecVals));
            //std::cout << tempJ << " " << tempBx << " " << tempBz << " " << tempEV << " " << std::endl << std::endl;
        }

    }



    private:
        
    bool hasCheckedEndian = false;

    void detectEndianess(){
        union {
            uint32_t i;
            char c[4];
            } bint = {0x01020304};
        this->isBigEndian = (bint.c[0] == 1); 
        this->hasCheckedEndian = true;
    }

    void swapEndianness(char * c, int size, bool isInputBigEndian){
        if(!this->hasCheckedEndian) this->detectEndianess();
        if((isInputBigEndian && !this->isBigEndian) || (!isInputBigEndian && this->isBigEndian)){
            for(int i = 0; i < size / 2; i++){
                char temp = c[i];
                c[i] = c[size - 1 - i];
                c[size - i - 1] = temp;
            }
        }
        
    }

};