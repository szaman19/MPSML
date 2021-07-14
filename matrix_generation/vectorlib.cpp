#pragma ONCE
#include<iostream>
#include<fstream>
#include<istream>
#include<ostream>
#include<vector>

/*  
File Format:
    - File Format version - integer
    - Is Big Endian - integer
    - Eigenvalue - double
    - Number of values in the vector - integer
    - Array of the values - doubles
*/

class PETSCVectorLoader{
    public:

        PETSCVectorLoader(){

        }

        void readPETSC(std::string fileName){
            //PETSC File format
            //PetscInt      PETSC_VEC_CLASSID
            //PetscInt      number of rows
            //PetscScalar*  array  

            std::ifstream inputFile(fileName, std::ios::out | std::ios::binary);

            int classid;
            int rows;
            inputFile.read((char*) &classid, sizeof(int));
            inputFile.read((char*) &rows, sizeof(int));
            swapEndianness((char*) &classid, sizeof(int), true);
            swapEndianness((char*) &rows, sizeof(int), true);
            for(int i = 0; i < rows; i++){
                double temp; 
                inputFile.read((char*) &temp, sizeof(double));
                swapEndianness((char*) &temp, sizeof(double), true);
                values.push_back(temp);
            }
            inputFile.close();

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
            swapEndianness((char*) &format, sizeof(int), false);
            swapEndianness((char*) &isBigEndianInt, sizeof(int), false);
            if(isBigEndianInt == 1) bigEndianMode = true;

            inputFile.read((char*) &eigenvalue, sizeof(double));
            swapEndianness((char*) &eigenvalue, sizeof(double), bigEndianMode);

            inputFile.read((char*) &numRows, sizeof(int));
            swapEndianness((char*) &numRows, sizeof(int), bigEndianMode);
            
            for(int i = 0; i < numRows; i++){
                double temp;
                inputFile.read((char*) &temp, sizeof(double));
                swapEndianness((char*) &temp, sizeof(double), bigEndianMode);
                values.push_back(temp);
            }
        }

        void write(std::string fileName){
            int format = 1;
            if(!hasCheckedEndian) detectEndianess();
            int isBigEndianInt = (isBigEndian) ? 1 : 0;

            swapEndianness((char*) &format, sizeof(int), isBigEndian);
            swapEndianness((char*) &isBigEndianInt, sizeof(int), isBigEndian);

            int numRows = values.size();

            std::ofstream outputFile(fileName, std::ios::out | std::ios::binary );
            outputFile.write((char*) &format, sizeof(int));
            outputFile.write((char*) &isBigEndianInt, sizeof(int));
            

            outputFile.write(reinterpret_cast<char*>(&eigenvalue), sizeof(double));
            outputFile.write((char*) &numRows, sizeof(int));
            for(int i = 0; i < numRows; i++){
                double temp = values[i];
                outputFile.write(reinterpret_cast<char*>(&temp), sizeof(double));
 
            }

            outputFile.close();
        }

        void setEigenval(double e){
            eigenvalue = e;
        }



    private:
        std::vector<double> values;
        double eigenvalue;
        bool isBigEndian;
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