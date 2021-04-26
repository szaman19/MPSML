#pragma once
#include<iostream>
#include<map>
#include <iomanip>
#include<sstream>
//A helper class for storing matrices, and getting them in a representation for MAGMA
//Column-major calculation (j * rows) + i, sparse representation

class DynamicMatrix{
    public:
        double * getForMagma();
        DynamicMatrix(int row, int col);
        DynamicMatrix(std::string input);
        DynamicMatrix(){
            cols = 1;
            rows = 1;
            set(0,0,1.0);
        }
        void set(int i, int j, double d);
        double get(int i, int j) const;
        bool isNotSparse(int i, int j) const;
        //Tensor Product
        DynamicMatrix  tensor(DynamicMatrix & other);
        DynamicMatrix & operator*(double d);
        DynamicMatrix & operator=(const DynamicMatrix& other);
        DynamicMatrix & operator+(const DynamicMatrix& other);
        DynamicMatrix & operator*(const DynamicMatrix& other);
        DynamicMatrix & operator-(const DynamicMatrix& other);
        DynamicMatrix(const DynamicMatrix & other);

        std::string printLatex();

        friend std::ostream& operator<<(std::ostream& os, const DynamicMatrix & matrix); 

    private:
        std::map<int, double> matrixEntries;
        int cols;
        int rows;

};