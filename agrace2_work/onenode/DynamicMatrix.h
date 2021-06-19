#pragma once
#include<iostream>
#include<map>
#include <iomanip>
#include<sstream>
#include<fstream>
//A helper class for storing matrices, and getting them in a representation for MAGMA
//Column-major calculation (j * rows) + i, sparse representation

class DynamicMatrix{
    public:

        void save(std::string filename);
        double * getForMagma();
        DynamicMatrix(long row, long col);
        DynamicMatrix(std::string input);
        DynamicMatrix(){
            cols = 1;
            rows = 1;
            set(0,0,1.0);
        }
        void set(long i, long j, double d);
        double get(long i, long j) const;
        bool isNotSparse(long i, long j) const;
        //Tensor Product
        DynamicMatrix  tensor(DynamicMatrix & other);
        DynamicMatrix operator*(double d);
        DynamicMatrix & operator=(const DynamicMatrix& other);
        DynamicMatrix operator+(const DynamicMatrix& other);
        DynamicMatrix operator*(const DynamicMatrix& other);
        DynamicMatrix operator-(const DynamicMatrix& other);
        DynamicMatrix(const DynamicMatrix & other);

        std::string printLatex();

        long getRows(){
            return rows;
        }
        long getCols(){
            return cols;
        }

        friend std::ostream& operator<<(std::ostream& os, const DynamicMatrix & matrix); 

    private:
        std::map<long, double> matrixEntries;
        long cols;
        long rows;

};