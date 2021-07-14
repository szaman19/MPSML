#pragma once
#include<iostream>
#include<unordered_map>
#include <iomanip>
#include<sstream>
#include<fstream>
#include<vector>
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
            hasCheckedEndian = false;
        }
        void set(long i, long j, double d);
        double get(long i, long j) const;
        bool isNotSparse(long i, long j) const;
        void tensorInPlace( DynamicMatrix& other );
        void addInPlace(DynamicMatrix& other );
        void addInPlace(DynamicMatrix other );
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
        void savePetsc(std::string filename);
        void checkEndian();
        void changeToBigEndian(char * c, int size);
        void clean();
    private:
        std::unordered_map<long, double> matrixEntries;
        long cols;
        long rows;
        bool hasCheckedEndian;
        bool bigEndian;

};