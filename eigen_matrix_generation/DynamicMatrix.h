#pragma once
#include<iostream>
#include<unordered_map>
#include<iomanip>
#include<sstream>
#include<vector>
#include<Eigen/sparse>
#include<algorithm>

class DynamicMatrix{
    public:
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
        void addInPlaceRef(DynamicMatrix& other );
        void addInPlace(DynamicMatrix other );

        void reserve(long numberValues);
        void setEigenSparseMatrix(Eigen::SparseMatrix<double> *mat);

        DynamicMatrix & operator=(const DynamicMatrix& other);
        DynamicMatrix(const DynamicMatrix & other);

        std::string printLatex();

        long getRows(){
            return rows;
        }
        long getCols(){
            return cols;
        }

        friend std::ostream& operator<<(std::ostream& os, const DynamicMatrix & matrix); 
        void clean();
    private:
        std::unordered_map<long, double> matrixEntries;
        long cols;
        long rows;
        bool hasCheckedEndian;
        bool bigEndian;

};