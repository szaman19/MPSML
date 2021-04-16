#include<iostream>
#include<map>
//A helper class for storing matrices, and getting them in a representation for MAGMA
//Column-major calculation (j * rows) + i, sparse representation

class DynamicMatrix{
    public:
        double * getForMagma();
        DynamicMatrix(int row, int col);
        DynamicMatrix(std::string input);
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


        friend std::ostream& operator<<(std::ostream& os, const DynamicMatrix & matrix); 

    private:
        std::map<int, double> matrixEntries;
        int cols;
        int rows;

};