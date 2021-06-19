#include"DynamicMatrix.h"


DynamicMatrix::DynamicMatrix(long row, long col){
    this->rows = row;
    this->cols = col;
}

DynamicMatrix::DynamicMatrix(std::string input){
    
}

DynamicMatrix::DynamicMatrix(const DynamicMatrix & other){
    this->rows = other.rows;
    this->cols = other.cols;
    matrixEntries.clear();
    for(auto const& x : other.matrixEntries){
        matrixEntries.insert(std::pair<long, double>(x.first, x.second));
    }

}

double DynamicMatrix::get(long i, long j) const{
    long elementPosition = j * rows + i;
     
    std::map<long, double>::const_iterator pos = matrixEntries.find(elementPosition);
    if (pos == matrixEntries.end()) {
    
    }
    else {
        return pos->second;

    }
    return 0.0;
}

void DynamicMatrix::set(long i, long j, double d){
    long elementPosition = j * rows + i;

    if(matrixEntries.find(elementPosition) != matrixEntries.end()){
        matrixEntries.erase(elementPosition);
    }
    if(d != 0.0){
        matrixEntries.insert(std::pair<long, double>(elementPosition, d));
    }

}

bool DynamicMatrix::isNotSparse(long i, long j) const{
    long elementPosition = j * rows + i;

    if(matrixEntries.find(elementPosition) != matrixEntries.end()){
        return true;
    }

   return false;

}

//Scalar multiplicaton
DynamicMatrix DynamicMatrix::operator*(double d){
    DynamicMatrix out(this->rows, this->cols);
    for(long i = 0; i < rows; i++){
        for(long j = 0; j < cols; j++){
            if(isNotSparse(i, j)){
                out.set(i, j, get(i,j) * d);
            }
        }
    }
    return out;

}

DynamicMatrix DynamicMatrix::operator*(const DynamicMatrix & other){
    if(this->cols == other.rows){
        DynamicMatrix out(this->rows, other.cols);
        for(long i = 0; i < this->rows; i++){
            for(long j = 0; j < other.cols; j++){
                double value = 0.0;
                for(long k = 0; k < this->cols; k++){
                    value += this->get( i, k) * other.get( k, j);
                }
                out.set(i,j,value);
            }
        }
        return out;
    }
    return *this;
}

DynamicMatrix& DynamicMatrix::operator=(const DynamicMatrix& other){
    if(this != &other){
        this->rows = other.rows;
        this->cols = other.cols;
        matrixEntries.clear();
        for(auto const& x : other.matrixEntries){
            matrixEntries.insert(std::pair<long, double>(x.first, x.second));
        }
    }
    return *this;
}

DynamicMatrix DynamicMatrix::operator+(const DynamicMatrix& other){
    DynamicMatrix output (this->rows, this->cols);
    if(other.rows == this->rows && other.cols == this->cols){
        for(long i = 0; i < this->rows; i++){
            for(long j = 0; j < this->cols; j++){
                double result = this->get(i,j)+ other.get(i,j);
                output.set(i,j, result);
            }

        }
    }
    return output;

}
DynamicMatrix DynamicMatrix::operator-(const DynamicMatrix& other){
    DynamicMatrix out(this->rows, this->cols);
    if(other.rows == this->rows && other.cols == this->cols){
        for(long i = 0; i < this->rows; i++){
            for(long j = 0; j < this->cols; j++){
                out.set(i,j, this->get(i,j)- other.get(i,j));
            }
        }
    }
    return out;

}

std::string DynamicMatrix::printLatex(){
    std::string output = "\\begin{bmatrix} \n" ;
    for(long i = 0; i < this->rows; i++){
        for(long j = 0; j < this->cols; j++){
            std::stringstream number;
            number.precision(2);
            number << this->get(i,j);
            output += number.str();
            if( j != this->cols - 1) output += " & ";
        }
        if(i != this->rows - 1) output += "\\\\ \n ";
    }
    output += "\\end{bmatrix}";

    return output;
}


DynamicMatrix DynamicMatrix::tensor(DynamicMatrix & other){
    /*
    define new matrix as
    A (x) B = 
    | a11 B ... a1n B |
    |  ...            | 
    |  ...            |
    | am1 B ... amn B |
    */

    long newRow = this->rows * other.rows;
    long newCol = this->cols * other.cols;
    DynamicMatrix output(newRow, newCol);

    for(long i = 0; i < this->rows; i++){
        for(long j = 0; j < this->cols; j++){
            double multiplier = this->get(i,j);
            for(long x = 0; x < other.rows; x++){
                for(long y = 0; y < other.cols; y++){
                    if(other.isNotSparse(x, y)){
                        output.set(other.rows * i + x, other.cols * j + y, other.get(x,y) * multiplier);
                    }
                }
            }
        }
    }

    return output;

}

std::ostream& operator<<(std::ostream& os, const DynamicMatrix& matrix){

    for(long i = 0; i < matrix.rows; i++){
        for(long j = 0; j < matrix.cols; j++){
            os << std::setfill(' ') << std::setw(5) << matrix.get(i,j) << "   ";
        }
        os << std::endl;
    }
    return os;
}

void DynamicMatrix::save(std::string filename){
    //Saves the matrix in a map format
    std::ofstream saveFile;
    saveFile.open(filename);
    saveFile << rows << std::endl;
    saveFile << cols << std::endl;
    for(auto const& x : matrixEntries){
        saveFile << x.first << " " << x.second <<std::endl;
    }
    saveFile.close();


}