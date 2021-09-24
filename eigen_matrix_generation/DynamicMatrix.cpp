#include"DynamicMatrix.h"


DynamicMatrix::DynamicMatrix(long row, long col){
    this->rows = row;
    this->cols = col;
    hasCheckedEndian = false;
}

DynamicMatrix::DynamicMatrix(std::string input){
    hasCheckedEndian = false;
}

DynamicMatrix::DynamicMatrix(const DynamicMatrix & other){
    this->rows = other.rows;
    this->cols = other.cols;
    matrixEntries.clear();
    hasCheckedEndian = false;
    matrixEntries =  other.matrixEntries;
}

double DynamicMatrix::get(long i, long j) const{
    long elementPosition = j * rows + i;
     
    std::unordered_map<long, double>::const_iterator pos = matrixEntries.find(elementPosition);
    if (pos != matrixEntries.end()) {
        return pos->second;
    }
    return 0.0;
}

void DynamicMatrix::set(long i, long j, double d){
    long elementPosition = j * rows + i;
    matrixEntries[elementPosition] = d;
}

bool DynamicMatrix::isNotSparse(long i, long j) const{
    long elementPosition = j * rows + i;
    if(matrixEntries.find(elementPosition) != matrixEntries.end()){
        if(matrixEntries.at(elementPosition) != 0.0)
        return true;
    }
   return false;
}

DynamicMatrix& DynamicMatrix::operator=(const DynamicMatrix& other){
    if(this != &other){
        this->rows = other.rows;
        this->cols = other.cols;
        matrixEntries.clear();
        matrixEntries =  other.matrixEntries;
    }
    return *this;
}


void DynamicMatrix::addInPlaceRef( DynamicMatrix& other){
    if(other.rows == this->rows && other.cols == this->cols){
       for(auto const& ome : other.matrixEntries){
            double existing = 0.0;
            if(this->matrixEntries.find(ome.first) != this->matrixEntries.end()){
                existing = this->matrixEntries[ome.first];
            }
            this->matrixEntries[ome.first] = existing + ome.second;
       }
    }

}

void DynamicMatrix::reserve(long numValues){
    this->matrixEntries.reserve(numValues);
}
void DynamicMatrix::addInPlace( DynamicMatrix other){
    if(other.rows == this->rows && other.cols == this->cols){
       for(auto const& ome : other.matrixEntries){
            double existing = 0.0;
            if(this->matrixEntries.find(ome.first) != this->matrixEntries.end()){
                existing = this->matrixEntries[ome.first];
            }
            this->matrixEntries[ome.first] = existing + ome.second;
       }
    }

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




void DynamicMatrix::tensorInPlace( DynamicMatrix & other){
    long oldRow = this->rows;
    long oldCol = this->cols;
    this->rows = this->rows * other.rows;
    this->cols = this->cols * other.cols;
    std::unordered_map<long,double> copyOfMatEntries(this->matrixEntries);
    this->matrixEntries.clear();
    for(auto const& thismatentry : copyOfMatEntries){
        //Decode location of this block
        long startx = thismatentry.first / oldCol;
        long starty = thismatentry.first % oldCol;
        startx *= other.rows;
        starty *= other.cols;
        for(auto const &othermatentry : other.matrixEntries){
            long localx = othermatentry.first / other.cols;
            long localy = othermatentry.first % other.cols;
            this->set(startx + localx, starty + localy, othermatentry.second * thismatentry.second);
        }
    }
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
    //Saves the matrix 
    std::ofstream saveFile;
    saveFile.open(filename);
    saveFile << rows << std::endl;
    saveFile << cols << std::endl;
    for(auto const& x : matrixEntries){
        saveFile << x.first << " " << x.second <<std::endl;
    }
    saveFile.close();
}

void DynamicMatrix::clean(){
    std::unordered_map<long, double>::iterator it = this->matrixEntries.begin();
    while(it != this->matrixEntries.end()){
        if(it->second == 0.0) it = this->matrixEntries.erase(it);
        else it++;
    }
}

void DynamicMatrix::setEigenSparseMatrix(Eigen::SparseMatrix<double> *mat){
    //form triplets
    std::vector<Eigen::Triplet<double>> triplets;
    for(auto const& x : matrixEntries){
        triplets.push_back(Eigen::Triplet<double>(x.first / cols, x.first % cols, x.second));
    }
    mat->setFromTriplets(triplets.begin(), triplets.end());


}