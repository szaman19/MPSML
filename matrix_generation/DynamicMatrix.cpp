#include"DynamicMatrix.h"
#include"petscsys.h"
#include<map>
#include<algorithm>

struct Entry{
    long idx;
    double value;
};
inline bool operator<(Entry a, Entry b){
    return a.idx < b.idx;
}

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

//Scalar multiplicaton
DynamicMatrix DynamicMatrix::operator*(double d){
    DynamicMatrix out(this->rows, this->cols);
    for(auto const& x : matrixEntries){
        long row = x.first / this->cols;
        long col = x.first % this->cols;
        out.set(row, col, x.second * d);
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
        matrixEntries =  other.matrixEntries;
    }
    return *this;
}

DynamicMatrix DynamicMatrix::operator+(const DynamicMatrix& other){
    DynamicMatrix output (this->rows, this->cols);
    if(other.rows == this->rows && other.cols == this->cols){
       output.matrixEntries = this->matrixEntries;
       for(auto const& ome : other.matrixEntries){
            double existing = 0.0;
            if(output.matrixEntries.find(ome.first) != output.matrixEntries.end()){
                existing = output.matrixEntries[ome.first];
            }
            output.matrixEntries[ome.first] = existing + ome.second;
       }
    }
    return output;

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
    long newRow = this->rows * other.rows;
    long newCol = this->cols * other.cols;
    DynamicMatrix output(newRow, newCol);
    for(auto const& thismatentry : this->matrixEntries){
        //Decode location of this block
        long startx = thismatentry.first / this->cols;
        long starty = thismatentry.first % this->cols;
        startx *= other.rows;
        starty *= other.cols;
        for(auto const &othermatentry : other.matrixEntries){
            long localx = othermatentry.first / other.cols;
            long localy = othermatentry.first % other.cols;
            output.set(startx + localx, starty + localy, othermatentry.second * thismatentry.second);
        }
    }
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

void DynamicMatrix::savePetsc(std::string filename){
    //Saves the matrix in PETSC format
    //PetscInt    MAT_FILE_CLASSID
    //PetscInt    number of rows
    //PetscInt    number of columns
    //PetscInt    total number of nonzeros
    //PetscInt    *number nonzeros in each row
    //PetscInt    *column indices of all nonzeros (starting index is zero)
    //PetscScalar *values of all nonzeros
    //this->clean();
    std::cout << "mark 1" << std::endl;
    int classID = 1211216;
    int numRows = this->rows;
    int numCols = this->cols;

    //int numNonzero =  matrixEntries.size();
    int numNonZero = 0;
    
    std::vector<int> nonZerosPerRow;
    std::cout << "attempting to allocate " << this->rows << "integers";
    nonZerosPerRow.reserve(this->rows);

    std::cout << "mark 1.1" << std::endl;

    std::vector<int> columnIndicesOfNonzeros;
    std::cout << "attempting to allocate " << this->matrixEntries.size() << "integers" <<std::endl;
    columnIndicesOfNonzeros.reserve(this->matrixEntries.size());

    std::cout << "mark 1.2" << std::endl;

    /*
    std::vector<double> values;
    std::cout << "attempting to allocate " << this->matrixEntries.size() << "doubles" << std::endl;
    values.reserve(this->matrixEntries.size());

    std::cout << "mark 1.3" << std::endl;

    std::map<long, double> orderedMatrixEntries(matrixEntries.begin(), matrixEntries.end());
    */

    

    for(int i = 0; i < this->rows; i++){
        nonZerosPerRow.push_back(0);
    }
    std::cout << "mark 2" << std::endl;

    std::vector<Entry> entries;
    for(auto const& x : matrixEntries){
        if(x.second != 0.0){
            double converted = x.second;
            this->changeToBigEndian((char *) &converted, sizeof(double));
            entries.push_back({x.first, converted});
            numNonZero++;

            nonZerosPerRow[x.first / this->cols]++;
        }
    }

    std::cout << "mark 2.5" << std::endl;

    std::sort(entries.begin(), entries.end());

    std::cout << "mark 3" << std::endl;

    for(long i = 0; i < entries.size(); i++){
        int x = entries[i].idx % cols;
        this->changeToBigEndian((char *) &x, sizeof(int));
        columnIndicesOfNonzeros.push_back(x);
    }
     std::cout << "mark 4" << std::endl;

    std::ofstream saveFile(filename, std::ios::out | std::ios::binary);

    this->changeToBigEndian((char *) &classID, sizeof(int));
    saveFile.write((char*) &classID, sizeof(int));
    
    this->changeToBigEndian((char *) &numRows, sizeof(int));
    saveFile.write((char*) &numRows, sizeof(int));
    
    this->changeToBigEndian((char *) &numCols, sizeof(int));
    saveFile.write((char*) &numCols, sizeof(int));
    
    this->changeToBigEndian((char *) &numNonZero, sizeof(int));
    saveFile.write((char*) &numNonZero, sizeof(int));
    std::cout << "mark 5" << std::endl;
    for(int i : nonZerosPerRow){
        this->changeToBigEndian((char *) &i, sizeof(int));
        saveFile.write((char*) &i, sizeof(int));
    }
    std::cout << "mark 6" << std::endl;

    for(int i : columnIndicesOfNonzeros){
        
        saveFile.write((char*) &i, sizeof(int));
    }
    std::cout << "mark 7" << std::endl;

    for(long i = 0; i < entries.size();i++){
        saveFile.write(reinterpret_cast<char*>(&entries[i].value), sizeof(double));
    }
    std::cout << "mark 8" << std::endl;
    saveFile.close();
}

void DynamicMatrix::checkEndian(){
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    this->bigEndian = (bint.c[0] == 1); 
}

inline void DynamicMatrix::changeToBigEndian(char * c, int size){
    if(!this->hasCheckedEndian) this->checkEndian();
    if(!this->bigEndian){
        for(int i = 0; i < size / 2; i++){
            char temp = c[i];
            c[i] = c[size - 1 - i];
            c[size - i - 1] = temp;
        }

    }
}

void DynamicMatrix::clean(){
    std::unordered_map<long, double>::iterator it = this->matrixEntries.begin();
    while(it != this->matrixEntries.end()){
        if(it->second == 0.0) it = this->matrixEntries.erase(it);
        else it++;
    }
}

