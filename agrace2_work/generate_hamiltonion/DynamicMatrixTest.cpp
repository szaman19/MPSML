#include"DynamicMatrix.cpp"

int main(){
    DynamicMatrix a(2, 2);
    a.set(0,0, 1.0);
    a.set(1,1, -1.0);

    
    DynamicMatrix b(2, 2);
    b.set(0,0,-1.0);
    b.set(0,1,1.0);
    b.set(1,0,1.0);
    b.set(1,1,-1.0);

    std::cout << b;
    
    DynamicMatrix c = a.tensor(b);
    std::cout << "set matrices" << std::endl;
    std::cout << c;

}