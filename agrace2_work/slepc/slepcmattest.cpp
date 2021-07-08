#include"SlepcMat.h"

/* Tester for SlepcMat.h */

static char help[] = "Help message";

int main(int argc, char *argv[]){

    SlepcInitialize(&argc, &argv, (char*)0, help);
    
    {
        /* Scoped to stop memory leaks */
        std::cout << "Slepc Intialized" << std::endl;
        SlepcMat a(2,2);
        std::cout << "Matrix created." <<std::endl;
        a.setValue(1,0,0);
        a.setValue(2,0,1);
        a.setValue(3,1,0);
        a.setValue(4,1,1);

        a.view();
        std::cout << std::endl << std::endl;

        SlepcMat b(2,2);
        b.setValue(1,0,0);
        b.setValue(0,0,1);
        b.setValue(0,1,0);
        b.setValue(1,1,1);

        b.view();
        std::cout << std::endl << std::endl;

        a.tensor_destroy_previous_matrix(b);

        a.view();
        std::cout <<std::endl;

    }

    SlepcFinalize();
    

    
}