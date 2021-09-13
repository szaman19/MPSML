#include<iostream>
#include<fstream>
#include<stdio.h>
int main(int argc, char* argv[]){
    //CSV Generator 

    //Argument 1: J value
    //Argument 2: Stops in Bx
    //Argument 3: Stops in Bz
    //Argument 4: Initial Bx
    //Argument 5: Stop Bx
    //Argument 6: Iitial Bz
    //Argument 7: Stop Bz;
    //Argument 8: file name

    double J, init_Bx, init_Bz, stop_Bx, stop_Bz;
    int stopsBx, stopsBz;

    std::string out;

    //Populate Arguments

    if(argc < 8 ){
        std::cout << "not enough argments provided." <<std::endl;
        std::cout << "usage: ./gencsv <J value> <Stops in Bx> <Stops in Bz> <Initial Bx> <Stop Bx> <Initial Bz> <Stop Bz> " << std::endl;
        return 0;
    }

    try{
        std::string str(argv[1]);
        J = std::stod(str);
    }
    catch(const std::invalid_argument & e){
        std::cout << "J is wrong." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "J is wrong." <<std::endl;
        return 0;
    }

    try{
        std::string str(argv[4]);
        init_Bx = std::stod(str);
    }
    catch(const std::invalid_argument & e){
        std::cout << "init_Bx is wrong." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "init_Bx is wrong." <<std::endl;
        return 0;
    }

    try{
        std::string str(argv[5]);
        init_Bz = std::stod(str);
    }
    catch(const std::invalid_argument & e){
        std::cout << "init_Bz is wrong." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "init_Bz is wrong." <<std::endl;
        return 0;
    }

    try{
        std::string str(argv[6]);
        stop_Bx = std::stod(str);
    }
    catch(const std::invalid_argument & e){
        std::cout << "stop_Bx is wrong." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "stop_Bx is wrong." <<std::endl;
        return 0;
    }

    try{
        std::string str(argv[7]);
        stop_Bz = std::stod(str);
    }
    catch(const std::invalid_argument & e){
        std::cout << "stop_Bz is wrong." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "stop_Bz is wrong." <<std::endl;
        return 0;
    }

    stopsBx = std::stoi(argv[2]);
    stopsBz = std::stoi(argv[3]);

    std::string fileName(argv[8]);

    std::cout << "Input parameters:" << std::endl;
    std::cout << "J: " << J << std::endl; 
    std::cout << "stopsBx: " << stopsBx << std::endl; 
    std::cout << "stopsBz: " << stopsBz << std::endl;
    std::cout << "init_Bx: " << init_Bx << std::endl; 
    std::cout << "init_Bz: " << init_Bz << std::endl; 
    std::cout << "stop_Bx: " << stop_Bx << std::endl; 
    std::cout << "stop_Bz: " << stop_Bz << std::endl; 

    std::ofstream output(fileName);
    if(output.is_open()){
        for(double i = init_Bx; i < stop_Bx; i += ((stop_Bx - init_Bx)/ stopsBx)){
            for(double j = init_Bz; j < stop_Bz; j += ((stop_Bz - init_Bz)/ stopsBz)){
                output << J << "," << i << "," << j << std::endl;
            }
        }

        output.close();
    }
    



}