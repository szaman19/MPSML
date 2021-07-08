FLAGS = -Wall -Wextra -std=c++11

all: driver

driver: driver.o IsingHamiltonian.o DynamicMatrix.o
		g++ driver.o IsingHamiltonian.o DynamicMatrix.o -o driver $(FLAGS)

driver.o:	driver.cpp
		g++ -c driver.cpp $(FLAGS)

IsingHamiltonian.o:	IsingHamiltonian.h IsingHamiltonian.cpp
		g++ -c IsingHamiltonian.cpp $(FLAGS)

DynamicMatrix.o: DynamicMatrix.h DynamicMatrix.cpp
		g++ -c DynamicMatrix.cpp $(FLAGS)

clean:
		rm *.o driver