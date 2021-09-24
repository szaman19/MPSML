CC := g++
CFLAGS := -Wno-sign-compare -g -std=c++17 -O2
INCLUDES := -I . -I EigenStuff/eigen -I EigenStuff/eigen/Eigen -I EigenStuff/eigen/unsupported

all: eigen_checker.cpp
	$(CC) eigen_checker.cpp -o eigencheck $(INCLUDES) $(CFLAGS)