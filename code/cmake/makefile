###############################################################################
#                                                                             #
#        Generic and Simple Makefile for this C++ file structure              #
#                                                                             #
###############################################################################

#------------ Compilers ------------#
CXX := g++
#CXX := icpc
#CXX := clang++

#------------ File Tree ------------#
SRC_DIR := src
OBJ_DIR := build
INC_DIR := include
TARGET := bin/phynet.x

#------------ Variables ------------#
SRCEXT := cpp
SRC := $(shell find $(SRC_DIR) -type f -name *.$(SRCEXT))
OBJ := $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(SRC:.$(SRCEXT)=.o))

#---------- Compile Flags ----------#
CXXFLAGS :=-std=c++11 #-g -march=native 

#---------- Include Flags ----------#
CXXINC :=-I$(INC_DIR) -I/usr/local/include/Eigen-3.3.7# -I${MKLROOT}/include# -I${MPIROOT}/include

#---------- Sequential Linker Flags -----------#
#LDLIBS := -Wl,--start-group 
#LDLIBS +=${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
#LDLIBS +=${MKLROOT}/lib/intel64/libmkl_sequential.a
#LDLIBS +=${MKLROOT}/lib/intel64/libmkl_core.a #-L$(MPIROOT)/lib -L$(MPIROOT)/lib/release -L$(MPIROOT)/libfabric/lib -lfabric
#LDLIBS +=-Wl,--end-group -liomp5 -lpthread -lm -ldl

#---------- Threaded Linker Flags -----------#
#LDLIBS := -Wl,--start-group 
#LDLIBS +=${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
#LDLIBS +=${MKLROOT}/lib/intel64/libmkl_intel_thread.a
#LDLIBS +=${MKLROOT}/lib/intel64/libmkl_core.a -L$(MPIROOT)/lib -L$(MPIROOT)/lib/release -L$(MPIROOT)/libfabric/lib -lfabric
#LDLIBS +=-Wl,--end-group -liomp5 -lpthread -lm -ldl

ifeq ($(CXX),g++)
	CXXFLAGS +=-O2 #-pedantic -Wall -g -pg #-Werror
	#CXXFLAGS +=-Wshadow -Wstrict-aliasing -Wextra -Wcast-align 
	#CXXFLAGS +=-Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization 
	#CXXFLAGS +=-Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations 
	#CXXFLAGS +=-Wmissing-include-dirs -Wnoexcept -Wold-style-cast 
	#CXXFLAGS +=-Woverloaded-virtual -Wredundant-decls -Wshadow 
	#CXXFLAGS +=-Wsign-conversion -Wsign-promo -Wstrict-null-sentinel 
	#CXXFLAGS +=-Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused
endif

ifeq ($(CXX),icpc)
	CXXFLAGS +=-O2 -ipo -inline-forceinline # -use-intel-optimized-headers -mkl=sequential -mtune=native# -qopenmp
endif

ifeq ($(CXX),clang++)
	CXXFLAGS += -Weverything
	LDLIBS += -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib
	LDLIBS += -L/usr/local/Cellar/boost/1.68.0/lib -lboost_system -lboost_filesystem
endif

###############################################################################
#                                                                             #
#                            End of Definitions                               #
#                                                                             #
###############################################################################

$(TARGET): $(OBJ)
	@echo " Linking...";
	@echo " $(CXX) $(CXXFLAGS) $^ -o $(TARGET) $(LDLIBS)"; $(CXX) $(CXXFLAGS) $^ -o $(TARGET) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.$(SRCEXT)
	@mkdir -p $(OBJ_DIR)
	@echo " $(CXX) $(CXXFLAGS) $(CXXINC) -c -o $@ $<"; $(CXX) $(CXXFLAGS) $(CXXINC) -c -o $@ $<

debug: $(OBJ)
	@echo " Linking...";
	@echo " $(CXX) $(LDLIBS) $^ -o bin/debug "; $(CXX) $(LDLIBS) $^ -o $(TARGET) 

clean:
	@echo " Cleaning..."; 	
	@echo " $(RM) -r $(OBJ_DIR)/* $(TARGET)"; $(RM) -r $(OBJ_DIR)/* $(TARGET)

docs:
	doxygen etc/Doxyfile; cd docs/latex ; pdflatex refman.tex; mv refman.pdf ../userguide.pdf; cd ../../
	#make -C docs/latex

.PHONY: clean docs debug
