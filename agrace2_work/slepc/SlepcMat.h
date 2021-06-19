/* Objective-oriented SLEPc and PETSc Wrapper for C++ */
/* Andrew Grace, 2021                                 */

#include<iostream>
#include<stdio.h>
#include"slepcsys.h"
#include<slepceps.h>
#include"petscsys.h"


class SlepcMat{
    private:
        

        /*  
        PETSC matrices have several stages to getting an actual usable matrix.
        This wrapper tries to hide that abstraction from the user as much as possible.
        See SlepcMat.cpp for more comments.
        My wrapper keeps the matrix assembled as much as possible
        */

        int rows;
        int cols;
        Mat mat;
        bool assembledState;
        
        /*
        My implementation has routines for creating the matrix and destroying the matrix.
        */
        void createMatrix();
        void destroy();

    public:
        /* Generates a 1x1 matrix */
        SlepcMat();

        /* Generates a rows x cols matrix */
        SlepcMat(int rows, int cols);

        /* Copy constructor */
        SlepcMat(const SlepcMat &other);

        /* Destructor which frees up space */
        ~SlepcMat();

        /* Assignment */
        SlepcMat& operator= (const SlepcMat& other);

        /* Addition */
        SlepcMat operator+ (SlepcMat other);

        /* Scalar multiplication */
        SlepcMat operator* (double x);

        /* Sets a value by creating a new matrix and adding it */
        void setValue(double value, int row, int col);

        /* Tensor product- useful for making hamiltonians.
           A future upgrade may change the behavior of this so that SLEPC does most of the work */
        void tensor_destroy_previous_matrix( SlepcMat & other);

        /* View using PETSC */
        void view();

        /* Use SLEPC to get the lowest eigenvalue and vector */
        void getLowestEigenpair(PetscScalar * vaReal, PetscScalar * vaImaginary, Vec * veReal, Vec * veImaginary);

        void shareMatrix(double *, Mat*, int, int);


};

