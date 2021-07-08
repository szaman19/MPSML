#include"SlepcMat.h"
#include"mpi.h"

class SlepcIsingHamiltonian{
    private:
        SlepcMat JMatrix;
        SlepcMat BxMatrix;
        SlepcMat BzMatrix;

        SlepcMat Identity;
        SlepcMat SigmaX;
        SlepcMat SigmaZ;

        int lattice_size;
        int N;

        bool verbose;

        SlepcMat generateSelfInteraction(char spin, int i){
            SlepcMat output(1, 1);
            output.setValue(1.0, 0, 0);
            for (int x = 1; x < N+1; x++)
            {
                if (x == i)
                {
                    if (spin == 'x')
                        output.tensor_destroy_previous_matrix(SigmaX);
                    else
                        output.tensor_destroy_previous_matrix(SigmaZ);
                }
                else
                    output.tensor_destroy_previous_matrix(Identity);
                }
            return output;
        }

        SlepcMat generateAdjacentInteractionZ(int i, int j){
            SlepcMat output(1, 1);
            output.setValue(1.0, 0, 0);
            for (int x = 1; x < N+1; x++)
            {
                if (x == i || x == j)
                    output.tensor_destroy_previous_matrix(SigmaZ);
                else
                    output.tensor_destroy_previous_matrix(Identity);
            }
            return output;
        }

        void generateJMatrix(){
            int matrixDim = 1;
            for(int i = 0; i < N; i++){
                matrixDim = 2 * matrixDim;
            }

            JMatrix= SlepcMat(matrixDim, matrixDim);
            for (int q = 1; q < N+1; q++)
            {
                if (q + 1 < (N+1) && (q + 2) % lattice_size != 0){   
 
                    JMatrix = JMatrix + generateAdjacentInteractionZ(q, q+1);

                }
                if (q + lattice_size < (N+1)){

                    JMatrix = JMatrix + generateAdjacentInteractionZ(q, q + lattice_size);

                }
            }

        }

        void generateBxMatrix(){
            BxMatrix = SlepcMat(1, 1);
            BxMatrix.setValue(1.0, 0, 0);
            for (int q = 1; q < N+1; q++)
            {
                if (q == 1)
                    BxMatrix = generateSelfInteraction('x', q);
                else
                    BxMatrix = BxMatrix + generateSelfInteraction('x', q);
            }

        }

        void generateBzMatrix(){
            BzMatrix = SlepcMat(1, 1);
            BzMatrix.setValue(1.0, 0, 0);
            for (int q = 1; q < N+1; q++)
            {
                if (q == 1)
                    BzMatrix = generateSelfInteraction('z', q);
                else
                    BzMatrix = BzMatrix + generateSelfInteraction('z', q);
            }

        }


    public:

        SlepcIsingHamiltonian(int lattice_size, bool verbose){
            this->lattice_size = lattice_size;
            this->verbose = verbose;

            /* Initialize Identity, SigmaX, SigmaZ */

            if(verbose) PetscPrintf(MPI_COMM_WORLD, "Initializing Pauli Matrices...");
            Identity = SlepcMat(2,2);
            Identity.setValue(1,0,0);
            Identity.setValue(0,0,1);
            Identity.setValue(0,1,0);
            Identity.setValue(1,1,1);

            SigmaX = SlepcMat(2,2);
            SigmaX.setValue(0,0,0);
            SigmaX.setValue(1,0,1);
            SigmaX.setValue(1,1,0);
            SigmaX.setValue(0,1,1);

            SigmaZ = SlepcMat(2,2);
            SigmaZ.setValue(1,0,0);
            SigmaZ.setValue(0,0,1);
            SigmaZ.setValue(0,1,0);
            SigmaZ.setValue(-1,1,1);

            this->N = lattice_size * lattice_size;
            if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0m");
            /* Generate the required matrices */
            if(verbose) PetscPrintf(MPI_COMM_WORLD, "Generating J Matrix...\n");
            generateJMatrix();
            if(verbose) PetscPrintf(MPI_COMM_WORLD,"\033[0;32mDone!\n\033[0mGenerating Bx Matrix\n");
            generateBxMatrix();
            if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0mGenerating Bz Matrix\n");
            generateBzMatrix();
            if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0m");

        }

        void getSolution(double J, double Bx, double Bz){
            if(verbose) PetscPrintf(MPI_COMM_WORLD,"Running scalar matrix multiplications...\n");
            SlepcMat Jmultiplied;
            SlepcMat Bxmultiplied; 
            SlepcMat Bzmultiplied;
            Jmultiplied = JMatrix * J;
            Bxmultiplied = BxMatrix * (Bx * -1.0);
            Bzmultiplied = BzMatrix * (Bz * -1.0);
            if(verbose)PetscPrintf(MPI_COMM_WORLD,"\033[0;32mDone!\n\033[0m\n");


            if(verbose){
            PetscPrintf(MPI_COMM_WORLD, "\nCalculated Matrices:\n\nJMul:\n");
            Jmultiplied.view();

            PetscPrintf(MPI_COMM_WORLD, "\nBxMul:\n");
            Bxmultiplied.view();

            PetscPrintf(MPI_COMM_WORLD,"\nBzMul:\n");
            Bzmultiplied.view();
            PetscPrintf(MPI_COMM_WORLD, "\n");
            }

            SlepcMat combined = Jmultiplied + Bxmultiplied + Bzmultiplied;


            if(verbose){
                PetscPrintf(MPI_COMM_WORLD, "\nCombined Matrix ( J + Bx + Bz ):\n");
                combined.view();
            }

            if(verbose) PetscPrintf(MPI_COMM_WORLD, "\nSetting up return vectors for eigenvalue calculation...\n" );
            Vec imaginary, real;
            PetscScalar im, re;
            if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0mCalculating lowest eigenpair...\n" );

            combined.getLowestEigenpair(&re, &im, &real, &imaginary);

            if(verbose) PetscPrintf(MPI_COMM_WORLD,"\033[0;32mDone!\n\033[0m");
            PetscPrintf(MPI_COMM_WORLD,"Lowest Eigenvalues:\nReal portion: %f\nImaginary Portion: %f\n" , re, im);
            PetscPrintf(MPI_COMM_WORLD, "Real Eigenvector:\n");
            VecView(real, PETSC_VIEWER_STDOUT_WORLD);
            PetscPrintf(MPI_COMM_WORLD, "Imaginary Eigenvector:\n ");
            VecView(imaginary, PETSC_VIEWER_STDOUT_WORLD);
            PetscPrintf(MPI_COMM_WORLD,"\n");
            VecDestroy(&imaginary);
            VecDestroy(&real);
        }



};