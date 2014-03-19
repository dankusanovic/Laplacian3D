//====================================================================================================================================
// FINITE ELEMENT PROGRAM FOR THE LAPLACE EQUATION
//====================================================================================================================================
// Syntax      : Main(ITER,PATH)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Runs the specified model located at PATH.
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : ITER        : Number of iterations                               [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : REFINE      : Number of iterations and Error                     [ITER,3]
//               Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               Solution    : Nodal Solution values                              [Nnodes,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      :
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Performs an analysis>
// Keywords    : <Main, Model, Program>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

#include "setAnalysis.hh"
#include "MUMPSSolver.hh"
#include "setModelData.hh"
#include "setModelDofs.hh"
#include "setTotalError.hh"
#include "saveModelData.hh"
#include "allocateModel.hh"
#include "setForceVector.hh"
#include "setMeshRefiner.hh"
#include "setStiffnessMatrix.hh"
#include <cstdlib>
#include <iostream>
int main(int argc, char** argv) {

    int ITER;
    double **REFINE;

    setAnalysis(argc,argv,ITER,REFINE);

    for(int i = 0; i < ITER; i++) {

        //Model Dimension:
        int Nnodes, Nelems, Nfaces, Ngauss, Nfree, Nzeros;

        //Model Variables:
        int    **Elements, **vecK, **Boundaries, **MeshRefine, *Restraints, *Dofs, *row, *col;
        double **GaussPoints, **Coordinates, *Stiffness, *Force, *etaK, etamax;

        //------------------------------------------------------------------------------------------------------------------------------
        // PRE-ANALYSIS :
        //------------------------------------------------------------------------------------------------------------------------------
        //Reads information:
        setModelData(Coordinates,Elements,vecK,Boundaries,MeshRefine,GaussPoints,Restraints,Nnodes,Nelems,Nfaces,Ngauss,etaK);
        setModelDofs(Elements,Restraints,Dofs,Nnodes,Nelems,Nfree,Nzeros);
        allocateModel(Stiffness,Force,row,col,Nfree,Nzeros);

        //------------------------------------------------------------------------------------------------------------------------------
        // RUN-ANALYSIS :
        //------------------------------------------------------------------------------------------------------------------------------
        //FEM Assembly:
        setForceVector(Coordinates,Elements,GaussPoints,Force,Dofs,Nelems,Ngauss);
        setStiffnessMatrix(Coordinates,Elements,Stiffness,row,col,Dofs,Nelems);

        //FEM Solution:
        MUMPSSolver(Nfree,Nzeros,row,col,Stiffness,Force,0);

        //------------------------------------------------------------------------------------------------------------------------------
        // POST-ANALYSIS:
        //------------------------------------------------------------------------------------------------------------------------------
        //Error Approximation:
        etamax = setTotalError(REFINE,Coordinates,Elements,GaussPoints,Force,Dofs,Nnodes,Nelems,Ngauss,i,etaK,vecK);

        //Mark the refining elemenents in MeshRefine, based on the etaK
        //markRefine(MeshRefine,etaK,Nelems,etamax);
        for (int j=0; j<Nelems; j++) {
            if (etaK[j]>=0.5*etamax) {
                MeshRefine[j][0]=1;
            }
            else {
                MeshRefine[j][0]=0;
            }
        }

        //Save the Solution:
        saveModelData(REFINE,ITER,Coordinates,Elements,Force,Dofs,Nnodes,Nelems);

        //Mesh Refinement:
        if(i < ITER - 1)
            setMeshRefiner(Coordinates,Elements,MeshRefine,Boundaries,Nnodes,Nelems,Nfaces);
        //Free Memory:
        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE DOFS VECTOR DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        delete[] Dofs;
        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE STIFFNESS MATRIX DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        delete[] row;
        delete[] col;
        delete[] Stiffness;
        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE FORCE VECTOR DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        delete[] Force;
        //freeModelData(Coordinates,Elements,GaussPoints,MeshRefine,Boundaries,Restraints,Nnodes,Nelems,Nfaces,Ngauss);
        //for (int i=0;i<Nelems;i++){
        //std::cout<< i << "\t"<<vecK[i][0]<< "\t"<<vecK[i][1]<< "\t"<<vecK[i][2]<< "\t"<<vecK[i][3]<<std::endl;
        //}

        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE COORDINATES DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        for(int i = 0; i < Nnodes; i++) {
            delete[] Coordinates[i];
        }
        delete[] Coordinates;

        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE ELEMENTS DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        for(int i = 0; i < Nelems; i++) {
            delete[] Elements[i];
        }
        delete[] Elements;

        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE GAUSS POINTS DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        for(int i = 0; i < Ngauss; i++) {
            delete[] GaussPoints[i];
        }
        delete[] GaussPoints;

        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE REFINEMENT DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        for(int i = 0; i < Nelems; i++) {
            delete[] MeshRefine[i];
        }
        delete[] MeshRefine;

        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE BOUNDARY DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        for(int i = 0; i < Nfaces; i++) {
            delete[] Boundaries[i];
        }
        delete[] Boundaries;

        //------------------------------------------------------------------------------------------------------------------------------------
        // FREE RESTRAINT DATA:
        //------------------------------------------------------------------------------------------------------------------------------------
        delete[] Restraints;

        delete vecK;
        delete etaK;
    }

    return 0;

}

//====================================================================================================================================
// EOF
//====================================================================================================================================

