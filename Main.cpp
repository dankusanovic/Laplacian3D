//====================================================================================================================================
// FINITE ELEMENT PROGRAM FOR THE LAPLACE EQUATION
//====================================================================================================================================
// Syntax      : Main(ITER,PATH)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Runs the specified model located at PATH. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : ITER        : Number of iterations                               [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Coordinates : List of coordinate values                          [Nnodes,3]
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

  #include <iostream> 
  #include "MUMPSSolver.hh"
  #include "setAnalysis.hh"
  #include "setModelData.hh"
  #include "setModelDofs.hh"
  #include "setTotalError.hh"
  #include "saveModelData.hh"
  #include "freeModelData.hh"
  #include "setForceVector.hh"
  #include "setMeshRefiner.hh" 
  #include "allocateModelData.hh"
  #include "setStiffnessMatrix.hh"

   int main(int argc, char** argv){

     int ITER;
     double **REFINE;
     setAnalysis(argc,argv,ITER,REFINE);

     for(int k = 0; k < ITER; k++){

       //Model Dimension:
         int Nnodes, Nelems, Nfaces, Ngauss, Nfree, Nzeros;

       //Model Variables:  
         int    **Elements, **Boundaries, **MeshRefine, *Restraints, *Dofs, *row, *col;
         double **GaussPoints, **Coordinates, *Stiffness, *Force; 

      //------------------------------------------------------------------------------------------------------------------------------
      // PRE-ANALYSIS :  
      //------------------------------------------------------------------------------------------------------------------------------
       //Reads information:
         setModelData(Coordinates,Elements,Boundaries,MeshRefine,GaussPoints,Restraints,Nnodes,Nelems,Nfaces,Ngauss);
         setModelDofs(Elements,Restraints,Dofs,Nnodes,Nelems,Nfree,Nzeros);
         allocateModelData(Stiffness,Force,row,col,Nfree,Nzeros);

      //------------------------------------------------------------------------------------------------------------------------------
      // RUN-ANALYSIS :  
      //------------------------------------------------------------------------------------------------------------------------------
       //FEM Assembly:
         setForceVector(Coordinates,Elements,GaussPoints,Force,Dofs,Nnodes,Nelems,Ngauss);
         setStiffnessMatrix(Coordinates,Elements,Stiffness,row,col,Dofs,Nelems);

       //FEM Solution:
         MUMPSSolver(Nfree,Nzeros,row,col,Stiffness,Force,0);

      //------------------------------------------------------------------------------------------------------------------------------
      // POST-ANALYSIS: 
      //------------------------------------------------------------------------------------------------------------------------------
       //Exact Solution:
         setTotalError(REFINE,Coordinates,Elements,GaussPoints,Force,Dofs,Nnodes,Nelems,Ngauss,k);

       //Store Solution:
         saveModelData(REFINE,ITER,Coordinates,Elements,Force,Dofs,Nnodes,Nelems);

       //Mesh Refinement:
         setMeshRefiner(Coordinates,Elements,MeshRefine,Boundaries,Nnodes,Nelems,Nfaces);

       //Free Memory:
         freeModelData(Coordinates,Elements,GaussPoints,MeshRefine,Boundaries,Restraints,Dofs,Stiffness,Force,row,col,Nnodes,Nelems,Nfaces,Ngauss);

     }

     return 0;

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================
   
