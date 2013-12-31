//====================================================================================================================================
// FINITE ELEMENT PROGRAM FOR THE LAPLACE EQUATION
//====================================================================================================================================
// Syntax      : Main(ITER,PATH)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Runs the specified model located at PATH. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : ITER        : Number of iterations                               [1,1]
//               DOFS        : Number of degree of freedom per node               [1,1]
//               PATH        : Path to the folder files                           [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Solution    : Nodal values                                       [nElem,1]
//               Mesh        : Vector of elements                                 [nElem,4]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Performs an analysis>
// Keywords    : <Main, Model>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <string>   
  #include <iostream> 
  #include "MUMPSSolver.hh"
  #include "setAnalysis.hh"
  #include "setModelData.hh"
//#include "saveModelData.hh"
  #include "freeModelData.hh"
  #include "setBoundaries.hh"
  #include "setForceVector.hh"
  #include "setStiffnessMatrix.hh"

   int main(int argc, char** argv){

     int ITER;
     std::string PATH = "";
     setAnalysis(argc,argv,PATH,ITER);

     for(int k = 0; k < ITER; k++){

       //Model Dimension:
         int Nnodes, Nelems, Ngauss, Nrestr, Nconst, Nzeros;

       //Model Variables:
         int    **Elements, **Restraints, **Constraints, *row, *col;
         double **GaussPoints, **Coordinates, *Stiffness, *Force;

      //------------------------------------------------------------------------------------------------------------------------------
      // PRE-ANALYSIS :  
      //------------------------------------------------------------------------------------------------------------------------------
       //Reads information:
         setModelData(PATH,Coordinates,Elements,GaussPoints,Restraints,Constraints,Stiffness,Force,row,col,Nnodes,Nelems,Ngauss,Nrestr,Nconst,Nzeros);

      //------------------------------------------------------------------------------------------------------------------------------
      // RUN-ANALYSIS :  
      //------------------------------------------------------------------------------------------------------------------------------
       //FEM Assembly:
         setForceVector(Coordinates,Elements,GaussPoints,Force,Nnodes,Nelems,Ngauss);
         setStiffnessMatrix(Coordinates,Elements,Stiffness,row,col,Nelems);
         setBoundaries(Coordinates,Restraints,Constraints,Force,Nrestr,Nconst);

       //FEM Solution:
         MUMPSSolver(Nnodes,Nzeros,row,col,Stiffness, Force,0);

      //------------------------------------------------------------------------------------------------------------------------------
      // POST-ANALYSIS: 
      //------------------------------------------------------------------------------------------------------------------------------
       //Save solution:
       //saveModelData(Coordinates,Elements,Restraints,Constraints,Force,Nnodes,Nelems,Nrestr,Nconst);

       //Free memory:
         freeModelData(Coordinates,Elements,GaussPoints,Restraints,Constraints,Stiffness,Force,row,col,Nnodes,Nelems,Ngauss,Nrestr,Nconst);

     }

     return 0;

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================
   
