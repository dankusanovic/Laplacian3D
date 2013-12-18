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
  #include <vector>
  #include <iostream> 
  #include "setAnalysis.hh"
  #include "setModelData.hh"
  #include "setForceVector.hh"
//#include "setStiffnessMatrix.hh"
  #include "setBoundaries.hh"
  #include "freeModelData.hh"

   int main(int argc, char** argv){

     int ITER;
     std::string PATH = "";

   //Model Dimension:
     int Nnodes, Nelems, Ngauss, Nrestr, Nconst;

   //Model Variables:
     int    **Elements, **Restraints, **Constraints;
     double **GaussPoints, **Coordinates, *Force;

  //----------------------------------------------------------------------------------------------------------------------------------
  // PRE-ANALYSIS : Reads information from PATH. 
  //----------------------------------------------------------------------------------------------------------------------------------
     setAnalysis(argc,argv,PATH,ITER);
     setModelData(PATH,Coordinates,Elements,GaussPoints,Restraints,Constraints,Force,Nnodes,Nelems,Ngauss,Nrestr,Nconst);

  //----------------------------------------------------------------------------------------------------------------------------------
  // RUN-ANALYSIS : Starts solver. 
  //----------------------------------------------------------------------------------------------------------------------------------
     setForceVector(Coordinates,Elements,GaussPoints,Restraints,Force,Nnodes,Nelems,Ngauss);
   //setStiffnessMatrix(Coordinates,Elements,GaussPoints,Nnodes,Nelems,Ngauss);
     setBoundaries(Coordinates,Restraints,Constraints,Force,Nrestr,Nconst);

  //----------------------------------------------------------------------------------------------------------------------------------
  // POST-ANALYSIS: Save information to PATH. 
  //----------------------------------------------------------------------------------------------------------------------------------
     freeModelData(Coordinates,Elements,GaussPoints,Restraints,Constraints,Force,Nnodes,Nelems,Ngauss,Nrestr,Nconst);

     return 0;

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================
   
