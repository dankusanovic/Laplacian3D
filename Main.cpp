//====================================================================================================================================
// FINITE ELEMENT PROGRAM FOR THE LAPLACE EQUATION
//====================================================================================================================================
// Syntax      : Main(ITER,PATH)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Runs the specified model located at PATH. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : ITER        : Number of iterations                               [1,1]
//               PATH        : Path to the folder files                           [1,1]
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

  #include <string>   
  #include <iostream> 
  #include "MUMPSSolver.hh"
  #include "setAnalysis.hh"
  #include "setModelData.hh"
  #include "setModelDofs.hh"
  #include "saveModelData.hh"
  #include "freeModelData.hh"
  #include "setForceVector.hh"
  #include "setRealSolution.hh"
  #include "allocateModelData.hh"
  #include "setStiffnessMatrix.hh"

   int main(int argc, char** argv){

     int ITER;
     std::string PATH = "";
     setAnalysis(argc,argv,PATH,ITER);

     for(int k = 0; k < ITER; k++){

       //Model Dimension:
         double TotalError;
         int Nnodes, Nelems, Ngauss, Nrestr, Nconst, Nfree, Nzeros;

       //Model Variables:
         int    **Elements, **Restraints, **Constraints, *Dofs, *row, *col;
         double **GaussPoints, **Coordinates, *Stiffness, *Force; 

      //------------------------------------------------------------------------------------------------------------------------------
      // PRE-ANALYSIS :  
      //------------------------------------------------------------------------------------------------------------------------------
       //Reads information:
          setModelData(PATH,Coordinates,Elements,GaussPoints,Restraints,Constraints,Dofs,Nnodes,Nelems,Ngauss,Nrestr,Nconst);
          setModelDofs(Elements,Restraints,Dofs,Nnodes,Nelems,Nrestr,Nfree,Nzeros);
          allocateModelData(Stiffness,Force,row,col,Nfree,Nzeros);

      //------------------------------------------------------------------------------------------------------------------------------
      // RUN-ANALYSIS :  
      //------------------------------------------------------------------------------------------------------------------------------
       //FEM Assembly:
         setForceVector(Coordinates,Elements,GaussPoints,Force,Dofs,Nnodes,Nelems,Ngauss);
         setStiffnessMatrix(Coordinates,Elements,Stiffness,row,col,Dofs,Nelems);

       //FEM Solution:
         MUMPSSolver(Nfree,Nzeros,row,col,Stiffness, Force,0);

      //------------------------------------------------------------------------------------------------------------------------------
      // POST-ANALYSIS: 
      //------------------------------------------------------------------------------------------------------------------------------
       //Exact Solution:
         setRealSolution(Coordinates,Elements,GaussPoints,Force,TotalError,Nelems,Ngauss);

	 std::cout << std::endl;
	 std::cout << "nDOfs: "<< Nfree << " Error: " << TotalError << std::endl;

       //Store Solution:
         saveModelData(PATH,Coordinates,Elements,Restraints,Constraints,Force,Dofs,Nnodes,Nelems,Nrestr,Nconst);

       //Mesh Refiner:
         //system("./Refinement");

       //Free Memory:
         freeModelData(Coordinates,Elements,GaussPoints,Restraints,Constraints,Dofs,Stiffness,Force,row,col,Nnodes,Nelems,Ngauss,Nrestr,Nconst);

     }

     return 0;

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================
   
