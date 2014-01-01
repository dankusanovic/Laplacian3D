//====================================================================================================================================
// IMPLEMENTATION FILE: "setAnalysis"
//====================================================================================================================================
// Syntax      : allocateModelData(Stiffness,Force,row,col,Nfree,Nzeros)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Memory allocation of model variables. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Stiffness   : List of Stiffness values                           [Nzeros,1]
//               Force       : Force vector                                       [Nfree,1]
//               row         : List of indeces i for the K matrix                 [Nzeros,1]
//               col         : List of indeces j for the K matrix                 [Nzeros,1]
//               Nfree       : Number of free nodes                               [1,1] 
//               Nzeros      : Number of nonzero values for stiffness matrix      [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Stiffness   : List of Stiffness values                           [Nzeros,1]
//               Force       : Force vector                                       [Nfree,1]
//               row         : List of indeces i for the K matrix                 [Nzeros,1]
//               col         : List of indeces j for the K matrix                 [Nzeros,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Sets model memory allocation>
// Keywords    : <stiffness, force, memory allocation>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, January 2014
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <cstdlib>
  
   void allocateModelData(double* &Stiffness, double* &Force, int* &row, int* &col, int Nfree, int Nzeros){ 

//------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATES STIFFNESS MATRIX DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     row       = new int[Nzeros];
     col       = new int[Nzeros];
     Stiffness = new double[Nzeros];

//------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATES FORCE VECTOR DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     Force = new double[Nfree];

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

