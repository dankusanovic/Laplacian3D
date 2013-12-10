//====================================================================================================================================
// FINITE ELEMENT PROGRAM FOR THE LAPLACIAN EQUATION
//====================================================================================================================================
// Syntax      : main(ITER,PATH)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Runs the specified model located at PATH. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : ITER        : Number of iterations                               [1,1]
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
// Written by Danilo S. Kusanovic, Dicember 2012
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <string>   
  #include <iostream> 
  #include <vector>


  #include <fstream>
  #include "setAnalysis.hh"
  #include "setModelData.hh"

   int main(int argc, char** argv){

     int ITER;
     std::string PATH = "";

     int Nnodes, Nelems, Ngauss;

     double **Coordinates;
     double **Elements;
     double **GaussPoints;

     setAnalysis(argc,argv,PATH,ITER);
     setModelData(PATH,Coordinates,Elements,Nnodes,Nelems);


//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     int cols;
     std::string COORDINATES = PATH;
     COORDINATES.append("/Coordinates.txt");

     std::ifstream INFILE(COORDINATES.c_str()); 

	INFILE >> Nnodes >> cols;

      //Memory Allocation for Coordinates:
        Coordinates = new double* [Nnodes];

        for(int i = 0; i < Nnodes; i++){
            Coordinates[i] = new double[cols]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
	    INFILE >> Coordinates[i][0] >> Coordinates[i][1] >> Coordinates[i][2];
	} 

     INFILE.close(); 

     for (int i = 0; i < Nnodes; ++i){
	  for (int j = 0; j < 3; ++j)
	       std::cout << Coordinates[i][j] << "\t";
          std::cout << std::endl;
     } 
/*
     for (int i = 0; i < Nelems; ++i){
	  for (int j = 0; j < 4; ++j)
	       std::cout << Elements[i][j] << "\t";
          std::cout << std::endl;
     } */

     std::cout << "Number of Iterations: " << ITER   << std::endl;
     std::cout << "Path of Files       : " << PATH   << std::endl;
     std::cout << "Number of Nodes     : " << Nnodes << std::endl;
     //std::cout << "Number of Elements  : " << Nelems << std::endl;

     for (int i = 0; i < Nnodes; i++){
          delete[] Coordinates[i]; 
     }
     delete[] Coordinates;

     for (int i = 0; i < Nelems; i++){
          delete[] Elements[i]; 
     }
     delete[] Elements;

     //delete[] GaussPoints;

     return 0;

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================
   
