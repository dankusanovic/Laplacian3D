//====================================================================================================================================
// IMPLEMENTATION FILE: "setAnalysis"
//====================================================================================================================================
// Syntax      : setModelData(PATH, Coordinates, Elements, GaussPoints)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Sets FEM Model of analysis. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : PATH        : Path to the folder files                           [1,1]
//               Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               GaussPoints : List of Gauss Integration values                   [Ngauss,3]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Coordinates : Updated list of coordinate values                  [Nnodes,3]
//               Elements    : Updated list of element values                     [Nelems,3]
//               GaussPoints : Updated list of Gauss Integration values           [Ngauss,3]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Sets model configuration>
// Keywords    : <Coordinates, Elements, GaussPoints>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, Dicember 2012
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <string>
  #include <cstdlib>
  #include <sstream>
  #include <iostream>
  #include <fstream>

   void setModelData(std::string PATH, double **Coordinates, double **Elements, int &Nnodes, int &Nelems){

     int cols;

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string COORDINATES = PATH;
     COORDINATES.append("/Coordinates.txt");

     std::ifstream INFILEa(COORDINATES.c_str()); 

	INFILEa >> Nnodes >> cols;

      //Memory Allocation for Coordinates:
        Coordinates = new double* [Nnodes];

        for(int i = 0; i < Nnodes; i++){
            Coordinates[i] = new double[cols]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
	    INFILEa >> Coordinates[i][0] >> Coordinates[i][1] >> Coordinates[i][2];
	} 

     INFILEa.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE ELEMENTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string ELEMENTS = PATH;
     ELEMENTS.append("/Elements.txt");

     std::ifstream INFILEb(ELEMENTS.c_str()); 

	INFILEb >> Nelems >> cols;

      //Memory Allocation for Coordinates:
        Elements = new double* [Nelems];

        for(int i = 0; i < Nelems; i++){
            Elements[i] = new double[cols]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Nelems; i++){
	    INFILEb >> Elements[i][0] >> Elements[i][1] >> Elements[i][2] >> Elements[i][3];
	} 

     INFILEb.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE GAUSS POINTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

