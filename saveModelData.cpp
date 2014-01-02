//====================================================================================================================================
// IMPLEMENTATION FILE: "setModelData"
//====================================================================================================================================
// Syntax      : saveModelData(PATH,Coordinates,Elements,GaussPoints,Restraints,Constraints,Dofs,Nnodes,Nelems,Ngauss,Nrestr,Nconst)
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
// Written by Danilo S. Kusanovic, January 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <string>
  #include <cstdlib>
  #include <sstream>
  #include <fstream>
  #include <iostream>
  
   void saveModelData(){ 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string COORDINATES = PATH;
     COORDINATES.append("/Coordinates.txt");

     std::ofstream OUTFILEcoord(COORDINATES.c_str()); 

	OUTFILEcoord << Nnodes << 3 << std::endl;

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
	    OUTFILEcoord << Coordinates[i][0] << " " << Coordinates[i][1] << " " << Coordinates[i][2] << std::endl;
	} 

     OUTFILEcoord.close(); 

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

