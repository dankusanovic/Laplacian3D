//====================================================================================================================================
// IMPLEMENTATION FILE: "saveModelData"
//====================================================================================================================================
// Syntax      : saveModelData(PATH,Coordinates,Elements,GaussPoints,Restraints,Constraints,Dofs,Nnodes,Nelems,Ngauss,Nrestr,Nconst)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Sets FEM Model of analysis. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : PATH        : Path to the folder files                           [1,1]
//               Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               GaussPoints : List of Gauss Integration values                   [Ngauss,3]
//               Restraints  : Restrained nodal values                            [Nrestr,3]
//               Force       : Force vector values                                [Nnodes,1]
//               Dofs        : Free degree of freedom numbering                   [Nnodes,1]
//               Nnodes      : Number of total nodes                              [1,1] 
//               Nelems      : Number of total elements                           [1,1] 
//               Nrestr      : Number of restrained nodes                         [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Coordinates : Text file with coordinate values                  [Nnodes,3]
//               Elements    : Text file with element values                     [Nelems,4]
//               Solution    : Text file with solution values                    [Nnodes,1]
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
  
   void saveModelData(std::string PATH, double** &Coordinates, int** &Elements, int** &Restraints, int** &Constraints, double* &Force, 
                      int* &Dofs, int Nnodes, int Nelems, int Nrestr, int Nconst){ 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string COORDINATES = PATH;
     COORDINATES.append("Results/Coordinates.txt");

     std::ofstream OUTFILEcoord(COORDINATES.c_str()); 
        OUTFILEcoord.precision(16);

	OUTFILEcoord << Nnodes << " " << 3 << std::endl;

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
	    OUTFILEcoord << Coordinates[i][0] << "\t" << Coordinates[i][1] << "\t" << Coordinates[i][2] << std::endl;
	} 

     OUTFILEcoord.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE ELEMENTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string ELEMENTS = PATH;
     ELEMENTS.append("Results/Elements.txt");

     std::ofstream OUTFILEelems(ELEMENTS.c_str()); 

	OUTFILEcoord << Nelems << " " << 4 << std::endl;

      //Save Values in Coordinates:
	for(int i = 0; i < Nelems; i++){
	    OUTFILEelems << Elements[i][0] << "\t" << Elements[i][1] << "\t" << Elements[i][2] << "\t" << Elements[i][3] << std::endl;
	} 

     OUTFILEelems.close();

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE SOLUTION DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string SOLUTION = PATH;
     SOLUTION.append("Results/Solution.txt");

     std::ofstream OUTFILEsolution(SOLUTION.c_str()); 
        OUTFILEsolution.precision(16);

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
            if(Dofs[i] != -1){
	       OUTFILEsolution << Force[Dofs[i]] << std::endl;
            }
            else{
               OUTFILEsolution << 0 << std::endl;
            }
	} 

     OUTFILEsolution.close();  

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

