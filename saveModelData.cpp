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

  #include <cstdlib>
  #include <sstream>
  #include <fstream>
  #include <iostream>
  
   void saveModelData(double** &REFINE, int ITER, double** &Coordinates, int** &Elements, double* &Force, int* &Dofs, int Nnodes, int Nelems){ 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ofstream OUTFILEcoord("Results/Coordinates.txt"); 
      
	OUTFILEcoord << Nnodes << std::endl;

      //Save Values in Coordinates:
        OUTFILEcoord.precision(20);
	for(int i = 0; i < Nnodes; i++){
	    OUTFILEcoord << Coordinates[i][0] << "\t" << Coordinates[i][1] << "\t" << Coordinates[i][2] << std::endl;
	} 

     OUTFILEcoord.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE ELEMENTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ofstream OUTFILEelems("Results/Elements.txt"); 

	OUTFILEelems << Nelems << std::endl;

      //Save Values in Coordinates:
	for(int i = 0; i < Nelems; i++){
	    OUTFILEelems << Elements[i][0] + 1 << "\t" << Elements[i][1] + 1 << "\t" << Elements[i][2] + 1 << "\t" 
                         << Elements[i][3] + 1 << std::endl;
	} 

     OUTFILEelems.close();

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE SOLUTION DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ofstream OUTFILEsolution("Results/Solution.txt"); 
        OUTFILEsolution.precision(20);

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

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE CONVERGENCE DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ofstream OUTFILEconvergence("Results/Convergence.txt"); 
        OUTFILEconvergence.precision(20);

      //Save Values in Coordinates:
	for(int i = 0; i < ITER; i++){
            OUTFILEconvergence << "Iteration : " << REFINE[i][0] << "\t NDofs: " << REFINE[i][1] << "\t Error: " << REFINE[i][2] << std::endl;
	} 

     OUTFILEconvergence.close();  

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

