//====================================================================================================================================
// IMPLEMENTATION FILE: "setMeshRefiner"
//====================================================================================================================================
// Syntax      : setMeshRefiner(Coordinates,Elements,MeshRefine,Boundaries,Nnodes,Nelems,Nfaces);
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Creates a File of Elements to be Refined. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               MeshRefine  : List of Mesh Information                           [Nelems,3]
//               Boundaries  : List of Restrained Nodes                           [Nfaces,3]
//               Nnodes      : Number of total nodes                              [1,1] 
//               Nelems      : Number of total elements                           [1,1] 
//               Nfaces      : Number of total restrained faces                   [1,1] 
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Coordinates : Text file with coordinate values                  [Nnodes,3]
//               Elements    : Text file with element values                     [Nelems,4]
//               Boundaries  : Text file with Restrained values                  [Nfaces,4]
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
  
   void setMeshRefiner(double** &Coordinates, int** &Elements, int** &MeshRefine, int** &Boundaries, int Nnodes, int Nelems, int Nfaces){ 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ofstream OUTFILE("3d_original.mesh3d"); 

	OUTFILE << Nnodes << std::endl;

      //Save Values in Coordinates:
        OUTFILE.precision(20);

	for(int i = 0; i < Nnodes; i++){
	    OUTFILE << Coordinates[i][0] << "\t" << Coordinates[i][1] << "\t" << Coordinates[i][2] << std::endl;
	} 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE ELEMENTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
	OUTFILE << Nelems << std::endl;

      //Save Values in Coordinates:
	for(int i = 0; i < Nelems; i++){
	    OUTFILE <<   1  
                    << "\t" << MeshRefine[i][1]   << "\t" << MeshRefine[i][2]   << "\t" << Elements[i][0] + 1
                    << "\t" << Elements[i][1] + 1 << "\t" << Elements[i][2] + 1 << "\t" << Elements[i][3] + 1 << std::endl;
	} 

//------------------------------------------------------------------------------------------------------------------------------------
// SAVE SOLUTION DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
        OUTFILE << Nfaces << std::endl;

      //Save Values in Coordinates:
	for(int i = 0; i < Nfaces; i++){
            OUTFILE << Boundaries[i][0] + 1 << "\t" << Boundaries[i][1] + 1 << "\t" 
                    << Boundaries[i][2] + 1 << "\t" << Boundaries[i][3] + 1 << std::endl;
        }

     OUTFILE.close();  

//------------------------------------------------------------------------------------------------------------------------------------
// MESH REFINEMENT PROCESS: 
//------------------------------------------------------------------------------------------------------------------------------------
     system("./Refinement");

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

