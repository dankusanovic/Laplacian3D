//====================================================================================================================================
// IMPLEMENTATION FILE: "setModelData"
//====================================================================================================================================
// Syntax      : setModelData(Coordinates,Elements,GaussPoints,Restraints,Constraints,Dofs,Nnodes,Nelems,Ngauss,Nrestr,Nconst)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Sets FEM Model of analysis. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               GaussPoints : List of Gauss Integration values                   [Ngauss,3]
//               Restraints  : Restrained nodal values                            [Nrestr,3]
//               Dofs        : Free degree of freedom numbering                   [Nnodes,1]
//               Nnodes      : Number of total nodes                              [1,1] 
//               Nelems      : Number of total elements                           [1,1] 
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Coordinates : Updated list of coordinate values                  [Nnodes,3]
//               Elements    : Updated list of element values                     [Nelems,3]
//               GaussPoints : Updated list of Gauss Integration values           [Ngauss,3]
//               Dofs        : Allocated memory for DOFs                          [Nnodes,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Sets model configuration>
// Keywords    : <Coordinates, Elements, GaussPoints>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <cstdlib>
  #include <sstream>
  #include <fstream>
  #include <iostream>
  
   void setModelData(double** &Coordinates, int** &Elements, int** &Boundaries, int** &MeshRefine, double ** &GaussPoints, 
                     int* &Restraints, int &Nnodes, int &Nelems, int &Nfaces, int &Ngauss){ 

//------------------------------------------------------------------------------------------------------------------------------------
// READ THE BOUNDARY FACES: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ifstream INFILEbound("Boundaries.txt"); 

	INFILEbound >> Nfaces;

        int BFaces[Nfaces];

      //Memory Allocation for Boundary Faces:
        for(int i = 0; i < Nfaces; i++){
            INFILEbound >> BFaces[i]; 
        }

     INFILEbound.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// READ COORDINATES INFORMATION: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::ifstream INFILEdata("3d_refined.mesh3d");

	INFILEdata >> Nnodes;

      //Memory Allocation for Coordinates:
        Coordinates = new double* [Nnodes];

        for(int i = 0; i < Nnodes; i++){
            Coordinates[i] = new double[3]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
	    INFILEdata >> Coordinates[i][0] >> Coordinates[i][1] >> Coordinates[i][2];
	} 

//------------------------------------------------------------------------------------------------------------------------------------
// READ ELEMENTS INFORMATION: 
//------------------------------------------------------------------------------------------------------------------------------------
	INFILEdata >> Nelems;

      //Memory Allocation for Elements:
        Elements   = new int* [Nelems];
        MeshRefine = new int* [Nelems];

        for(int i = 0; i < Nelems; i++){
            Elements[i]   = new int[4]; 
            MeshRefine[i] = new int[3];
        }

      //Save Values in Elements:
	for(int i = 0; i < Nelems; i++){
	    INFILEdata >> MeshRefine[i][0] >> MeshRefine[i][1] >> MeshRefine[i][2] 
                       >>   Elements[i][0] >>   Elements[i][1] >>   Elements[i][2] >> Elements[i][3];
  
          //C++ Format index:
            Elements[i][0]--; Elements[i][1]--; Elements[i][2]--; Elements[i][3]--;
	} 

//------------------------------------------------------------------------------------------------------------------------------------
// READ BOUNDARY INFORMATION: 
//------------------------------------------------------------------------------------------------------------------------------------
	INFILEdata >> Nfaces;

      //Memory Allocation for Boundaries:
        Boundaries = new int* [Nfaces];

        for(int i = 0; i < Nfaces; i++){
            Boundaries[i] = new int[4]; 
        }

      //Memory Allocation for Restraints:
        Restraints = new int[Nnodes];

      //Initialization of Restraints Vector:
        for(int i = 0; i < Nnodes; i++){
            Restraints[i] = 0;
        }

      //Save Values in Boundaries:
	for(int i = 0; i < Nfaces; i++){
	    INFILEdata >> Boundaries[i][0] >> Boundaries[i][1] >> Boundaries[i][2] >> Boundaries[i][3];
  
          //C++ Format index:
            Boundaries[i][0]--; Boundaries[i][1]--; Boundaries[i][2]--; Boundaries[i][3]--;

	    switch(BFaces[Boundaries[i][0]]){
		case 1: {
		         Restraints[Boundaries[i][1]] = 1;
		         Restraints[Boundaries[i][2]] = 1;
		         Restraints[Boundaries[i][3]] = 1;
		         break;
			}

		case 2: {
		       //Not yet Implemented.
		         break;
			}
	     }
	} 

   INFILEdata.close();

//------------------------------------------------------------------------------------------------------------------------------------
// READ INTEGRATION DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
   std::ifstream INFILEgauss("GaussPoints.txt"); 

	INFILEgauss >> Ngauss;

      //Memory Allocation for Coordinates:
        GaussPoints = new double* [Ngauss];

        for(int i = 0; i < Ngauss; i++){
            GaussPoints[i] = new double[5]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Ngauss; i++){
	    INFILEgauss >> GaussPoints[i][0] >> GaussPoints[i][1] >> GaussPoints[i][2] >> GaussPoints[i][3] >> GaussPoints[i][4];
	} 

   INFILEgauss.close(); 

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

