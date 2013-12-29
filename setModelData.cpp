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
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <string>
  #include <cstdlib>
  #include <sstream>
  #include <fstream>
  #include <iostream>
  
   void setModelData(std::string PATH, double** &Coordinates, int** &Elements, double ** &GaussPoints, int** &Restraints, int** &Constraints, 
                     double* &Stiffness, double* &Force, int* &row, int* &col, int &Nnodes, int &Nelems, int &Ngauss, int &Nrestr, int &Nconst, int &Nzeros){ 

     int cols;

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE COORDINATES DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string COORDINATES = PATH;
     COORDINATES.append("/Coordinates.txt");

     std::ifstream INFILEcoord(COORDINATES.c_str()); 

	INFILEcoord >> Nnodes >> cols;

      //Memory Allocation for Coordinates:
        Coordinates = new double* [Nnodes];

        for(int i = 0; i < Nnodes; i++){
            Coordinates[i] = new double[cols]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Nnodes; i++){
	    INFILEcoord >> Coordinates[i][0] >> Coordinates[i][1] >> Coordinates[i][2];
	} 

     INFILEcoord.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE ELEMENTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string ELEMENTS = PATH;
     ELEMENTS.append("/Elements.txt");

     std::ifstream INFILEelems(ELEMENTS.c_str()); 

	INFILEelems >> Nelems >> cols;

      //Memory Allocation for Coordinates:
        Elements = new int* [Nelems];

        for(int i = 0; i < Nelems; i++){
            Elements[i] = new int[cols]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Nelems; i++){
	    INFILEelems >> Elements[i][0] >> Elements[i][1] >> Elements[i][2] >> Elements[i][3];
	} 

     INFILEelems.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE GAUSS POINTS DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     std::string INTEGRATION = PATH;
     INTEGRATION.append("/GaussPoints.txt");

     std::ifstream INFILEgauss(INTEGRATION.c_str()); 

	INFILEgauss >> Ngauss >> cols;

      //Memory Allocation for Coordinates:
        GaussPoints = new double* [Ngauss];

        for(int i = 0; i < Ngauss; i++){
            GaussPoints[i] = new double[cols]; 
        }

      //Save Values in Coordinates:
	for(int i = 0; i < Ngauss; i++){
	    INFILEgauss >> GaussPoints[i][0] >> GaussPoints[i][1] >> GaussPoints[i][2] >> GaussPoints[i][3] >> GaussPoints[i][4];
	} 

     INFILEgauss.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// READ AND SAVE RESTRAINT AND CONSTRAINT DATA:
//------------------------------------------------------------------------------------------------------------------------------------
     std::string BOUNDARY = PATH;
     BOUNDARY.append("/Boundary.txt");

     std::ifstream INFILEboundary(BOUNDARY.c_str()); 

	INFILEboundary >> Nrestr >> Nconst >> cols;

      //Memory Allocation for Restraints:
        Restraints = new int* [Nrestr];

        for(int i = 0; i < Nrestr; i++){
            Restraints[i] = new int[cols]; 
        }

      //Memory Allocation for Constraints:
        Constraints = new int* [Nconst];

        for(int i = 0; i < Nconst; i++){
            Constraints[i] = new int[cols]; 
        }

      //Differentiate between Restraints (1) and Constraints (2): 
        int j = 0, k = 0, Nboundary = Nrestr + Nconst;

        for(int i = 0; i < Nboundary; i++){
            int cond, dumm;
            INFILEboundary >> cond;

            if(cond == 1){
               INFILEboundary >> dumm >> Restraints[j][0] >> Restraints[j][1] >> Restraints[j][2];
               j = j + 1;
            }
            else if(cond == 2){ 
               INFILEboundary >> dumm >> Constraints[k][0] >> Constraints[k][1] >> Constraints[k][2];
               k = k + 1;
            }
        }

     INFILEboundary.close(); 

//------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATES STIFFNESS MATRIX DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     Nzeros    = 16*Nelems;
     row       = new int[Nzeros];
     col       = new int[Nzeros];
     Stiffness = new double[Nzeros];

//------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATES FORCE VECTOR DATA: 
//------------------------------------------------------------------------------------------------------------------------------------
     Force = new double[Nelems];

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

