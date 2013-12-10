//====================================================================================================================================
// IMPLEMENTATION FILE: "setAnalysis"
//====================================================================================================================================
// Syntax      : setAnalysis(argc, argv, PATH, ITER)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Sets parameters of analysis. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : ITER        : Number of iterations                               [1,1]
//               PATH        : Path to the folder files                           [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : ITER        : Number of iterations                               [1,1]
//               PATH        : Path to the folder files                           [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Sets parameters>
// Keywords    : <Path, Iteration, Method>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, Dicember 2012
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <string>
  #include <cstdlib>
  #include <sstream>
  #include <iostream>

   void setAnalysis(int argc, char** argv, std::string &PATH, int &ITER){

//------------------------------------------------------------------------------------------------------------------------------------
// INPUTS FILES TO BE LOADED: 
//------------------------------------------------------------------------------------------------------------------------------------
     switch(argc){
        case 1:{
	         ITER = 8;
                 PATH = "./";
                 break;
	       }
        case 2:{
    	         ITER = atoi(argv[1]);
                 PATH = "./";
                 break;
	       }
        case 3:{
	         ITER = atoi(argv[1]);
                 PATH = std::string(argv[2]);
                 break;
	       }

        default:{
                 break;
 	        }
     }

//------------------------------------------------------------------------------------------------------------------------------------
// ANALYSIS TO BE CARRIED OUT: 
//------------------------------------------------------------------------------------------------------------------------------------
     

   }

//====================================================================================================================================
// EOF
//====================================================================================================================================

