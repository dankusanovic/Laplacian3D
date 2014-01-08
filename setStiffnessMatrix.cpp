//====================================================================================================================================
// IMPLEMENTATION FILE: "setStiffnessMatrix"
//====================================================================================================================================
// Syntax      : setStiffnessMatrix(Coordinates, Elements, Stiffness,row,col,Nelems)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Computes the stiffness matrix. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               Stiffness   : List of Stiffness values                           [Nzeros,1]
//               row         : List of indeces i for the K matrix                 [Nelems,1]
//               col         : List of indeces j for the K matrix                 [Nelems,1]
//               Dofs        : Free degree of freedom numbering                   [Nnodes,1]
//               Nelems      : Number of elements                                 [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Stiffness   : Stiffness matrix values                            [Nelems,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Stiffness Matrix>
// Keywords    : <Stiffness,Assembly>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <cmath>

   void setStiffnessMatrix(double** &Coordinates, int** &Elements, double* &Stiffness, int* &row, int* &col, int* &Dofs, int Nelems){
                                     
   //Element properties:
     int    DOFi, DOFj;
     double Volume, nx[4], ny[4], nz[4];
     double ux, uy, uz, vx, vy, vz, wx, wy, wz;

     int count = 0;
     for(int k = 0; k < Nelems; k++){

      //------------------------------------------------------------------------------------------------------------------------------
      // LOCAL NORMAL VECTOR COMPUTATION:
      //------------------------------------------------------------------------------------------------------------------------------
         ux = Coordinates[Elements[k][3]][0] - Coordinates[Elements[k][1]][0];
         uy = Coordinates[Elements[k][3]][1] - Coordinates[Elements[k][1]][1];
         uz = Coordinates[Elements[k][3]][2] - Coordinates[Elements[k][1]][2];

         vx = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][1]][0];
         vy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][1]][1];
         vz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][1]][2];

         nx[0] = (uy*vz - vy*uz)/2; ny[0] = (uz*vx - vz*ux)/2; nz[0] = (ux*vy - vx*uy)/2;

         ux = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][0]][0];
         uy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][0]][1];
         uz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][0]][2];

         vx = Coordinates[Elements[k][3]][0] - Coordinates[Elements[k][0]][0];
         vy = Coordinates[Elements[k][3]][1] - Coordinates[Elements[k][0]][1];
         vz = Coordinates[Elements[k][3]][2] - Coordinates[Elements[k][0]][2];

         nx[1] = (uy*vz - vy*uz)/2; ny[1] = (uz*vx - vz*ux)/2; nz[1] = (ux*vy - vx*uy)/2;

         ux = Coordinates[Elements[k][3]][0] - Coordinates[Elements[k][0]][0];
         uy = Coordinates[Elements[k][3]][1] - Coordinates[Elements[k][0]][1];
         uz = Coordinates[Elements[k][3]][2] - Coordinates[Elements[k][0]][2];

         vx = Coordinates[Elements[k][1]][0] - Coordinates[Elements[k][0]][0];
         vy = Coordinates[Elements[k][1]][1] - Coordinates[Elements[k][0]][1];
         vz = Coordinates[Elements[k][1]][2] - Coordinates[Elements[k][0]][2];

         nx[2] = (uy*vz - vy*uz)/2; ny[2] = (uz*vx - vz*ux)/2; nz[2] = (ux*vy - vx*uy)/2;

         ux = Coordinates[Elements[k][1]][0] - Coordinates[Elements[k][0]][0];
         uy = Coordinates[Elements[k][1]][1] - Coordinates[Elements[k][0]][1];
         uz = Coordinates[Elements[k][1]][2] - Coordinates[Elements[k][0]][2];

         vx = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][0]][0];
         vy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][0]][1];
         vz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][0]][2];

         nx[3] = (uy*vz - vy*uz)/2; ny[3] = (uz*vx - vz*ux)/2; nz[3] = (ux*vy - vx*uy)/2;

      //------------------------------------------------------------------------------------------------------------------------------
      // LOCAL ELEMENT VOLUME COMPUTATION:
      //------------------------------------------------------------------------------------------------------------------------------
         ux = Coordinates[Elements[k][1]][0] - Coordinates[Elements[k][0]][0];
         uy = Coordinates[Elements[k][1]][1] - Coordinates[Elements[k][0]][1];
         uz = Coordinates[Elements[k][1]][2] - Coordinates[Elements[k][0]][2];

         vx = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][0]][0];
         vy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][0]][1];
         vz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][0]][2];

         wx = Coordinates[Elements[k][3]][0] - Coordinates[Elements[k][0]][0];
         wy = Coordinates[Elements[k][3]][1] - Coordinates[Elements[k][0]][1];
         wz = Coordinates[Elements[k][3]][2] - Coordinates[Elements[k][0]][2];

         Volume = 1.5*fabs(ux*(vy*wz - vz*wy) - uy*(vx*wz - vz*wx) + uz*(vx*wy - vy*wx));

      //------------------------------------------------------------------------------------------------------------------------------
      // COMPUTES LOCAL ELEMENT STIFFNESS MATRIX:
      //------------------------------------------------------------------------------------------------------------------------------
         for(int i = 0; i < 4; i++){
             DOFi = Dofs[Elements[k][i]]; 

             for(int j = 0; j < 4; j++){
                 DOFj= Dofs[Elements[k][j]];

                 if(DOFi != -1 &&  DOFj != -1){
	          //COO Stiffness Matix format:
		    row[count]       = DOFi + 1;
		    col[count]       = DOFj + 1;
		    Stiffness[count] = 1.0/Volume*(nx[i]*nx[j] + ny[i]*ny[j] + nz[i]*nz[j]);

		  //Increases index for COO format:
		    count = count + 1;
                 }
             }
         }

     }

   } 

//====================================================================================================================================
// EOF
//====================================================================================================================================

