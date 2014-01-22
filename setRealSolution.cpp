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
  #include <iostream> 

//------------------------------------------------------------------------------------------------------------------------------------
// ANALYTICAL FLUX FUNCTION: 
//------------------------------------------------------------------------------------------------------------------------------------
   double get_ux(double x, double y, double z){

     double ux = -y*z*(2*x - 1)*(y - 1)*(z - 1);

     return ux;
   }

   double get_uy(double x, double y, double z){

     double uy = -x*z*(2*y - 1)*(x - 1)*(z - 1);

     return uy;
   }

   double get_uz(double x, double y, double z){

     double uz = -x*y*(2*z - 1)*(x - 1)*(y - 1);

     return uz;
   }

//------------------------------------------------------------------------------------------------------------------------------------
// COMPUTES THE ERROR (L2-NORM): 
//------------------------------------------------------------------------------------------------------------------------------------
   void setRealSolution(double** &Coordinates, int** &Elements, double** &GaussPoints, double* &uh, double & TotalError, int Nelems, int Ngauss){
                                     
   //Element properties:
     double Volume, nx[4], ny[4], nz[4];
     double ux, uy, uz, vx, vy, vz, wx, wy, wz;
     double uhx, uhy, uhz, xi, eta, zeta, error;

     TotalError = 0.0;

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

         Volume = fabs(ux*(vy*wz - vz*wy) - uy*(vx*wz - vz*wx) + uz*(vx*wy - vy*wx))/6;

      //------------------------------------------------------------------------------------------------------------------------------
      // COMPUTES THE ERROR:
      //------------------------------------------------------------------------------------------------------------------------------       
	 error = 0.0; 

       //Numerical Flux Values:
	 uhx   = (-nx[0]*uh[Elements[k][0]] - nx[1]*uh[Elements[k][1]] - nx[2]*uh[Elements[k][2]] - nx[3]*uh[Elements[k][3]])/(3*Volume);
	 uhy   = (-ny[0]*uh[Elements[k][0]] - ny[1]*uh[Elements[k][1]] - ny[2]*uh[Elements[k][2]] - ny[3]*uh[Elements[k][3]])/(3*Volume);
	 uhz   = (-nz[0]*uh[Elements[k][0]] - nz[1]*uh[Elements[k][1]] - nz[2]*uh[Elements[k][2]] - nz[3]*uh[Elements[k][3]])/(3*Volume);

         for(int j = 0; j < Ngauss; j++){

           //Transformed triangular Coordinates:
	     xi   = Coordinates[Elements[k][0]][0]*GaussPoints[j][1] + Coordinates[Elements[k][1]][0]*GaussPoints[j][2] + 
                    Coordinates[Elements[k][2]][0]*GaussPoints[j][3] + Coordinates[Elements[k][3]][0]*GaussPoints[j][4];
             eta  = Coordinates[Elements[k][0]][1]*GaussPoints[j][1] + Coordinates[Elements[k][1]][1]*GaussPoints[j][2] + 
                    Coordinates[Elements[k][2]][1]*GaussPoints[j][3] + Coordinates[Elements[k][3]][1]*GaussPoints[j][4];
             zeta = Coordinates[Elements[k][0]][2]*GaussPoints[j][1] + Coordinates[Elements[k][1]][2]*GaussPoints[j][2] + 
                    Coordinates[Elements[k][2]][2]*GaussPoints[j][3] + Coordinates[Elements[k][3]][2]*GaussPoints[j][4];

           //Error function evaluated at each triangular coordinates:
	     error = error + GaussPoints[j][0]*(pow(get_ux(xi,eta,zeta) - uhx,2) + pow(get_uy(xi,eta,zeta) - uhy,2) + pow(get_uz(xi,eta,zeta) - uhz,2));

	 }

         TotalError  = TotalError + Volume*error;

     }

     TotalError = sqrt(fabs(TotalError));

   } 

//====================================================================================================================================
// EOF
//====================================================================================================================================

