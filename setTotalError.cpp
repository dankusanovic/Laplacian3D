//====================================================================================================================================
// IMPLEMENTATION FILE: "setTotalError"
//====================================================================================================================================
// Syntax      : setTotalError(Coordinates, Elements, Stiffness,row,col,Nelems)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Computes the Total error due to the approximation. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               GaussPoints : List of Gauss Integration values                   [Ngauss,3]
//               Dofs        : Free degree of freedom numbering                   [Nnodes,1]
//               Nelems      : Number of elements                                 [1,1]
//               Ngauss      : Number of integration points                       [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Error       : Total error made during the approximation          [1,1]
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
  #include <cstdlib>
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
   void setTotalError(double** &REFINE,double** &Coordinates, int** &Elements, double** &GaussPoints, double* &Solution, int* &Dofs,
                      int Nnodes, int Nelems, int Ngauss, int i){

     double Error = 0.0;

   //Total Vector Solution:
     int  DOFi;
     double *uh;
     uh   = new double[Nnodes];

     for(int j = 0; j < Nnodes; j++){
         DOFi = Dofs[j];

         if(DOFi != -1){
            uh[j] = Solution[DOFi];
         }
         else {
            uh[j] = 0.0;
         } 
     }
                                     
   //Element properties:
     double Volume, nx[4], ny[4], nz[4];
     double ux, uy, uz, vx, vy, vz, wx, wy, wz;
     double uhx, uhy, uhz, xi, eta, zeta, error;

     for(int k = 0; k < Nelems; k++){
         
         error = 0.0; 

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
	 ux = Coordinates[Elements[k][1]][0] - Coordinates[Elements[k][3]][0];
	 uy = Coordinates[Elements[k][1]][1] - Coordinates[Elements[k][3]][1];
	 uz = Coordinates[Elements[k][1]][2] - Coordinates[Elements[k][3]][2];

	 vx = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][3]][0];
	 vy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][3]][1];
	 vz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][3]][2];

         wx = Coordinates[Elements[k][0]][0] - Coordinates[Elements[k][3]][0];
         wy = Coordinates[Elements[k][0]][1] - Coordinates[Elements[k][3]][1];
         wz = Coordinates[Elements[k][0]][2] - Coordinates[Elements[k][3]][2];

         Volume = (wx*(uy*vz - vy*uz) + wy*(uz*vx - vz*ux) + wz*(ux*vy - vx*uy))/6;

      //------------------------------------------------------------------------------------------------------------------------------
      // COMPUTES THE ERROR:
      //------------------------------------------------------------------------------------------------------------------------------      
          
       //Numerical Flux Values:
	 uhx = (-nx[0]*uh[Elements[k][0]] - nx[1]*uh[Elements[k][1]] - nx[2]*uh[Elements[k][2]] - nx[3]*uh[Elements[k][3]])/(3*Volume);
	 uhy = (-ny[0]*uh[Elements[k][0]] - ny[1]*uh[Elements[k][1]] - ny[2]*uh[Elements[k][2]] - ny[3]*uh[Elements[k][3]])/(3*Volume);
	 uhz = (-nz[0]*uh[Elements[k][0]] - nz[1]*uh[Elements[k][1]] - nz[2]*uh[Elements[k][2]] - nz[3]*uh[Elements[k][3]])/(3*Volume);

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

         Error  = Error + Volume*error;

     }

     delete[] uh;
     Error = sqrt(fabs(Error));

   //L2 Norm Information:
     REFINE[i][0] = (double) i + 1; 
     REFINE[i][1] = (double) Nnodes;
     REFINE[i][2] = Error;

   } 

//====================================================================================================================================
// EOF
//====================================================================================================================================

