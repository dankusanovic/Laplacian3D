//====================================================================================================================================
// IMPLEMENTATION FILE: "setForceVector"
//====================================================================================================================================
// Syntax      : setStiffnessMatrix(Coordinates, Elements, GaussPoints,Force,Nelems,Ngauss)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Computes the stiffness matrix. 
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Coordinates : List of coordinate values                          [Nnodes,3]
//               Elements    : List of elements values                            [Nelems,3]
//               GaussPoints : List of Gauss Integration values                   [Ngauss,3]
//               Nelems      : Number of elements                                 [1,1]
//               Ngauss      : Number of integration points                       [1,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Force       : Force vector                                       [Nnodes,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      : 
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Stiffness Matrix>
// Keywords    : <Force,RHS>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

  #include <cmath>

//------------------------------------------------------------------------------------------------------------------------------------
// COMPUTES THE STIFFNESS MATRIX: 
//------------------------------------------------------------------------------------------------------------------------------------
   void setStiffnessMatrix(double** &Coordinates, int** &Elements, double** &GaussPoints, int** &Restraints, double* &Force, int Nnodes, 
                           int Nelems, int Ngauss){
        
   //Element Coordinates:
     double xi, eta, zeta;  

   //Element properties: 
     double Volume, Value, Integral[4];

   //Side vector components:
     double ux, uy, uz, vx, vy, vz, wx, wy, wz;

     for(int i = 0; i < Nelems; i++){

      //------------------------------------------------------------------------------------------------------------------------------
      // Computes Tetrahedron volume given four coordinates:
      //------------------------------------------------------------------------------------------------------------------------------
         ux = Coordinates[Elements[i][1]][0] - Coordinates[Elements[i][0]][0];
         uy = Coordinates[Elements[i][1]][1] - Coordinates[Elements[i][0]][1];
         uz = Coordinates[Elements[i][1]][2] - Coordinates[Elements[i][0]][2];

         vx = Coordinates[Elements[i][2]][0] - Coordinates[Elements[i][0]][0];
         vy = Coordinates[Elements[i][2]][1] - Coordinates[Elements[i][0]][1];
         vz = Coordinates[Elements[i][2]][2] - Coordinates[Elements[i][0]][2];

         wx = Coordinates[Elements[i][3]][0] - Coordinates[Elements[i][0]][0];
         wy = Coordinates[Elements[i][3]][1] - Coordinates[Elements[i][0]][1];
         wz = Coordinates[Elements[i][3]][2] - Coordinates[Elements[i][0]][2];

         Volume = fabs(ux*(vy*wz - vz*wy) - uy*(vx*wz - vz*wx) + uz*(vx*wy - vy*wx))/6;

      //------------------------------------------------------------------------------------------------------------------------------
      // Performs Gaussian Integration over the element:
      //------------------------------------------------------------------------------------------------------------------------------
         Integral[0] = 0.0; Integral[1] = 0.0; 
         Integral[2] = 0.0; Integral[3] = 0.0;

         for(int j = 0; j < Ngauss; j++){

          //Transformed triangular Coordinates:
             xi   = Coordinates[Elements[i][0]][0]*GaussPoints[j][1] + Coordinates[Elements[i][1]][0]*GaussPoints[j][2] + 
                    Coordinates[Elements[i][2]][0]*GaussPoints[j][3] + Coordinates[Elements[i][3]][0]*GaussPoints[j][4];
             eta  = Coordinates[Elements[i][0]][1]*GaussPoints[j][1] + Coordinates[Elements[i][1]][1]*GaussPoints[j][2] + 
                    Coordinates[Elements[i][2]][1]*GaussPoints[j][3] + Coordinates[Elements[i][3]][1]*GaussPoints[j][4];
             zeta = Coordinates[Elements[i][0]][2]*GaussPoints[j][1] + Coordinates[Elements[i][1]][2]*GaussPoints[j][2] + 
                    Coordinates[Elements[i][2]][2]*GaussPoints[j][3] + Coordinates[Elements[i][3]][2]*GaussPoints[j][4];

          //Force function evaluated at each triangular coordinates:
            Value = GaussPoints[j][0]*getForce(xi,eta,zeta);

          //Integration evaluated at each triangular Coordinates:
            Integral[0] += Value*(1.0 - GaussPoints[j][2] - GaussPoints[j][3] - GaussPoints[j][4]);
            Integral[1] += Value*(GaussPoints[j][2]);
            Integral[2] += Value*(GaussPoints[j][3]);
            Integral[3] += Value*(GaussPoints[j][4]);
            
         }

      //------------------------------------------------------------------------------------------------------------------------------
      // Force value at each degree of freedom:
      //------------------------------------------------------------------------------------------------------------------------------
         for(int j = 0; j < 4; j++){
             Force[Elements[i][j]] += Volume*Integral[j];
         }

     }

   } 

//====================================================================================================================================
// EOF
//====================================================================================================================================
