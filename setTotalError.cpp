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
double get_ux(double x, double y, double z) {

    double ux = -y*z*(2*x - 1)*(y - 1)*(z - 1);

    return ux;
}

double get_uy(double x, double y, double z) {

    double uy = -x*z*(2*y - 1)*(x - 1)*(z - 1);

    return uy;
}

double get_uz(double x, double y, double z) {

    double uz = -x*y*(2*z - 1)*(x - 1)*(y - 1);

    return uz;
}

//------------------------------------------------------------------------------------------------------------------------------------
// ANALYTICAL RHS FUNCTION:
//------------------------------------------------------------------------------------------------------------------------------------
double getForce(double x, double y, double z);

//------------------------------------------------------------------------------------------------------------------------------------
// COMPUTES THE ERROR (L2-NORM):
//------------------------------------------------------------------------------------------------------------------------------------
double setTotalError(double** &REFINE,double** Coordinates, int** Elements, double** GaussPoints, double* Solution, int* Dofs,
                     int Nnodes, int Nelems, int Ngauss, int i, double* &etaK, int** vecK) {
    double etamax=0;
    double Error = 0.0;

    //Total Vector Solution:
    int  DOFi;
    double *uh;
    uh   = new double[Nnodes];

    for(int j = 0; j < Nnodes; j++) {
        DOFi = Dofs[j];

        if(DOFi != -1) {
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

    for(int k = 0; k < Nelems; k++)
    {
        double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;

        x0 = Coordinates[Elements[k][0]][0];
        y0 = Coordinates[Elements[k][0]][1];
        z0 = Coordinates[Elements[k][0]][2];
        x1 = Coordinates[Elements[k][1]][0];
        y1 = Coordinates[Elements[k][1]][1];
        z1 = Coordinates[Elements[k][1]][2];
        x2 = Coordinates[Elements[k][2]][0];
        y2 = Coordinates[Elements[k][2]][1];
        z2 = Coordinates[Elements[k][2]][2];
        x3 = Coordinates[Elements[k][3]][0];
        y3 = Coordinates[Elements[k][3]][1];
        z3 = Coordinates[Elements[k][3]][2];

        double largosL[6];

        largosL[0]=sqrt(pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2));
        largosL[1]=sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));
        largosL[2]=sqrt(pow(x2-x0,2)+pow(y2-y0,2)+pow(z2-z0,2));
        largosL[3]=sqrt(pow(x3-x0,2)+pow(y3-y0,2)+pow(z3-z0,2));
        largosL[4]=sqrt(pow(x3-x1,2)+pow(y3-y1,2)+pow(z3-z1,2));
        largosL[5]=sqrt(pow(x3-x2,2)+pow(y3-y2,2)+pow(z3-z2,2));

        double h_k=0;

        for (int jj=0; jj<6; jj++)
        {
            if (h_k<largosL[jj])
            {
                h_k=largosL[jj];
            }
        }

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

        nx[0] = (uy*vz - vy*uz)/2;
        ny[0] = (uz*vx - vz*ux)/2;
        nz[0] = (ux*vy - vx*uy)/2;

        ux = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][0]][0];
        uy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][0]][1];
        uz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][0]][2];

        vx = Coordinates[Elements[k][3]][0] - Coordinates[Elements[k][0]][0];
        vy = Coordinates[Elements[k][3]][1] - Coordinates[Elements[k][0]][1];
        vz = Coordinates[Elements[k][3]][2] - Coordinates[Elements[k][0]][2];

        nx[1] = (uy*vz - vy*uz)/2;
        ny[1] = (uz*vx - vz*ux)/2;
        nz[1] = (ux*vy - vx*uy)/2;

        ux = Coordinates[Elements[k][3]][0] - Coordinates[Elements[k][0]][0];
        uy = Coordinates[Elements[k][3]][1] - Coordinates[Elements[k][0]][1];
        uz = Coordinates[Elements[k][3]][2] - Coordinates[Elements[k][0]][2];

        vx = Coordinates[Elements[k][1]][0] - Coordinates[Elements[k][0]][0];
        vy = Coordinates[Elements[k][1]][1] - Coordinates[Elements[k][0]][1];
        vz = Coordinates[Elements[k][1]][2] - Coordinates[Elements[k][0]][2];

        nx[2] = (uy*vz - vy*uz)/2;
        ny[2] = (uz*vx - vz*ux)/2;
        nz[2] = (ux*vy - vx*uy)/2;

        ux = Coordinates[Elements[k][1]][0] - Coordinates[Elements[k][0]][0];
        uy = Coordinates[Elements[k][1]][1] - Coordinates[Elements[k][0]][1];
        uz = Coordinates[Elements[k][1]][2] - Coordinates[Elements[k][0]][2];

        vx = Coordinates[Elements[k][2]][0] - Coordinates[Elements[k][0]][0];
        vy = Coordinates[Elements[k][2]][1] - Coordinates[Elements[k][0]][1];
        vz = Coordinates[Elements[k][2]][2] - Coordinates[Elements[k][0]][2];

        nx[3] = (uy*vz - vy*uz)/2;
        ny[3] = (uz*vx - vz*ux)/2;
        nz[3] = (ux*vy - vx*uy)/2;

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

        //Numerical Flux Values:
        uhx = (-nx[0]*uh[Elements[k][0]] - nx[1]*uh[Elements[k][1]] - nx[2]*uh[Elements[k][2]] - nx[3]*uh[Elements[k][3]])/(3*Volume);
        uhy = (-ny[0]*uh[Elements[k][0]] - ny[1]*uh[Elements[k][1]] - ny[2]*uh[Elements[k][2]] - ny[3]*uh[Elements[k][3]])/(3*Volume);
        uhz = (-nz[0]*uh[Elements[k][0]] - nz[1]*uh[Elements[k][1]] - nz[2]*uh[Elements[k][2]] - nz[3]*uh[Elements[k][3]])/(3*Volume);


        //------------------------------------------------------------------------------------------------------------------------------
        // COMPUTES THE ERROR:
        //------------------------------------------------------------------------------------------------------------------------------

        double jump=0;
        //double residual=0;
        for (int tt=0; tt<4; tt++)
        {

            if (vecK[k][tt]>=0)
            {
                double Volumev, nxv[4], nyv[4], nzv[4];
                double uxv, uyv, uzv, vxv, vyv, vzv, wxv, wyv, wzv;
                double uhxv, uhyv, uhzv;//, xiv, etav, zetav, errorv;

                uxv = Coordinates[Elements[vecK[k][tt]][3]][0] - Coordinates[Elements[vecK[k][tt]][1]][0];
                uyv = Coordinates[Elements[vecK[k][tt]][3]][1] - Coordinates[Elements[vecK[k][tt]][1]][1];
                uzv = Coordinates[Elements[vecK[k][tt]][3]][2] - Coordinates[Elements[vecK[k][tt]][1]][2];

                vxv = Coordinates[Elements[vecK[k][tt]][2]][0] - Coordinates[Elements[vecK[k][tt]][1]][0];
                vyv = Coordinates[Elements[vecK[k][tt]][2]][1] - Coordinates[Elements[vecK[k][tt]][1]][1];
                vzv = Coordinates[Elements[vecK[k][tt]][2]][2] - Coordinates[Elements[vecK[k][tt]][1]][2];

                nxv[0] = (uyv*vzv - vyv*uzv)/2;
                nyv[0] = (uzv*vxv - vzv*uxv)/2;
                nzv[0] = (uxv*vyv - vxv*uyv)/2;

                uxv = Coordinates[Elements[vecK[k][tt]][2]][0] - Coordinates[Elements[vecK[k][tt]][0]][0];
                uyv = Coordinates[Elements[vecK[k][tt]][2]][1] - Coordinates[Elements[vecK[k][tt]][0]][1];
                uzv = Coordinates[Elements[vecK[k][tt]][2]][2] - Coordinates[Elements[vecK[k][tt]][0]][2];

                vxv = Coordinates[Elements[vecK[k][tt]][3]][0] - Coordinates[Elements[vecK[k][tt]][0]][0];
                vyv = Coordinates[Elements[vecK[k][tt]][3]][1] - Coordinates[Elements[vecK[k][tt]][0]][1];
                vzv = Coordinates[Elements[vecK[k][tt]][3]][2] - Coordinates[Elements[vecK[k][tt]][0]][2];

                nxv[1] = (uyv*vzv - vyv*uzv)/2;
                nyv[1] = (uzv*vxv - vzv*uxv)/2;
                nzv[1] = (uxv*vyv - vxv*uyv)/2;

                uxv = Coordinates[Elements[vecK[k][tt]][3]][0] - Coordinates[Elements[vecK[k][tt]][0]][0];
                uyv = Coordinates[Elements[vecK[k][tt]][3]][1] - Coordinates[Elements[vecK[k][tt]][0]][1];
                uzv = Coordinates[Elements[vecK[k][tt]][3]][2] - Coordinates[Elements[vecK[k][tt]][0]][2];

                vxv = Coordinates[Elements[vecK[k][tt]][1]][0] - Coordinates[Elements[vecK[k][tt]][0]][0];
                vyv = Coordinates[Elements[vecK[k][tt]][1]][1] - Coordinates[Elements[vecK[k][tt]][0]][1];
                vzv = Coordinates[Elements[vecK[k][tt]][1]][2] - Coordinates[Elements[vecK[k][tt]][0]][2];

                nxv[2] = (uyv*vzv - vyv*uzv)/2;
                nyv[2] = (uzv*vxv - vzv*uxv)/2;
                nzv[2] = (uxv*vyv - vxv*uyv)/2;

                uxv = Coordinates[Elements[vecK[k][tt]][1]][0] - Coordinates[Elements[vecK[k][tt]][0]][0];
                uyv = Coordinates[Elements[vecK[k][tt]][1]][1] - Coordinates[Elements[vecK[k][tt]][0]][1];
                uzv = Coordinates[Elements[vecK[k][tt]][1]][2] - Coordinates[Elements[vecK[k][tt]][0]][2];

                vxv = Coordinates[Elements[vecK[k][tt]][2]][0] - Coordinates[Elements[vecK[k][tt]][0]][0];
                vyv = Coordinates[Elements[vecK[k][tt]][2]][1] - Coordinates[Elements[vecK[k][tt]][0]][1];
                vzv = Coordinates[Elements[vecK[k][tt]][2]][2] - Coordinates[Elements[vecK[k][tt]][0]][2];

                nxv[3] = (uyv*vzv - vyv*uzv)/2;
                nyv[3] = (uzv*vxv - vzv*uxv)/2;
                nzv[3] = (uxv*vyv - vxv*uyv)/2;

                //------------------------------------------------------------------------------------------------------------------------------
                // LOCAL ELEMENT VOLUME COMPUTATION:
                //------------------------------------------------------------------------------------------------------------------------------
                uxv = Coordinates[Elements[vecK[k][tt]][1]][0] - Coordinates[Elements[vecK[k][tt]][3]][0];
                uyv = Coordinates[Elements[vecK[k][tt]][1]][1] - Coordinates[Elements[vecK[k][tt]][3]][1];
                uzv = Coordinates[Elements[vecK[k][tt]][1]][2] - Coordinates[Elements[vecK[k][tt]][3]][2];

                vxv = Coordinates[Elements[vecK[k][tt]][2]][0] - Coordinates[Elements[vecK[k][tt]][3]][0];
                vyv = Coordinates[Elements[vecK[k][tt]][2]][1] - Coordinates[Elements[vecK[k][tt]][3]][1];
                vzv = Coordinates[Elements[vecK[k][tt]][2]][2] - Coordinates[Elements[vecK[k][tt]][3]][2];

                wxv = Coordinates[Elements[vecK[k][tt]][0]][0] - Coordinates[Elements[vecK[k][tt]][3]][0];
                wyv = Coordinates[Elements[vecK[k][tt]][0]][1] - Coordinates[Elements[vecK[k][tt]][3]][1];
                wzv = Coordinates[Elements[vecK[k][tt]][0]][2] - Coordinates[Elements[vecK[k][tt]][3]][2];

                Volumev = (wxv*(uyv*vzv - vyv*uzv) + wyv*(uzv*vxv - vzv*uxv) + wzv*(uxv*vyv - vxv*uyv))/6;

                //Numerical Flux Values:
                uhxv = (-nxv[0]*uh[Elements[vecK[k][tt]][0]] - nxv[1]*uh[Elements[vecK[k][tt]][1]] - nxv[2]*uh[Elements[vecK[k][tt]][2]] - nxv[3]*uh[Elements[vecK[k][tt]][3]])/(3*Volumev);
                uhyv = (-nyv[0]*uh[Elements[vecK[k][tt]][0]] - nyv[1]*uh[Elements[vecK[k][tt]][1]] - nyv[2]*uh[Elements[vecK[k][tt]][2]] - nyv[3]*uh[Elements[vecK[k][tt]][3]])/(3*Volumev);
                uhzv = (-nzv[0]*uh[Elements[vecK[k][tt]][0]] - nzv[1]*uh[Elements[vecK[k][tt]][1]] - nzv[2]*uh[Elements[vecK[k][tt]][2]] - nzv[3]*uh[Elements[vecK[k][tt]][3]])/(3*Volumev);

                double areav=sqrt( nx[tt]*nx[tt] + ny[tt]*ny[tt] +nz[tt]*nz[tt] );

                jump+= pow( (uhx-uhxv)*nx[tt] + (uhy-uhyv)*ny[tt] + (uhz-uhzv)*nz[tt] , 2 )*(1/areav)*(h_k);
            }
        }

        //------------------------------------------------------------------------------------------------------------------------------
        // Performs Gaussian Integration over the element:
        //------------------------------------------------------------------------------------------------------------------------------
        double Value=0;


        for(int j = 0; j < Ngauss; j++)
        {

            //Transformed triangular Coordinates:
            xi   = Coordinates[Elements[k][0]][0]*GaussPoints[j][1] + Coordinates[Elements[k][1]][0]*GaussPoints[j][2] +
                   Coordinates[Elements[k][2]][0]*GaussPoints[j][3] + Coordinates[Elements[k][3]][0]*GaussPoints[j][4];
            eta  = Coordinates[Elements[k][0]][1]*GaussPoints[j][1] + Coordinates[Elements[k][1]][1]*GaussPoints[j][2] +
                   Coordinates[Elements[k][2]][1]*GaussPoints[j][3] + Coordinates[Elements[k][3]][1]*GaussPoints[j][4];
            zeta = Coordinates[Elements[k][0]][2]*GaussPoints[j][1] + Coordinates[Elements[k][1]][2]*GaussPoints[j][2] +
                   Coordinates[Elements[k][2]][2]*GaussPoints[j][3] + Coordinates[Elements[k][3]][2]*GaussPoints[j][4];

            Value+=GaussPoints[j][0]*pow(getForce(xi,eta,zeta),2);

            //Error function evaluated at each triangular coordinates:
            error = error + GaussPoints[j][0]*(pow(get_ux(xi,eta,zeta) - uhx,2) + pow(get_uy(xi,eta,zeta) - uhy,2) + pow(get_uz(xi,eta,zeta) - uhz,2));
        }


        Value*=Volume*h_k*h_k;

        etaK[k]=Value+jump;

        Error  = Error + Volume*error;
        if (etaK[k]>etamax)
            etamax=etaK[k];
    }
    /*
        double contador=0;

        for (int ii=0; ii<Nelems; ii++)
        {
            contador+=etaK[ii];
        }
        //etamax=sqrt(contador);
        etamax=contador;*/

    delete[] uh;
    Error = sqrt(fabs(Error));

    //L2 Norm Information:
    REFINE[i][0] = (double) i + 1;
    REFINE[i][1] = (double) Nnodes;
    REFINE[i][2] = Error;
    REFINE[i][3] = etamax;
    std::cout<< Nnodes << "  " << Error << "  " << etamax<<std::endl;
    return etamax;
}




//====================================================================================================================================
// EOF
//====================================================================================================================================

