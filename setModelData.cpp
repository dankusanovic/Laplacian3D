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
#include <vector>
#define ADJAC1 NodeElemAdjac[Elements[i][(j+1)%4]]
#define ADJAC2 NodeElemAdjac[Elements[i][(j+2)%4]]
#define ADJAC3 NodeElemAdjac[Elements[i][(j+3)%4]]
using std::vector;

void setModelData(double** &Coordinates, int** &Elements, int** &vecK, int** &Boundaries, int** &MeshRefine, double ** &GaussPoints,
                  int* &Restraints, int &Nnodes, int &Nelems, int &Nfaces, int &Ngauss, double* & etaK) {

//------------------------------------------------------------------------------------------------------------------------------------
// READ THE BOUNDARY FACES:
//------------------------------------------------------------------------------------------------------------------------------------
    std::ifstream INFILEbound("Boundaries.txt");

    INFILEbound >> Nfaces;

    int BFaces[Nfaces];

    for(int i = 0; i < Nfaces; i++) {
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

    //Save Values in Coordinates:
    for(int i = 0; i < Nnodes; i++) {
        Coordinates[i] = new double[3];
        INFILEdata >> Coordinates[i][0] >> Coordinates[i][1] >> Coordinates[i][2];
    }

//------------------------------------------------------------------------------------------------------------------------------------
// READ ELEMENTS INFORMATION:
//------------------------------------------------------------------------------------------------------------------------------------
    INFILEdata >> Nelems;

    //Memory Allocation for Elements:
    Elements   = new int* [Nelems];
    MeshRefine = new int* [Nelems];

    for(int i = 0; i < Nelems; i++) {
        Elements[i]   = new int[4];
        MeshRefine[i] = new int[3];
    }

    //Save Values in Elements:
    for(int i = 0; i < Nelems; i++) {
        INFILEdata >> MeshRefine[i][0] >> MeshRefine[i][1] >> MeshRefine[i][2]
                   >>   Elements[i][0] >>   Elements[i][1] >>   Elements[i][2] >> Elements[i][3];

        //C++ Format index:
        Elements[i][0]--;
        Elements[i][1]--;
        Elements[i][2]--;
        Elements[i][3]--;
    }

//------------------------------------------------------------------------------------------------------------------------------------
// Get elements adjacency
//------------------------------------------------------------------------------------------------------------------------------------
    vecK = new int* [Nelems];
    std::vector< std::vector<int> > NodeElemAdjac(Nnodes);
    for (int i=0; i<Nelems; i++)
    {
        /*Initialization of vector of e-e adjacency*/
        vecK[i] = new int[4];
        vecK[i][0]=-1;
        vecK[i][1]=-1;
        vecK[i][2]=-1;
        vecK[i][3]=-1;
        /**/
        NodeElemAdjac[Elements[i][0]].push_back(i);
        NodeElemAdjac[Elements[i][1]].push_back(i);
        NodeElemAdjac[Elements[i][2]].push_back(i);
        NodeElemAdjac[Elements[i][3]].push_back(i);
    }
    for (int i=0; i<Nelems; i++)
    {
        for (int j=0; j<4; j++)
        {
            int k1 = 0;
            int k2 = 0;
            int k3 = 0;

            int size1=ADJAC1.size();
            int size2=ADJAC2.size();
            int size3=ADJAC3.size();


            while (k1 <  size1 && k2 < size2 && k3 < size3) {
                if (ADJAC1[k1]<ADJAC2[k2])
                    ++k1;
                else if (ADJAC2[k2]<ADJAC1[k1])
                    ++k2;
                if (ADJAC1[k1]<ADJAC3[k3])
                    ++k1;
                else if (ADJAC3[k3]<ADJAC1[k1])
                    ++k3;
                if (ADJAC3[k3]==ADJAC2[k2]&&ADJAC3[k3]==ADJAC1[k1]) {
                    if (ADJAC1[k1]!=i) {
                        vecK[i][j]=ADJAC1[k1];
                        break;
                    }
                    else {
                        ++k1;
                        ++k2;
                        ++k3;
                    }
                }

            }

        }
    }
    /*
    for (int i=0;i<Nnodes;i++){
    std::cout<<i;
    for (int j=0;j<NodeElemAdjac[i].size();j++)
        	std::cout<< "\t"<<NodeElemAdjac[i][j];
    std::cout<<std::endl;
    }
    */

//------------------------------------------------------------------------------------------------------------------------------------
// READ BOUNDARY INFORMATION:
//------------------------------------------------------------------------------------------------------------------------------------
    INFILEdata >> Nfaces;

    //Memory Allocation for Boundaries:
    Boundaries = new int* [Nfaces];

    for(int i = 0; i < Nfaces; i++) {
        Boundaries[i] = new int[4];
    }

    //Memory Allocation for Restraints:
    Restraints = new int[Nnodes];

    //Initialization of Restraints Vector:
    for(int i = 0; i < Nnodes; i++) {
        Restraints[i] = 0;
    }

    //Save Values in Boundaries:
    for(int i = 0; i < Nfaces; i++) {
        INFILEdata >> Boundaries[i][0] >> Boundaries[i][1] >> Boundaries[i][2] >> Boundaries[i][3];

        //C++ Format index:
        Boundaries[i][0]--;
        Boundaries[i][1]--;
        Boundaries[i][2]--;
        Boundaries[i][3]--;

        switch(BFaces[Boundaries[i][0]]) {
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

    for(int i = 0; i < Ngauss; i++) {
        GaussPoints[i] = new double[5];
    }

    //Save Values in Coordinates:
    for(int i = 0; i < Ngauss; i++) {
        INFILEgauss >> GaussPoints[i][0] >> GaussPoints[i][1] >> GaussPoints[i][2] >> GaussPoints[i][3] >> GaussPoints[i][4];
    }

    INFILEgauss.close();


    //Memory Allocation for Coordinates:
    etaK = new double [Nelems];
}



//====================================================================================================================================
// EOF
//====================================================================================================================================

