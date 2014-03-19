//====================================================================================================================================
// IMPLEMENTATION FILE: "MUMPSSolver"
//====================================================================================================================================
// Syntax      : MUMPSSolver(Nfree, Nzeros, row, col, Stiffness, Force, sym)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Computes the vector solution for K*x = F.
//------------------------------------------------------------------------------------------------------------------------------------
// Input       : Nfree       : Number of unknowns                                 [1,1]
//               Nzeros      : Number of components in COO format                 [1,1]
//               row         : Row vector components                              [Nzeros,1]
//               col         : Column vector components                           [Nzeros,1]
//               Stiffness   : Stiffness matrix vector values                     [Nzeros,1]
//               Force       : Force vector values                                [Nfree,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Output      : Force       : Solution vector values                             [Nfree,1]
//------------------------------------------------------------------------------------------------------------------------------------
// Folder      :
//------------------------------------------------------------------------------------------------------------------------------------
// Description : <Solves linear system>
// Keywords    : <solver,Mumps,configuration>
// Version     : <>
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

#include <dmumps_c.h>
#include "MUMPSSolver.hh"

void MUMPSSolver(int Nfree, int Nzeros, int* &row, int* &col, double* &Stiffness, double* &Force, unsigned int sym) {

//------------------------------------------------------------------------------------------------------------------------------------
// INITIALIZES MUMPS SOLVER:
//------------------------------------------------------------------------------------------------------------------------------------
    DMUMPS_STRUC_C Model;

    Model.par = 1;
    Model.sym = sym;
    Model.comm_fortran = USE_COMM_WORLD;

    //Initialize a MUMPS instance
    Model.job = JOB_INIT;
    dmumps_c(&Model);

//------------------------------------------------------------------------------------------------------------------------------------
// DEFINES THE PROBLEM ON THE HOST:
//------------------------------------------------------------------------------------------------------------------------------------
    Model.n   = Nfree;
    Model.nz  = Nzeros;
    Model.irn = row;
    Model.jcn = col;
    Model.a   = Stiffness;
    Model.rhs = Force;

    //No outputs Messages: Erros, warnings.
    Model.ICNTL(1)  = 0;
    Model.ICNTL(2)  = 0;
    Model.ICNTL(3)  = 0;
    Model.ICNTL(4)  = 0;

    //Percentage increase in the estimated working space.
    Model.ICNTL(14) = 50;

    //Calling the MUMPS package:
    Model.job = 6;
    dmumps_c(&Model);

    //Terminate a MUMPS instance:
    Model.job = JOB_END;
    dmumps_c(&Model);

}

//====================================================================================================================================
// EOF
//====================================================================================================================================

