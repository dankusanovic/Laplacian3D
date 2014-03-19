//====================================================================================================================================
// HEADER FILE : "MUMPSSolver"
//====================================================================================================================================
// Syntax      : MUMPSSolver(int, int, int*, int*, double*, double*, unsigned int, unsigned int)
//------------------------------------------------------------------------------------------------------------------------------------
// Purpose     : Computes the vector solution for K*x = F.
//====================================================================================================================================
// Written by Danilo S. Kusanovic, December 2013
// Last revised by D.S Kusanovic.
//====================================================================================================================================

#ifndef MUMPSSOLVER_HH
#define MUMPSSOLVER_HH

#define JOB_INIT 			-1
#define JOB_END 			-2
#define USE_COMM_WORLD 		-987654
#define UNSYMMETRIC			 0
#define SYMMETRIC_POSITIVE_DIFINITE	 1
#define GENERAL_SYMMETRIC		 2
#define ICNTL(I)                       icntl[(I)-1]

void MUMPSSolver(int, int, int* &, int* &, double* &, double* &, unsigned int);

#endif

//====================================================================================================================================
// EOF
//====================================================================================================================================

