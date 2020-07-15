#ifndef _FermiFockSpace_h
#define _FermiFockSpace_h



/****   AUTHOR INFORMATION

 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : July/01/2020

 COMMENTS

 Basic routines to setup the multiconfiguration space. Functions to setup a
 *hashing table*, convert configuration to index and index to configuration

 A relevant paper for the theoretical background is

    "A perfect Hashing function for exact diagonalization of
    many-body systems of identical particles", Shoudan Liang
    https://doi.org/10.1016/0010-4655(95)00108-R
 
****/



#include <limits.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

// Vector of real numbers
typedef double * Rarray;

// Vector of complex numbers
typedef double complex * Carray;

// Matrix of complex numbers
typedef double complex ** Cmatrix;

// Vector of integer numbers
typedef int * Iarray;





Rarray rarrDef(unsigned int n)
{
    double * ptr;

    ptr = (double * ) malloc( n * sizeof(double) );

    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for double. ");
        printf("Size required : %d\n\n",n);
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Carray carrDef(unsigned int n)
{
    double complex * ptr;

    ptr = (double complex * ) malloc( n * sizeof(double complex) );

    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for complex. ");
        printf("Size required : %d\n\n",n);
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Iarray iarrDef(unsigned int n)
{
    int * ptr;

    ptr = (int * ) malloc( n * sizeof(int) );

    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for integer.");
        printf(" Size requested : %ld\n\n", n * sizeof(int));
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Cmatrix cmatDef(unsigned int m, unsigned int n)
{

    int i;

    double complex ** ptr;

    ptr = (double complex ** ) malloc( m * sizeof(double complex *) );

    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for (complex *)\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = carrDef(n);

    return ptr;
}



int NC(int N, int M)
{

/** Number of Configurations(NC) of 'N' fermions in 'M' individual
    particle state. A derivation can be consulted in the paper:
 
    "A perfect Hashing function for exact diagonalization of
    many-body systems of identical particles", Shoudan Liang
    https://doi.org/10.1016/0010-4655(95)00108-R         **/

    long
        j,
        i,
        n;

    // Pauli exclusion principle forbid more particles than IPS
    if (N > M)  return 0;
    if (M == 0) return 0;

    n = 1;
    j = 2;

    if (M > 2*N)
    {
        for (i = M; i > M - N; i--)
        {
            if (n > INT_MAX / i)
            {
                printf("\n\n\nINTEGER SIZE ERROR : overflow occurred");
                printf(" representing the number of configurations as ");
                printf("integers of 32 bits. Config. space without any ");
                printf("restrictions is too large for Npar = %d and ",N);
                printf("Norb = %d\n\n",M);
                exit(EXIT_FAILURE);
            }
            n = n * i;
            if (n % j == 0 && j <= N)
            {
                n = n / j;
                j = j + 1;
            }
        }

        for (i = j; i <= N; i++) n = n / i;

        return ((int) n);
    }

    for (i = M; i > N; i--)
    {
        if (n > INT_MAX / i)
        {
            printf("\n\n\nINTEGER SIZE ERROR : overflow occurred");
            printf(" representing the number of configurations as ");
            printf("integers of 32 bits. Config. space without any ");
            printf("restrictions is too large for Npar = %d and ",N);
            printf("Norb = %d\n\n",M);
            exit(EXIT_FAILURE);
        }
        n = n * i;
        if (n % j == 0 && j <= M - N)
        {
            n = n / j;
            j = j + 1;
        }
    }

    for (i = j; i <= M - N; i++) n = n / i;

    return ((int) n);
}



void IndexToFock(int k, int N, int M, Iarray v)
{

/** Given an integer index 0 < k < NC(N,M) setup on v
    the corresponding Fock vector with v[j] being the
    occupation number on state j, that is 0 or 1  **/

    int
        i,
        m;

    m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    while ( k > 0 )
    {
        // Check if can 'pay' the cost to put the particle
        // in current state.  If not,  try a 'cheaper' one
        while (k - NC(N,m) < 0) m = m - 1;

        k = k - NC(N,m); // subtract cost
        v[m] = 1;        // One more particle in orbital m
        N = N - 1;       // Less one particle to setup
    }

    // with k zero, allocate the rest of particles sequentially
    for (i = N-1; i >= 0; i--) v[i] = 1;
}



int FockToIndex(int N, int M, Iarray v)
{

/** Convert an occupation vector v to a integer number
    from 0 to the NC(N,M) - 1,  index of configuration **/

    int
        i,
        k;

    k = 0;

    // Empty IPS adding the cost to assemble the particles
    for (i = M - 1; i > 0; i--)
    {
        if (v[i] > 0)
        {
            k = k + NC(N,i); // number of combinations needed
            N = N - 1;       // decrease total # of particles
        }
    }

    return k;
}



Iarray setupFocks(int N, int M)
{

/** A hashing table for the Fock states ordering. It stores index of
    configurations as its row numbers and the occupation along cols.
    The structure is allocated as an 1D array and   ItoFock[j + k*M]
    gives the occupation number of configuration k in the orbital  j **/

    int
        k,
        nc;

    Iarray
        ItoFock;

    nc = NC(N,M);

    if (nc > INT_MAX / M)
    {
        printf("\n\nMEMORY ERROR : Because of the size of the");
        printf(" hashing table it can't be indexed by 32-bit integers\n\n");
        exit(EXIT_FAILURE);
    }

    ItoFock = iarrDef(nc*M);

    for (k = 0; k < nc; k++)
    {
        IndexToFock(k,N,M,&ItoFock[M*k]);
    }

    return ItoFock;
}



#endif
