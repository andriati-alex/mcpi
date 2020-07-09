#ifndef _BoseFockSpace_h
#define _BoseFockSpace_h

/****   AUTHOR INFORMATION
 
NAME : Alex Valerio Andriati
AFFILIATION : University of Sao Paulo - Brazil

Last update : July/02/2019

    * COMMENTS

Basic routines to setup a multiconfiguration problem for bosons.  The
ordering of the config. follows the cost function based on the number
of possibilities to fit indistinguishable balls in different boxes.

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
        printf("\n\n\nMEMORY ERROR : malloc fail for double\n\n");
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
        printf("\n\n\nMEMORY ERROR : malloc fail for complex\n\n");
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
        printf("\n\n\nMEMORY ERROR : malloc fail for integer.");
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
        printf("\n\n\nMEMORY ERROR : malloc fail for (complex *)\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = carrDef(n);

    return ptr;
}



int NC(int N, int M)
{

/** Number of Configurations(NC) of N particles in M states **/

    long
        j,
        i,
        n;

    n = 1;
    j = 2;

    if  (M > N)
    {
        for (i = N + M - 1; i > M - 1; i --)
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

        for (i = j; i < N + 1; i++) n = n / i;

        return ((int) n);
    }

    for (i = N + M - 1; i > N; i --)
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
        if (n % j == 0 && j <= M - 1)
        {
            n = n / j;
            j = j + 1;
        }
    }

    for (i = j; i < M; i++) n = n / i;

    return ((int) n);
}



Iarray setupNCmat(int N, int M)
{

/** Matrix of all possible outcomes from NC function  with
    NCmat[i + N*j] = NC(i,j), where i <= N and j <= M, the
    number of particles and states respectively.

    This is an auxiliar structure to avoid calls of  NC
    function many times when converting Fock states  to
    indexes **/

    int
        i,
        j;

    Iarray
        NCmat;

    NCmat = iarrDef((M + 1) * (N + 1));

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i + (N+1)*0] = 0;
        NCmat[i + (N+1)*1] = 1;
        for (j = 2; j < M + 1; j++) NCmat[i + (N+1)*j] = NC(i,j);
    }

    return NCmat;
}



void IndexToFock(int k, int N, int M, Iarray v)
{

/** Given an integer index 0 < k < NC(N,M) setup on v
    the corresponding Fock vector with v[j] being the
    occupation number on state j.

    This routine corresponds to an implementation of Algorithm 2
    of the article **/

    int
        i,
        m;

    m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    while ( k > 0 )
    {
        // Check if can 'pay' the cost to put  the  particle
        // in current state. If not, try a 'cheaper' one
        while ( k - NC(N,m) < 0 ) m = m - 1;

        k = k - NC(N,m); // subtract cost
        v[m] = v[m] + 1; // One more particle in orbital m
        N = N - 1;       // Less one particle to setup
    }

    // with k zero put the rest in the first state
    if ( N > 0 ) v[0] = v[0] + N;
}



int FockToIndex(int N, int M, Iarray NCmat, Iarray v)
{

/** Convert an occupation vector v to a integer number from
    0 to the NC(N,M) - 1. It uses the  NCmat  structure  to
    avoid calls of NC function, see setupNCmat function above

    This routines is an implementation of algorithm 1 of the article **/

    int
        i,
        k,
        n,
        col;

    k = 0;
    col = N + 1;

    // Empty one by one orbital starting from the last one
    for (i = M - 1; i > 0; i--)
    {
        n = v[i]; // Number of particles in a given orbital

        while (n > 0)
        {
            k = k + NCmat[N + col*i]; // number of combinations needed
            N = N - 1; // decrease total # of particles
            n = n - 1; // decrease # of particles in current state
        }
    }

    return k;
}



Iarray setupFocks(int N, int M)
{

/** A hashing table for the config. ordering.  It stores for each index
    of configurations the occupation numbers of the corresponding  Fock
    state, thus, replacing the usage of IndexToFock routine by a memory
    access.

    ItoFock[j + k*M]  gives the occupation number of configuration k in
    the orbital j **/

    int
        k,
        nc;

    Iarray
        ItoFock;

    nc = NC(N,M);

    if (nc > INT_MAX / M)
    {
        printf("\n\n\nMEMORY ERROR : Because of the size of the ");
        printf("hashing table it cannot be indexed by 32-bit integers\n\n");
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
