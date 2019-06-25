#ifndef _confMap_h
#define _confMap_h

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>


// Vector of real numbers
typedef double * Rarray;


// Vector of complex numbers
typedef double complex * Carray;


// Vector of integer numbers
typedef int * Iarray;





Rarray rarrDef(int n)
{
    double * ptr;

    ptr = (double * ) malloc( n * sizeof(double) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for double\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Carray carrDef(int n)
{
    double complex * ptr;

    ptr = (double complex * ) malloc( n * sizeof(double complex) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Iarray iarrDef(int n)
{
    int * ptr;

    ptr = (int * ) malloc( n * sizeof(int) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integer\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



int fac(int n)
{
    int
        i,
        nfac;

    nfac = 1;

    for (i = 1; i < n; i++) nfac = nfac * (i + 1);

    return nfac;
}



int NC(int N, int M)
{

/** Number of Configurations(NC) of N particles in M states **/

    int
        i,
        n;

    n = 1;

    if  (M > N)
    {
        for (i = N + M - 1; i > M - 1; i --) n = n * i;
        return n / fac(N);
    }

    for (i = N + M - 1; i > N; i --) n = n * i;
    return n / fac(M - 1);
}



Iarray MountNCmat(int N, int M)
{

/** Matrix of all possible outcomes form NC function with
  * NCmat[i + N*j] = NC(i,j).
**/

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
  * the corresponding Fock vector with v[j] being the
  * occupation number on state j.
**/

    int
        i,
        m,
        x;

    m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    // Put the particles in orbitals while has combinations to spend
    while ( k > 0 )
    {
        while ( k - NC(N,m) < 0 ) m = m - 1;

        x = k - NC(N,m);
        while ( x >= 0 )
        {
            v[m] = v[m] + 1; // One more particle in orbital m
            N = N - 1;       // Less one particle to setup
            k = x;
            x = x - NC(N,m);
        }
    }

    // with k zero put the rest in the first state
    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}



int FockToIndex(int N, int M, Iarray NCmat, Iarray v)
{

/** Convert an occupation vector v to a integer number
  * from 0 to the NC(N,M) - 1.
**/

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



Iarray MountFocks(int N, int M)
{

/** All possible occupation vectors in a matrix **/

    int
        i,
        k,
        nc;

    Iarray
        ItoFock,
        occ;

    nc = NC(N,M);
    occ = iarrDef(M);
    ItoFock = iarrDef(nc * M);

    for (k = 0; k < nc; k++)
    {
        IndexToFock(k, N, M, occ);
        for (i = 0; i < M; i++)
        {
            ItoFock[k + nc*i] = occ[i];
        }
    }

    free(occ);
    return ItoFock;
}



#endif
