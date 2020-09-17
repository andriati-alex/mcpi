#ifndef _FermiFockSpace_h
#define _FermiFockSpace_h

/****   CONFIGURATIONAL SPACE FOR FERMIONS WITHOUT RESTRICTIONS

 Basic routines to setup the multiconfiguration space. Functions to setup a
 *hashing table*, convert configuration to index and index to configuration

 A relevant paper for the theoretical background is

    "A perfect Hashing function for exact diagonalization of
    many-body systems of identical particles", Shoudan Liang
    https://doi.org/10.1016/0010-4655(95)00108-R
 
****/

#include "StandardAuxiliarLib.h"



int NCF(int N, int M)
{

/** Number of Configurations(NC) of 'N' (Fe)rmions in 'M' single
    particle state.  A derivation can be consulted in the paper:
 
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



void FindexToConfig(int k, int N, int M, Farray v)
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
        while (k - NCF(N,m) < 0) m = m - 1;
        k = k - NCF(N,m); // subtract cost
        v[m] = 1;        // One more particle in orbital m
        N = N - 1;       // Less one particle to setup
    }

    // with k zero, allocate the rest of particles sequentially
    for (i = N-1; i >= 0; i--) v[i] = 1;
}



int FconfigToIndex(int N, int M, Farray v)
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
            k = k + NCF(N,i); // number of combinations needed
            N = N - 1;       // decrease total # of particles
        }
    }
    return k;
}



Farray FsetupConfigHT(int N, int M)
{

/** A hashing table for the Fock states ordering. It stores index of
    configurations as its row numbers and the occupation along cols.
    The structure is allocated as an 1D array and   ItoFock[j + k*M]
    gives the occupation number of configuration k in the orbital  j **/

    int
        k,
        nc;

    Farray
        ItoFock;

    nc = NCF(N,M);

    if (nc > INT_MAX / M)
    {
        printf("\n\nMEMORY ERROR : Because of the size of the");
        printf(" hashing table it can't be indexed by 32-bit integers\n\n");
        exit(EXIT_FAILURE);
    }

    ItoFock = farrDef(nc*M);

    for (k = 0; k < nc; k++)
    {
        FindexToConfig(k,N,M,&ItoFock[M*k]);
    }

    return ItoFock;
}



#endif
