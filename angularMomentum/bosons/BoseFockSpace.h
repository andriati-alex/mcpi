#ifndef _BoseFockSpace_h
#define _BoseFockSpace_h

/****   AUTHOR INFORMATION
 
NAME : Alex Valerio Andriati
AFFILIATION : University of Sao Paulo - Brazil

Last update : July/02/2019

    COMMENTS

Basic routines to setup a multiconfiguration problem for bosons.  The
ordering of the config. follows the cost function based on the number
of possibilities to fit indistinguishable balls in different boxes.

****/

#include "StandardAuxiliarLib.h"



int NC(int N, int M)
{

/** Number of Configurations(NC) of N particles in M states. It is
    equivalent to the problem of the number of ways to fit N balls
    in M different boxes. **/

    long
        j,
        i,
        n;

    n = 1;
    j = 2;

    // No place to put the particles
    if (M == 0)
    {
        printf("\n\nERROR : 0 single particle states given !\n\n");
        exit(EXIT_FAILURE);
    }

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



void indexToConfig(int k, int N, int M, Iarray v)
{

/** Given an integer index  0 < k < NC(N,M)  setup on v
    the corresponding configuration with v[j] being the
    occupation number on single particle state j.   **/

    int
        i,
        m;

    m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    while ( k > 0 )
    {
        // Check if can 'pay' the cost to put the particle
        // in the current state. If not, try a cheaper one
        while ( k - NC(N,m) < 0 ) m = m - 1;

        k = k - NC(N,m); // subtract cost
        v[m] = v[m] + 1; // One more particle in orbital m
        N = N - 1;       // Less one particle to setup
    }

    // with k zero put the rest in the first state
    if ( N > 0 ) v[0] = v[0] + N;
}



int configToIndex(int N, int M, Iarray v)
{

/** Convert an occupation vector v to a integer number between
    0 and NC(N,M) - 1. It uses the NCmat structure to avoid
    calls of NC function, see setupNCmat function above    **/

    int
        i,
        k,
        n;

    k = 0;

    // Empty one by one orbital starting from the last one
    for (i = M - 1; i > 0; i--)
    {
        n = v[i]; // Number of particles in a given orbital

        while (n > 0)
        {
            k = k + NC(N,i); // number of combinations needed
            N = N - 1; // decrease total # of particles
            n = n - 1; // decrease # of particles in current state
        }
    }

    return k;
}



Iarray setupConfigHT(int N, int M)
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
        printf("\n\nMEMORY ERROR : Because of the size of the ");
        printf("hashing table it cannot be indexed by 32-bit integers\n\n");
        exit(EXIT_FAILURE);
    }

    ItoFock = iarrDef(nc*M);

    for (k = 0; k < nc; k++)
    {
        indexToConfig(k,N,M,&ItoFock[M*k]);
    }

    return ItoFock;
}



#endif
