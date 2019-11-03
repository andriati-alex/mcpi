#ifndef _configurationsMap_cuh
#define _configurationsMap_cuh

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <stdlib.h>
#include <stdio.h>

/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : November/02/2019

----------------------------------------------------------------------------

 COMMENTS

 Almost all routines are present in the ordinary CPU code. These routines
 have much less description since it is required the familiarity with the
 algorithms.

 The number of threads per block is chosen accordingly to Nvidia CUDA
 developer guide.

****/

#define THREADS_PER_BLOCK 256


typedef double * Rarray;

typedef cuDoubleComplex * Carray;

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



void cuda_rarrDef(int n, Rarray * ptr_address)
{

    // Error code to check return values for CUDA calls
    cudaError_t
        err;

    err = cudaSuccess;

    err = cudaMalloc( (void **) ptr_address, n * sizeof(double));

    if (err != cudaSuccess)
    {
        printf("\n\nFailed to allocate device double vector - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}



Carray carrDef(int n)
{
    cuDoubleComplex * ptr;

    ptr = (cuDoubleComplex * ) malloc( n * sizeof(cuDoubleComplex) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



void cuda_carrDef(int n, Carray * ptr_address)
{

    // Error code to check return values for CUDA calls
    cudaError_t
        err;

    err = cudaSuccess;

    err = cudaMalloc( (void **) ptr_address, n * sizeof(cuDoubleComplex));

    if (err != cudaSuccess)
    {
        printf("\n\nFailed to allocate device complex vector - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}



Iarray iarrDef(int n)
{
    int * ptr;

    ptr = (int * ) malloc( n * sizeof(int) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integer.");
        printf(" Size requested : %ld\n\n", n * sizeof(int));
        exit(EXIT_FAILURE);
    }

    return ptr;
}



void cuda_iarrDef(int n, Iarray * ptr_address)
{

    // Error code to check return values for CUDA calls
    cudaError_t
        err;

    err = cudaSuccess;

    err = cudaMalloc( (void **) ptr_address, n * sizeof(int));

    if (err != cudaSuccess)
    {
        printf("\n\nFailed to allocate device integer vector - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

}



long fac(int n)
{
    long
        i;

    long
        nfac;

    nfac = 1;

    for (i = 1; i < n; i++) nfac = nfac * (i + 1);

    return nfac;
}



int NC(int N, int M)
{

    long
        i,
        n;

    n = 1;

    if  (M > N)
    {
        for (i = N + M - 1; i > M - 1; i --) n = n * i;
        return (int) (n / fac(N));
    }

    for (i = N + M - 1; i > N; i --) n = n * i;
    return (int) (n / fac(M - 1));
}



Iarray setupNCmat(int N, int M)
{

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

    int
        k,
        nc;

    Iarray
        ItoFock;

    nc = NC(N,M);
    ItoFock = iarrDef(nc * M);

    for (k = 0; k < nc; k++)
    {
        IndexToFock(k,N,M,&ItoFock[M*k]);
    }

    return ItoFock;
}



Iarray OneOneMap(int N, int M, Iarray NCmat, Iarray IF)
{

    int i,
        q,
        k,
        l,
        nc;

    Iarray
        v,
        Map;

    nc = NC(N,M);
    
    v = iarrDef(M);

    Map = iarrDef(M * M * nc);

    for (i = 0; i < nc * M * M; i++) Map[i] = -1;

    for (i = 0; i < nc; i++)
    {
        // Copy the occupation vector from C[i] coeff.

        for (q = 0; q < M; q++) v[q] = IF[q + M*i];

        for (k = 0; k < M; k++)
        {
            // Take one particle from k state
            if (v[k] < 1) continue;

            for (l = 0; l < M; l++)
            {
                // Put one particle in l state
                v[k] -= 1;
                v[l] += 1;
                Map[i + k * nc + l * M * nc] = FockToIndex(N,M,NCmat,v);
                v[k] += 1;
                v[l] -= 1;
            }
        }
    }

    free(v);

    return Map;
}



Iarray allocTwoTwoMap(int nc, int M, Iarray IF)
{
    int
        i,
        k,
        s,
        chunks;

    Iarray
        Map;

    chunks = 0;

    for (i = 0; i < nc; i++)
    {

        for (k = 0; k < M; k++)
        {

            if (IF[k + i * M] < 1) continue;

            for (s = k + 1; s < M; s++)
            {

                if (IF[s + i * M] < 1) continue;

                chunks++;
            }
        }
    }

    Map = iarrDef(chunks * M * M);

    for (i = 0; i < M * M * chunks; i++) Map[i] = -1;

    return Map;
}



Iarray TwoTwoMap(int N, int M, Iarray NCmat, Iarray IF, Iarray strideC)
{

    int
        i,
        k,
        s,
        l,
        q,
        nc,
        chunksO,
        chunksC,
        strideO,
        mapIndex;

    Iarray
        occ,
        Map;

    nc = NC(N,M);
    occ = iarrDef(M);

    Map = allocTwoTwoMap(nc,M,IF);

    chunksC = 0;

    for (i = 0; i < nc; i++)
    {
        strideC[i] = chunksC * (M * M);

        for (k = 0; k < M; k++) occ[k] = IF[k + M * i];

        chunksO = 0;

        for (k = 0; k < M; k++)
        {

            if (occ[k] < 1) continue;

            for (s = k + 1; s < M; s++)
            {

                if (occ[s] < 1) continue;

                strideO = chunksO * M * M;

                for (l = 0; l < M; l++)
                {
                    // put one particle in l state
                    for (q = 0; q < M; q++)
                    {
                        // Put one particle in q state
                        mapIndex = strideC[i] + strideO + (l + q * M);
                        occ[k] -= 1;
                        occ[s] -= 1;
                        occ[l] += 1;
                        occ[q] += 1;
                        Map[mapIndex] = FockToIndex(N,M,NCmat,occ);
                        occ[k] += 1;
                        occ[s] += 1;
                        occ[l] -= 1;
                        occ[q] -= 1;
                    }
                }

                chunksO++;
                chunksC++;
            }
        }
    }

    free(occ);

    return Map;
}



Iarray allocOneTwoMap(int nc, int M, Iarray IF)
{

    int
        i,
        k,
        chunks;

    Iarray
        Map;

    chunks = 0;

    for (i = 0; i < nc; i++)
    {

        for (k = 0; k < M; k++)
        {

            if (IF[k + i * M] < 2) continue;

            chunks++;

        }
    }

    Map = iarrDef(chunks * M * M);

    for (i = 0; i < M * M * chunks; i++) Map[i] = -1;

    return Map;
}



Iarray OneTwoMap(int N, int M, Iarray NCmat, Iarray IF, Iarray strideC)
{
    
    int
        i,
        q,
        k,
        l,
        nc,
        chunksO,
        chunksC,
        strideO,
        mapIndex;

    Iarray
        occ,
        Map;

    occ = iarrDef(M);
    nc = NC(N,M);

    Map = allocOneTwoMap(nc,M,IF);

    chunksC = 0;

    for (i = 0; i < nc; i++)
    {

        strideC[i] = chunksC * M * M;

        // Copy the occupation vector from C[i] coeff.
        for (k = 0; k < M; k++) occ[k] = IF[k + i * M];

        chunksO = 0;

        for (k = 0; k < M; k++)
        {

            // Must be able to remove two particles
            if (occ[k] < 2) continue;

            strideO = chunksO * M * M;

            for (l = 0; l < M; l++)
            {
                // put one particle in l state
                for (q = 0; q < M; q++)
                {
                    mapIndex = strideC[i] + strideO + (l + q * M);
                    // Put one particle in q state
                    occ[k] -= 2;
                    occ[l] += 1;
                    occ[q] += 1;
                    Map[mapIndex] = FockToIndex(N,M,NCmat,occ);
                    occ[k] += 2;
                    occ[l] -= 1;
                    occ[q] -= 1;
                }
            }
            chunksO++;
            chunksC++;
        }
    }

    free(occ);

    return Map;
}



#endif
