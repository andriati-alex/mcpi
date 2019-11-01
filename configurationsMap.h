#ifndef _configurationsMap_h
#define _configurationsMap_h

/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : 10/31/2019
 
****/

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
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integer.");
        printf(" Size requested : %ld\n\n", n * sizeof(int));
        exit(EXIT_FAILURE);
    }

    return ptr;
}



long fac(int n)
{
    long
        i,
        nfac;

    nfac = 1;

    for (i = 1; i < n; i++) nfac = nfac * (i + 1);

    return nfac;
}



int NC(int N, int M)
{

/** Number of Configurations(NC) of N particles in M states **/

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

/** Matrix of all possible outcomes form NC function  with
  * NCmat[i + N*j] = NC(i,j), where i <= N and j <= M, the
  * number of particles and states respectively.
  *
  * This is an auxiliar structure to avoid calls of  NC
  * function many times when converting Fock states  to
  * indexes **/

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
  *
  * This routine corresponds to an implementation of Algorithm 2
  * of the article **/

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
  * 0 to the NC(N,M) - 1. It uses the  NCmat  structure  to
  * avoid calls of NC function, see setupNCmat function above
  *
  * This routines is an implementation of algorithm 1 of the article **/

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

/** A hashing table for the Fock states ordering. Stores for each index
  * of configurations the occupation numbers of the corresponding  Fock
  * state, thus, replacing the usage of IndexToFock routine by a memory
  * access.
  *
  * ItoFock[j + k*M]  gives the occupation number of configuration k in
  * the orbital j **/

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

/** Given the first configuration index, map it in  a second  one  which
  * the occupation vector differs from the first by a jump of a particle
  * from one state to another.
  *
  * Thus given the first index 'i' of a Fock state :
  *
  * Map[i + k * nc + l * nc * M] = index of another Fock state which
  * have one particle less in k that has been added in l         **/

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

    // The structure consider that for any configuration, there are
    // M^2 possible jumps among the individual particle states.  In
    // spite of the wasted elements from forbidden transitions that
    // are based on states that are empty, this is no problem compared
    // to the routines that maps double jumps.
    Map = iarrDef(M * M * nc);

    // Forbidden transitions are marked with -1
    for (i = 0; i < nc * M * M; i++) Map[i] = -1;

    for (i = 0; i < nc; i++)
    {
        // Copy the occupation vector from C[i] coeff.

        for (q = 0; q < M; q++) v[q] = IF[q + M*i];

        for (k = 0; k < M; k++)
        {
            // check if there is at least a particle to remove
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

/** Structure allocation of mapping between two different Fock states
  * whose the occupation numbers are related by jumps of  2 particles
  * from different individual particle states
  *
  * Given an non-empty orbital k, look for the next non-empty s > k
  * When found such a combination, it is necessary to allocate  M^2
  * new elements corresponding to the particles destiny, those that
  * were removed from states k and s  **/

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

/** Structure to direct map a configuration to another by replacing
  * particle from two different orbitals. To build such a structure
  * in a vector of integers it looks in each configuration how many
  * different possibilities are to remove two particles  from  two
  * different states, and for  each time  it happens there are M^2
  * different places to put those particles. Thus this function also
  * has as output the last argument, vector strideC which  for  each
  * enumerated configuration i store the integer number, a index  of
  * the mapping where those possibilites to remove two particle starts.
  *
  * EXAMPLE : Given a configuration i, find configuration j which has
  * a particle less in states 'k' ans 's' (s > k),  and  two  more on
  * states 'q' and 'l'.
  *
  * SOL : Using the map returned by this structure we start  by  the
  * stride from the configuration i, excluding the mapped index from
  * all previous configurations. Then, walk in chunks of size M^2 until
  * reach the orbitals desired to remove the particles.
  *
  * m = strideC[i];
  *
  * for h = 0 ... k - 1
  * {
  *     for g = h + 1 ... M - 1
  *     {
  *         if occupation on h and g are greater than 1 then
  *         {
  *             m = m + M * M;
  *         }
  *     }
  * }
  *
  * for g = k + 1 ... s - 1
  * {
  *     if occupation on k and g are greater than 1 then
  *     {
  *         m = m + M * M;
  *     }
  * }
  *
  * after adding the orbital strides proportional to the M^2 possibilities
  * where the particles removed can be placed, finish adding
  *
  * m = q + l * M;
  *
  * j = MapTT[m]; **/

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
        // strideC[i] is the index where jumps based on Fock state  i begins
        // chunksC counts how many possibilities were found for the previous
        // Fock states, and for each possibility M^2 is the size of the chunk
        // because of the number of possible destinies for removed particles
        strideC[i] = chunksC * (M * M);

        for (k = 0; k < M; k++) occ[k] = IF[k + M * i];

        chunksO = 0;

        for (k = 0; k < M; k++)
        {

            if (occ[k] < 1) continue;

            for (s = k + 1; s < M; s++)
            {

                if (occ[s] < 1) continue;

                // If it is possible to remove particles from k and s
                // we need to setup a chunk that corresponds  to  all
                // possible destinies for the particles that are going
                // to be removed from k and s orbitals

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

                // New chunk set up. Update how many chunks have been done
                chunksO++; // for current configuration
                chunksC++; // at all
            }
        }
    }

    free(occ);

    // the final size of this mapping is given by the strideC[nc-1], since
    // the last configuration cannot have a double jump of particles  from
    // two different states.

    return Map;
}



Iarray allocOneTwoMap(int nc, int M, Iarray IF)
{

/** Analogously to allocTwoTwoMap create a structure for mapping between
  * two different Fock states, though here the occupation numbers are
  * related by jumps of 2 particles from the same orbital
  *
  * Given an orbital k that has at least 2 particles of a  configuration
  * there are M^2 possible orbitals to put these particle removed from k
  * For each time it happens among all configurations add a chunck of M^2
  * elements to the total size of the mapping array **/

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

/** Configure the mapping array of jump of 2 particle from the same orbital.
  * In contrast to the TwoTwoMap, for each configuration, for each orbital
  * that has more than 2 particles, store the index of configurations with
  * all possible destinies for the 2 removed particles.
  *
  * EXAMPLE : Given a configuration i, find configuration j which has
  * two particle less in state 'k', and place them in states 'q' and 'l'
  *
  * m = strideC[i];
  *
  * for h = 0 ... k - 1
  * {
  *     if occupation on h are greater than 2 then
  *     {
  *         m = m + M * M;
  *     }
  * }
  *
  * j = MapTT[m + q + l*M]; **/

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

        // stride for where the transitions of configuration i starts
        strideC[i] = chunksC * M * M;

        // Copy the occupation vector from C[i] coeff.
        for (k = 0; k < M; k++) occ[k] = IF[k + i * M];

        chunksO = 0;

        for (k = 0; k < M; k++)
        {

            // Must be able to remove two particles
            if (occ[k] < 2) continue;

            // jump a stride of previous transitions in this same conf.
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

    // The size of this mapping array is given by strideC[nc-1]  plus the
    // number of possible double jumps from the last configurations. From
    // the last configuration we can only remove particles from the  last
    // orbital, them the total size is strideC[nc-1] + M^2

    return Map;
}



#endif
