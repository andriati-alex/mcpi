#ifndef _LBoseBoseFockSpace_h
#define _LBoseBoseFockSpace_h

#include "LBoseFockSpace.h"



CompoundSpace AllocCompBasis(int Na, int Nb, int lmaxA, int lmaxB, int L)
{

/** Na : Number of particles of species A
    Nb : Number of particles of species B
    lmaxA : maximum momentum for species A
    lmaxB : maximum momentum for species B
    L : Total (combined) momentum from both species
    --------------------------------------------------------------------
    The compound basis is build assembling an array of subspaces,  which
    are characterized by one of the possible combination of  the  spaces
    with fixed momentum La and Lb for each species, such that the sum of
    both give us L. The subspaces are the tensor product of these  fixed
    momentum spaces for every possible combination of La and Lb.     **/

    int
        i,
        j,
        La,
        Lb,
        size,
        n_sub,
        LboundA;

    long
        MemReq;

    Iarray
        nc_a,
        nc_b;

    struct _CompoundSpace
        * S;

    S = (CompoundSpace) malloc(sizeof(struct _CompoundSpace));
    S->L = L;
    S->Na = Na;
    S->Nb = Nb;
    S->lmaxA = lmaxA;
    S->lmaxB = lmaxB;
    LboundA = Na*lmaxA;

    // array to store the number of config.(size of individual basis)
    // for each species when sweeping through the possible values  of
    // the momentum. For bosons the values for 'La'  are in the range
    // - Na * lmax , ... , 0 , ... , Na * lmax. Once 'La' is fixed
    // we must have Lb = L - La
    nc_a = iarrDef(2*LboundA+1);
    nc_b = iarrDef(2*LboundA+1);

    i = 0;
    size = 0;
    n_sub = 0;
    MemReq = 0;
    // COMPUTE THE SIZE OF THE SPACE  ENSURING
    // THAT IT DO NOT EXCEED SYSTEM CAPABILITY
    for (La = -LboundA; La <= LboundA; La++)
    {
        // Search for the possible combinations of momenta of spaces A and B
        Lb = L-La;
        nc_a[i] = BFixedMom_mcsize(Na,lmaxA,La);
        nc_b[i] = BFixedMom_mcsize(Nb,lmaxB,Lb);
        if (nc_a[i]*nc_b[i] > 0)
        {
            // When the separate spaces with momenta 'La' and 'Lb'
            // coexist, they configure a subspace of the compound
            // system through the tensor product
            n_sub++; // Add counter of subspaces found

            // Assess if it is possible the enumeration of the basis
            if (size > INT_MAX - nc_a[i]*nc_b[i])
            {
                printf("\n\nINDEX ERROR : the structure for the basis ");
                printf("of the compound multiconfig. basis cannot be ");
                printf("indexed by 32-bit integers\n\n");
                printf("Program aborted at function 'AllocCompBasis' ");
                printf("in file 'LBoseBoseFockSpace.h'\n\n");
                exit(EXIT_FAILURE);
            }
            size = size + nc_a[i]*nc_b[i]; // size of compound basis updated

            MemReq += sizeof(struct _TensorProd_ht);
            MemReq += (nc_a[i]*(2*lmaxA+1) + nc_b[i]*(2*lmaxB+1))*sizeof(int);
            if (MemReq > MEMORY_TOL)
            {
                printf("\n\nPROGRAM ABORTED : Estimated memory required ");
                printf("to setup the compound config. basis exceeded the ");
                printf("tolerance %.1lf(GB)\n\n",((double) MEMORY_TOL)/1E9);
                exit(EXIT_FAILURE);
            }
        }
        i++;
    }

    S->sub = (TensorProd_ht) malloc(n_sub*sizeof(struct _TensorProd_ht));
    S->strides = iarrDef(n_sub);
    S->n_sub = n_sub;   // Number of subspaces each one with 'La/Lb' fixed
    S->size = size;     // number of config. for compound basis

    j = 0;
    i = 0;
    size = 0;
    for (La = -LboundA; La <= LboundA; La++)
    {
        Lb = L-La;
        if (nc_a[i]*nc_b[i] > 0)
        {
            // configure the subspace
            S->strides[j] = size;
            S->sub[j].La = La;
            S->sub[j].Lb = Lb;
            S->sub[j].sizeA = nc_a[i];
            S->sub[j].sizeB = nc_b[i];
            S->sub[j].hta = BAssembleHT(Na,lmaxA,La,nc_a[i]);
            S->sub[j].htb = BAssembleHT(Nb,lmaxB,Lb,nc_b[i]);
            size = size + nc_a[i]*nc_b[i];
            j++;
        }
        i++;
    }

    free(nc_a);
    free(nc_b);

    return S;
}



unsigned long EstimateMemory(CompoundSpace S)
{

/** Compute approximately the memory required to set up the ccompound space **/

    unsigned long
        j,
        sub_mem,
        mem_req,
        lmaxA = S->lmaxA,
        lmaxB = S->lmaxB;

    sub_mem = 0;
    for (j = 0; j < S->n_sub; j++)
    {
        sub_mem = sub_mem + S->sub[j].sizeA*(2*lmaxA+1)*sizeof(int);
        sub_mem = sub_mem + S->sub[j].sizeB*(2*lmaxB+1)*sizeof(int);
    }

    mem_req = sub_mem + S->n_sub*(sizeof(struct _TensorProd_ht)+sizeof(int));
    return mem_req;
}



void freeCompSpace(CompoundSpace S)
{
    int
        j,
        i;

    free(S->strides);

    // Free the Hashing table for each subspace
    for (j = 0; j < S->n_sub; j++)
    {
        for (i = 0; i < S->sub[j].sizeA; i++)
        {
            free(S->sub[j].hta[i]);
        }
        free(S->sub[j].hta);
        for (i = 0; i < S->sub[j].sizeB; i++)
        {
            free(S->sub[j].htb[i]);
        }
        free(S->sub[j].htb);
    }

    free(S->sub);
    free(S);
}



void BBgetConfigs(int k, CompoundSpace S, Iarray occA, Iarray occB)
{

/** Set up in 'occA' and 'occB' the occupation numbers corresponding to
    configuration 'k' of the compound space 'S'                     **/

    int
        i,
        j,
        La,
        Lb,
        Nla,
        Nlb,
        indexA,
        indexB,
        subIndex,
        SubSizeA,
        SubSizeB;

    if (k > S->size)
    {
        printf("\n\nERROR : The index %d of the configurations ",k);
        printf("exceed the size of the compound space %d\n\n",S->size);
        exit(EXIT_FAILURE);
    }

    // Total number of IPS for each species
    Nla = 2 * S->lmaxA + 1;
    Nlb = 2 * S->lmaxB + 1;

    // First, find the index of the subspace comparing to strides
    // related to number of states in each subspace
    for (j = 1; j < S->n_sub; j++)
    {
        if (k < S->strides[j]) break;
    }
    subIndex = j-1;             // subspace number
    k = k - S->strides[j-1];    // subtract cost of accessing the subspace

    SubSizeA = S->sub[subIndex].sizeA;
    SubSizeB = S->sub[subIndex].sizeB;
    La = S->sub[subIndex].La;
    Lb = S->sub[subIndex].Lb;

    // Assert k has the right range in the given subspace
    // it must be in the interval [0,sizeA*sizeB)
    if (k >= SubSizeA*SubSizeB)
    {
        printf("\n\nERROR : wrong subspace index found %d - ",k);
        printf("Subspace with La = %d and Lb = %d ",La,Lb);
        printf("has %d elements\n\n",SubSizeA*SubSizeB);
        exit(EXIT_FAILURE);
    }

    indexB = k / SubSizeA;
    indexA = k % SubSizeA;

    // Setup the occupations as output parameters
    for (i = 0; i < Nla; i++) occA[i] = S->sub[subIndex].hta[indexA][i];
    for (i = 0; i < Nlb; i++) occB[i] = S->sub[subIndex].htb[indexB][i];
}



int BBsubIndex(CompoundSpace S, Iarray occA)
{

/** From occupation numbers in 'occA' compute the angular momentum 'La'
    and return the index of the subspace. The subspaces are sorted  as
    increasing values of 'La'                                      **/

    int
        i,
        n,
        La,
        Nla,
        upper,
        lower;

    Nla = 2 * S->lmaxA + 1;
    // Compute the momentum provided by config. A
    La = 0;
    for (i = 0; i < Nla; i++)
    {
        La = La + occA[i] * (i - S->lmaxA);
    }
    // find the subspace index
    lower = 0;
    upper = S->n_sub;
    n = (upper + lower) / 2;
    while(La != S->sub[n].La)
    {
        if (La > S->sub[n].La) lower = n;
        else                   upper = n;
        n = (upper + lower)/2;
    }
    return n;
}



int BBgetIndex(CompoundSpace S, Iarray occA, Iarray occB)
{

/** Return the index of the tensor product of 'occA' with 'occB' **/

    int
        i,
        n,
        ia,
        ib,
        La,
        Lb,
        Nla,
        Nlb,
        upper,
        lower;

    Nla = 2 * S->lmaxA + 1;
    Nlb = 2 * S->lmaxB + 1;
    // assert the 'occA' and 'occB' belong to the input multiconfig. space
    if (!assert_BoseNocc(Nla,S->Na,occA))
    {
        printf("\n\nERROR: In operations computing Hamiltonian action the ");
        printf("following occupation numbers for Bosons are not valid \n");
        for (i = 0; i < Nla; i++) printf(" %d",occA[i]);
        printf("\n\nProgram aborted at function 'BBgetIndex'\n\n");
        exit(EXIT_FAILURE);
    }
    if (!assert_BoseNocc(Nlb,S->Nb,occB))
    {
        printf("\n\nERROR: In operations computing Hamiltonian action the ");
        printf("following occupation numbers for Bosons are not valid \n");
        for (i = 0; i < Nlb; i++) printf(" %d",occB[i]);
        printf("\n\nProgram aborted at function 'BBgetIndex'\n\n");
        exit(EXIT_FAILURE);
    }

    // Compute the momentum provided by config. A
    La = 0;
    for (i = 0; i < Nla; i++)
    {
        La = La + occA[i] * (i - S->lmaxA);
    }
    // Compute the momentum provided by config. B
    Lb = 0;
    for (i = 0; i < Nlb; i++)
    {
        Lb = Lb + occB[i] * (i - S->lmaxB);
    }
    // Assert the combined momentum give the correct result
    if (La + Lb != S->L)
    {
        printf("\n\nERROR : in function BBgetIndex the configurations ");
        printf("provided break the momentum conservation\n\n");
        exit(EXIT_FAILURE);
    }

    // First find the subspace index
    lower = 0;
    upper = S->n_sub;
    n = (upper + lower) / 2;
    while(La != S->sub[n].La)
    {
        if (La > S->sub[n].La) lower = n;
        else                   upper = n;
        n = (upper + lower)/2;
    }

    // Find the index in each individual hashing table
    ia = BgetIndex(S->lmaxA,S->sub[n].sizeA,S->sub[n].hta,occA);
    ib = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,occB);

    return S->strides[n] + ib*S->sub[n].sizeA + ia;
}



int BBgetIndex_A(CompoundSpace S, int n, int ib, Iarray occA)
{

/** Return the index of configuration using the  occupation  numbers
    of species A with the subspace  index  'n'  and  index  of  some
    configuration for species B. Useful for operators that act  only
    in space A and conserves separately the momentum of both species **/

    int
        ia;

    // Find the index in A hashing table
    ia = BgetIndex(S->lmaxA,S->sub[n].sizeA,S->sub[n].hta,occA);
    return S->strides[n] + ib*S->sub[n].sizeA + ia;
}



int BBgetIndex_B(CompoundSpace S, int n, int ia, Iarray occB)
{

/** Similar to 'BBgetIndex_A' **/

    int
        ib;

    // Find the index in B hashing table
    ib = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,occB);
    return S->strides[n] + ib*S->sub[n].sizeA + ia;
}



void PrintConfig(CompoundSpace S)
{

    int
        i,
        j,
        Ma,
        Mb;

    Iarray
        occA,
        occB;

    Ma = 2 * S->lmaxA + 1;
    Mb = 2 * S->lmaxB + 1;
    occA = iarrDef(Ma);
    occB = iarrDef(Mb);

    printf("\n\nConfig. Index   [-lmaxA, ..., lmaxA] x [-lmaxB, ..., lmaxB]");
    printf("\n___________________________________________________________\n");

    for (i = 0; i < S->size; i++)
    {
        printf("\n%8d  ",i);
        BBgetConfigs(i,S,occA,occB);
        printf("[ ");
        for (j = 0; j < Ma; j++) printf("%2d ",occA[j]);
        printf(" ] x ");
        printf("[ ");
        for (j = 0; j < Mb; j++) printf("%2d ",occB[j]);
        printf(" ] ");
        j = BBgetIndex(S,occA,occB);
        printf("%d",j);
        if (i != j)
        {
            printf("\n\nERROR : The Hashing function is wrong ! ");
            printf("Contact developer\n\n");
            exit(EXIT_FAILURE);
        }
    }
}



#endif
