#include "LBoseFockSpace.h"



struct _TensorProd_ht
{
    int
        La,     // momentum of config. space A
        Lb,     // momentum of config. space B
        sizeA,  // multiconfig. space size A
        sizeB;  // multiconfig. space size B

    Iarray
        * hta,  // hashing table for config. in space A
        * htb;  // hashing table for config. in space B
};

typedef struct _TensorProd_ht * TensorProd_ht;

struct _CompoundSpace
{

/** The compound space is design by selecting the momenta in  integer
    units separately for species A and B, such that  La + Lb = L.  In
    this way, it is possible to define an array of subspaces for each
    possible combination of angular momenta. **/

    int
        L,     // total angular momentum (from both species)
        size,  // total number of elements in compound basis
        n_sub, // number of subspaces - size of arrays below
        lmaxA, // max. momentum for IPS for species A
        lmaxB; // max. momentum for IPS for species B

    // In a contiguous enumeration of the product space,
    // the strides mark in which index  a  new  subspace
    // begins, which are related to momenta of  A and  B
    Iarray
        strides;

    struct _TensorProd_ht
        * sub; // subspaces all possible Tensor product spaces
};

typedef struct _CompoundSpace * CompoundSpace;



CompoundSpace AllocCompBasis(int Na, int Nb, int lmaxA, int lmaxB, int L)
{

/** Na : Number of particles of species A
    Nb : Number of particles of species B
    lmaxA : maximum momentum for species A
    lmaxB : maximum momentum for species B
    L : Total (combined) momentum from both species **/

    int
        i,
        j,
        La,
        Lb,
        n_sub,
        size,
        MemReq;

    Iarray
        nc_a,
        nc_b;

    struct _CompoundSpace
        * S;

    S = (CompoundSpace) malloc(sizeof(struct _CompoundSpace));
    S->L = L;
    S->lmaxA = lmaxA;
    S->lmaxB = lmaxB;

    // array to store the number of config.(size of individual basis)
    // for each species when sweeping through the possible values  of
    // the momentum. For bosons the values for 'La'  are in the range
    // - Na * lmax , ... , 0 , ... , Na * lmax. Once 'La' is fixed
    // we must have Lb = L - La
    nc_a = iarrDef(2*Na*lmaxA+1);
    nc_b = iarrDef(2*Na*lmaxA+1);

    i = 0;
    size = 0;
    n_sub = 0;
    MemReq = 0;
    for (La = -Na*lmaxA; La <= Na*lmaxA; La++)
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
            size = size + nc_a[i]*nc_b[i]; // size of compound basis updated

            if (size > INT_MAX - nc_a[i]*nc_b[i])
            {
                printf("\n\nINDEX ERROR : the structure for the basis ");
                printf("of the compound multiconfig. basis cannot be ");
                printf("indexed with 32-bit integers\n\n");
                printf("Program aborted at function 'AllocCompBasis' ");
                printf("in file 'mixtures.h'\n\n");
                exit(EXIT_FAILURE);
            }

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
    for (La = -Na*lmaxA; La <= Na*lmaxA; La++)
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



unsigned int EstimateMemory(CompoundSpace S)
{

/** Compute approximately the memory required to setup the config. space **/

    unsigned int
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

/** Setup in 'occA' and 'occB' the occupation numbers corresponding to 'k' **/

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
    subIndex = j-1;
    k = k - S->strides[j-1];

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



int BBgetIndex(CompoundSpace S, Iarray occA, Iarray occB)
{

/** Return the index of the tensor product of 'occA' with 'occB' **/

    int
        i,
        n,
        ia,
        ib,
        La,
        Nla,
        upper,
        lower;

    Nla = 2 * S->lmaxA + 1;

    La = 0;
    for (i = 0; i < Nla; i++)
    {
        La = La + occA[i] * (i - S->lmaxA);
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
            printf("\n\nERROR : The Hashing function are wrong ! ");
            printf("Contact developer\n\n");
            exit(EXIT_FAILURE);
        }
    }
}