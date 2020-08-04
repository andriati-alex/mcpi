#ifndef _LBoseFermiFockSpace_h
#define _LBoseFermiFockSpace_h

#include "LBoseFockSpace.h"
#include "LFermiFockSpace.h"



BFCompoundSpace BoseFermiBasis(int Nb, int Nf, int lmaxB, int lmaxF, int L)
{

/** Nb : Number of bosons
    Nf : Number of fermions
    lmaxB : maximum single particle momentum for bosons
    lmaxF : maximum single particle momentum for fermions
    L : Total (combined) momentum from both species
    --------------------------------------------------------------------
    The compound basis is build assembling an array of subspaces,  which
    are characterized by one of the possible combination of  the  spaces
    with fixed momentum Lb and Lf for each species, such that the sum of
    both give us L. The subspaces are the tensor product of these  fixed
    momentum spaces for every possible combination of Lb and Lf.     **/

    int
        i,
        j,
        Lb,
        Lf,
        size,
        n_sub,
        LboundB;
    long
        MemReq;
    Iarray
        nc_b,
        nc_f;
    struct _BFCompoundSpace
        * S;

    S = (BFCompoundSpace) malloc(sizeof(struct _BFCompoundSpace));
    S->L = L;
    S->Nb = Nb;
    S->Nf = Nf;
    S->lmaxB = lmaxB;
    S->lmaxF = lmaxF;
    LboundB = Nb*lmaxB;

    // array to store the number of config.(size of individual basis)
    // for each species when sweeping through the possible values  of
    // the momentum. For bosons the values for 'Lb'  are in the range
    // - Nb * lmax , ... , 0 , ... , Nb * lmax. Once 'Lb' is fixed
    // we must have Lf = L - Lb
    nc_b = iarrDef(2*LboundB+1);
    nc_f = iarrDef(2*LboundB+1);

    i = 0;
    size = 0;
    n_sub = 0;
    MemReq = 0;
    // COMPUTE THE SIZE OF THE SPACE  ENSURING
    // THAT IT DO NOT EXCEED SYSTEM CAPABILITY
    for (Lb = -LboundB; Lb <= LboundB; Lb++)
    {
        // Search for the possible combinations of momenta of spaces B and F
        Lf = L-Lb;
        nc_b[i] = BFixedMom_mcsize(Nb,lmaxB,Lb);
        nc_f[i] = FFixedMom_mcsize(Nf,lmaxF,Lf);
        if (nc_b[i]*nc_f[i] > 0)
        {
            // When the separate spaces with momenta 'Lb' and 'Lf'
            // coexist, they configure a subspace of the compound
            // system through the tensor product
            n_sub++; // Add counter of subspaces found

            // Assess if it is possible the enumeration of the basis
            if (size > INT_MAX - nc_b[i]*nc_f[i])
            {
                printf("\n\nINDEX ERROR : the structure for the basis ");
                printf("of the compound multiconfig. basis cannot be ");
                printf("indexed by 32-bit integers\n\n");
                printf("Program aborted at function 'AllocCompBasis' ");
                printf("in file 'LBoseBoseFockSpace.h'\n\n");
                exit(EXIT_FAILURE);
            }
            size = size + nc_b[i]*nc_f[i]; // size of compound basis updated

            MemReq += sizeof(struct _BFTensorProd_ht);
            MemReq += (nc_b[i]*(2*lmaxB+1) + nc_f[i]*(2*lmaxF+1))*sizeof(int);
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

    S->sub = (BFTensorProd_ht) malloc(n_sub*sizeof(struct _BFTensorProd_ht));
    S->strides = iarrDef(n_sub);
    S->n_sub = n_sub;   // Number of subspaces each one with 'La/Lb' fixed
    S->size = size;     // number of config. for compound basis

    j = 0;
    i = 0;
    size = 0;
    for (Lb = -LboundB; Lb <= LboundB; Lb++)
    {
        Lf = L-Lb;
        if (nc_b[i]*nc_f[i] > 0)
        {
            // set up the subspace
            S->strides[j] = size;
            S->sub[j].Lb = Lb;
            S->sub[j].Lf = Lf;
            S->sub[j].sizeB = nc_b[i];
            S->sub[j].sizeF = nc_f[i];
            S->sub[j].htb = BAssembleHT(Nb,lmaxB,Lb,nc_b[i]);
            S->sub[j].htf = FAssembleHT(Nf,lmaxF,Lf,nc_f[i]);
            size = size + nc_b[i]*nc_f[i];
            j++;
        }
        i++;
    }

    free(nc_b);
    free(nc_f);

    return S;
}



unsigned long EstimateMemoryBF(BFCompoundSpace S)
{

/** Compute approximately the memory required to set up the compound space **/

    unsigned long
        j,
        sub_mem,
        mem_req,
        lmaxB = S->lmaxB,
        lmaxF = S->lmaxF;

    sub_mem = 0;
    for (j = 0; j < S->n_sub; j++)
    {
        sub_mem = sub_mem + S->sub[j].sizeB*(2*lmaxB+1)*sizeof(int);
        sub_mem = sub_mem + S->sub[j].sizeF*(2*lmaxF+1)*sizeof(char);
    }

    mem_req = sub_mem + S->n_sub*(sizeof(struct _BFTensorProd_ht)+sizeof(int));
    return mem_req;
}



void freeBoseFermiSpace(BFCompoundSpace S)
{
    int
        j,
        i;

    free(S->strides);
    // Free the Hashing table for each subspace
    for (j = 0; j < S->n_sub; j++)
    {
        for (i = 0; i < S->sub[j].sizeB; i++) free(S->sub[j].htb[i]);
        free(S->sub[j].htb);
        for (i = 0; i < S->sub[j].sizeF; i++) free(S->sub[j].htf[i]);
        free(S->sub[j].htf);
    }
    free(S->sub);
    free(S);
}



void BFgetConfigs(int k, BFCompoundSpace S, Iarray occB, Farray occF)
{

/** Set up in 'occA' and 'occB' the occupation numbers corresponding to
    configuration 'k' of the compound space 'S'                     **/

    int
        i,
        j,
        Lb,
        Lf,
        Nlb,
        Nlf,
        indexB,
        indexF,
        subIndex,
        SubSizeB,
        SubSizeF;

    if (k > S->size)
    {
        printf("\n\nERROR : The index %d of the configurations ",k);
        printf("exceed the size of the compound space %d\n\n",S->size);
        exit(EXIT_FAILURE);
    }

    // Total number of IPS for each species
    Nlb = 2 * S->lmaxB + 1;
    Nlf = 2 * S->lmaxF + 1;

    // First, find the index of the subspace comparing to strides
    // related to number of states in each subspace
    for (j = 1; j < S->n_sub; j++)
    {
        if (k < S->strides[j]) break;
    }
    subIndex = j-1;             // subspace number
    k = k - S->strides[j-1];    // subtract cost of accessing the subspace

    SubSizeB = S->sub[subIndex].sizeB;
    SubSizeF = S->sub[subIndex].sizeF;
    Lb = S->sub[subIndex].Lb;
    Lf = S->sub[subIndex].Lf;

    // Assert k has the right range in the given subspace
    // it must be in the interval [0,sizeA*sizeB)
    if (k >= SubSizeB*SubSizeF)
    {
        printf("\n\nERROR : wrong subspace index found %d - ",k);
        printf("Subspace with Lb = %d and Lf = %d ",Lb,Lf);
        printf("has %d elements\n\n",SubSizeB*SubSizeF);
        exit(EXIT_FAILURE);
    }

    indexF = k / SubSizeB;
    indexB = k % SubSizeB;

    // copy the occupations to the output parameters
    for (i = 0; i < Nlb; i++) occB[i] = S->sub[subIndex].htb[indexB][i];
    for (i = 0; i < Nlf; i++) occF[i] = S->sub[subIndex].htf[indexF][i];
}



int BFsubIndex(BFCompoundSpace S, Iarray occB)
{

/** From occupation numbers in 'occB' compute the angular momentum 'Lb'
    and return the index of the subspace. The subspaces are sorted  as
    increasing values of 'Lb' (bosons angular momentum)            **/

    int
        i,
        n,
        Lb,
        Nlb,
        upper,
        lower;

    Nlb = 2 * S->lmaxB + 1;
    // Compute the momentum provided by config. A
    Lb = 0;
    for (i = 0; i < Nlb; i++)
    {
        Lb = Lb + occB[i] * (i - S->lmaxB);
    }
    // find the subspace index
    lower = 0;
    upper = S->n_sub;
    n = (upper + lower) / 2;
    while(Lb != S->sub[n].Lb)
    {
        if (Lb > S->sub[n].Lb) lower = n;
        else                   upper = n;
        n = (upper + lower)/2;
    }
    return n;
}



int BFgetIndex(BFCompoundSpace S, Iarray occB, Farray occF)
{

/** Return the index of the tensor product of 'occB' with 'occF' **/

    int
        i,
        n,
        indexB,
        indexF,
        Lb,
        Lf,
        Nlb,
        Nlf,
        upper,
        lower;

    Nlb = 2 * S->lmaxB + 1;
    Nlf = 2 * S->lmaxF + 1;
    if (!assert_FermiNocc(Nlf,S->Nf,occF))
    {
        printf("\n\nERROR: In operations computing Hamiltonian action the ");
        printf("following occupation numbers for Fermions are not valid \n");
        for (i = 0; i < Nlf; i++) printf(" %d",occF[i]);
        printf("\n\nProgram aborted at function 'BFgetIndex'\n\n");
        exit(EXIT_FAILURE);
    }
    if (!assert_BoseNocc(Nlb,S->Nb,occB))
    {
        printf("\n\nERROR: In operations computing Hamiltonian action the ");
        printf("following occupation numbers for Bosons are not valid \n");
        for (i = 0; i < Nlb; i++) printf(" %d",occB[i]);
        printf("\n\nProgram aborted at function 'BFgetIndex'\n\n");
        exit(EXIT_FAILURE);
    }

    // Compute the momentum provided by config. A
    Lb = 0;
    for (i = 0; i < Nlb; i++)
    {
        Lb = Lb + occB[i] * (i - S->lmaxB);
    }
    // Compute the momentum provided by config. B
    Lf = 0;
    for (i = 0; i < Nlf; i++)
    {
        Lf = Lf + occF[i] * (i - S->lmaxF);
    }
    // Assert the combined momentum give the correct result
    if (Lb + Lf != S->L)
    {
        printf("\n\nERROR : in function BFgetIndex the configurations ");
        printf("provided break the momentum conservation\n\n");
        exit(EXIT_FAILURE);
    }

    // First find the subspace index
    lower = 0;
    upper = S->n_sub;
    n = (upper + lower) / 2;
    while(Lb != S->sub[n].Lb)
    {
        if (Lb > S->sub[n].Lb) lower = n;
        else                   upper = n;
        n = (upper + lower)/2;
    }

    // Find the index in each individual hashing table
    indexB = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,occB);
    indexF = FgetIndex(S->lmaxF,S->sub[n].sizeF,S->sub[n].htf,occF);

    return S->strides[n] + indexF*S->sub[n].sizeB + indexB;
}



int BFgetIndex_B(BFCompoundSpace S, int n, int indexF, Iarray occB)
{

/** Return the index of configuration using the occupation numbers
    of bosons with the subspace index 'n' and index a Fermionic
    configuration indexF. Useful for operators that act only on
    bosons and conserves separately the momentum of both species **/

    int
        indexB;

    // Find the index in A hashing table
    indexB = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,occB);
    return S->strides[n] + indexF*S->sub[n].sizeB + indexB;
}



int BFgetIndex_F(BFCompoundSpace S, int n, int indexB, Farray occF)
{

/** Similar to 'BFgetIndex_B' but using a config. of fermions 'occF' **/

    int
        indexF;

    // Find the index in B hashing table
    indexF = FgetIndex(S->lmaxF,S->sub[n].sizeF,S->sub[n].htf,occF);
    return S->strides[n] + indexF*S->sub[n].sizeB + indexB;
}



void PrintConfigBF(BFCompoundSpace S)
{

    int
        i,
        j,
        Mb,
        Mf;
    Iarray
        occB;
    Farray
        occF;

    Mb = 2 * S->lmaxB + 1;
    Mf = 2 * S->lmaxF + 1;
    occB = iarrDef(Mb);
    occF = farrDef(Mf);

    printf("\n\nConfig. Index   [-lmaxB, ..., lmaxB] x [-lmaxF, ..., lmaxF]");
    printf("\n___________________________________________________________\n");

    for (i = 0; i < S->size; i++)
    {
        printf("\n%8d  ",i);
        BFgetConfigs(i,S,occB,occF);
        printf("[ ");
        for (j = 0; j < Mb; j++) printf("%2d ",occB[j]);
        printf(" ] x ");
        printf("[ ");
        for (j = 0; j < Mf; j++) printf("%d ",occF[j]);
        printf("] ");
        j = BFsubIndex(S,occB);
        printf(" Lbos = %4d  -  LFer = %4d",S->sub[j].Lb,S->sub[j].Lf);
        j = BFgetIndex(S,occB,occF);
        if (i != j)
        {
            printf("\n\nERROR : The Hashing function is wrong ! ");
            printf("Contact developer\n\n");
            exit(EXIT_FAILURE);
        }
    }
}



#endif
