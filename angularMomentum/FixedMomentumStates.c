#include <time.h>
#include "HMatrixSetup.h"

#define PI 3.141592653589793



int NaiveSetup(unsigned int Npar, unsigned int Morb, int L)
{

/** setup the hashing table of the multiconf. space with constrained
    momentum 'L' in a naive way, which first assemble the hashing
    table of the unconstrained (arbitrary 'L') problem. The total
    number of orbitals is Morb = 2*lmax+1 where 'lmax' is the max.
    orbital momentum, assuming them order as  -lmax , ... ,  lmax
 
    RETURN : size of constrained multiconfig. space with momentum 'L' **/

    int
        i,
        n,
        nc,
        count,
        totalMom;

    Iarray
        general_ht,
        HTindexes,
        ht;

    count = 0;
    nc = NC(Npar,Morb);
    HTindexes = iarrDef(nc);

    // first setup a hashing table for the general problem without
    // fixing the angular momentum required in the argument 'L'
    general_ht = setupFocks(Npar,Morb);

    for (n = 0; n < nc; n++)
    {
        totalMom = 0;
        for (i = 0; i < Morb; i++)
        {
            totalMom = totalMom + (i - Morb/2)*general_ht[i+n*Morb];
        }

        // Update the numbering of Fock states
        // which satisfy the momentum required
        if (totalMom == L)
        {
            // record index whose Fock state has momentum 'L'
            HTindexes[count] = n;
            count++;
        }
    }

    ht = iarrDef(count*Morb);

    // Copy the configurations that satisfy the momentum required
    for (n = 0; n < count; n++)
    {
        for (i = 0; i < Morb; i++)
        {
            ht[i + n*Morb] = general_ht[i + HTindexes[n]*Morb];
        }
    }

    free(HTindexes); free(ht); free(general_ht);
    return count;
}



int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        nc,
        nnz,
        Npar, // number of particles
        lmax, // number of orbitals
        totalL,
        mcSize;

    double
        l,
        time_H,
        time_naive,
        time_direct;

    clock_t
        start,
        end;

    Iarray
        NNZrow,
        * ht;

    Carray
        Ho;

    HConfMat
        H;

    if (argc != 4)
    {
        printf("\n\nERROR: Need three integer numbers from command line");
        printf("\n\t1. Number of particles");
        printf("\n\t2. max. IPS angular momentum");
        printf("\n\t3. Total angular momentum\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);     // Total number of particles
    sscanf(argv[2],"%d",&lmax);     // maximum momentum of IPS
    sscanf(argv[3],"%d",&totalL);   // Constrained momentum of config.

    // COMPUTE DIRECTLY THE SPACE WITH MOMENTUM 'totalL'
    start = clock(); // trigger to measure time
    mcSize = BFixedMom_mcsize(Npar,lmax,totalL);
    ht = BAssembleHT(Npar,lmax,totalL,mcSize);
    end = clock();   // finish time measure
    time_direct = ((double) (end - start)) / CLOCKS_PER_SEC;

    // COMPUTE IN A NAIVE WAY FIRST ASSEMBLING THE
    // THE HASHING TABLE WITHOUT RESTRICTIOS
    nc = NC(Npar,2*lmax+1); // Size of config. without restrictions
    if (nc*(2*lmax+1)*sizeof(int) < 1E9)
    {
        start = clock(); // trigger to measure time
        k = NaiveSetup(Npar,2*lmax+1,totalL);
        end = clock();   // finish time measure
        time_naive = ((double) (end - start)) / CLOCKS_PER_SEC;

        // Self-consistency check between methos to generate config.
        if (k != mcSize)
        {
            printf("\n\nERROR : The naive method does not agree with ");
            printf("the direct one in the multiconfig. size. ");
            printf("Critical routine error\n\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        printf("\n\nWARNING : The naive setup of Conf. space is not ");
        printf("being evaluated because it exceeds memory tolerance to ");
        printf("store the general hashing table ( > 1GB )\n\n");
    }



    // SHOW ON THE SCREEN SOME CONFIG.
    if (mcSize < 1000)
    {
        printf("\n\nConfigurations with ang. momentum L = %d and",totalL);
        printf(" max(l) = %d\n",lmax);
        printf("=========================================================\n\n");

        printf("conf_index  | ");
        for (i = 0; i < 2*lmax; i++) printf("%4d  , ",i-lmax);
        printf("%4d  | get_index",i-lmax);
        printf("\n");
        for (i = 0; i < mcSize; i++)
        {
            printf("\n %6d     [ ",i+1);
            for (j = 0; j < 2*lmax; j++) printf("%4d  , ",ht[i][j]);
            printf("%4d  ]",ht[i][j]);

            k = BgetIndex(lmax,mcSize,ht,ht[i]);
            if (i != k)
            {
                printf("\n\nFATAL ERROR : FUNCTION TO GET INDEX IS WRONG !!");
                printf("\n\n");
                exit(EXIT_FAILURE);
            }
            printf("  %6d",k+1);
        }

        printf("\n\n=========================================================");

    }

    NNZrow = iarrDef(mcSize); // Number of non-zero entries per row
    nnz = NNZ_PerRow(Npar,lmax,mcSize,ht,NNZrow);

    k = 0;
    j = 0;
    // Extract maximum number of non-zero entries in a same row
    for (i = 0; i < mcSize; i++)
    {
        if (k < NNZrow[i]) k = NNZrow[i];
        j = j + NNZrow[i];
    }
    // Self-consistency check for function that
    // computes number of non-zero entries in H
    if (j != nnz)
    {
        printf("\n\nFATAL ERROR : FUNCTION TO NONZERO ENTRIES IS WRONG\n\n");
        exit(EXIT_FAILURE);
    }

    printf("\n\nThere are %d Fock states with L = %d",mcSize,totalL);
    printf("\nthere are %d Fock states with arbitrary L",nc);
    printf("\nNumber of nonzero entries in H : %d | ",nnz);
    printf("Spartisty = %.5lf",1.0 - ((double) nnz)/mcSize/mcSize);
    printf("\nMaximum number of nonzero entries in a same row : %d",k);
    printf("\n\nMemory required for the:\n");
    printf("\tmulticonfig. space constraining the momentum : %.1lf(Mb)\n",
            ((double) mcSize*(2*lmax+1)*sizeof(int))/1E6);
    printf("\tmulticonfig. space naively using all config. : %.1lf(Mb)\n",
            ((double) (mcSize + nc)*(2*lmax+1) + nc)*sizeof(int)/1E6);
    l = ((double) (nnz+mcSize)*sizeof(int) + nnz*sizeof(double complex));
    printf("\tto store sparse 'Hamiltonian' matrix : %.1lf(Mb)\n",l/1E6);

    Ho = carrDef(2*lmax+1);
    for (i = 0; i < 2*lmax+1; i++)
    {
        l = 2 * PI * (i-lmax);
        Ho[i] = 0.5*l*l;
    }

    start = clock(); // trigger to measure time
    H = assembleH(Npar,lmax,mcSize,ht,Ho,1.0);
    end = clock();   // finish time measure
    time_H = ((double) (end - start)) / CLOCKS_PER_SEC;

    // TIME REQUIRED IN EACH IMPORTANT DATA STRUCTURE SETUP
    printf("\n\nTime to setup multiconfig. constraining the momentum : ");
    if (time_direct < 0.1) printf("%.1lf(ms)",time_direct*1000);
    else printf("%.1lf(s)",time_direct);

    printf("\nTime to setup multiconfig. naively using all config. : ");
    if (time_naive < 0.1) printf("%.1lf(ms)",time_naive*1000);
    else printf(" %.1lf(s)",time_naive);

    printf("\nTime to setup H. matrix : ");
    if (time_H < 0.1)
    {
        printf("%.1lf(ms)",time_H*1000);
    }
    else
    {
        if (time_H < 60) printf("%.1lf(s)",time_H);
        else             printf("%.1lf(min)",time_H/60.0);
    }

    printf("\n\nDone.\n\n");

    for (i = 0; i < mcSize; i++) free(ht[i]);
    free(ht);
    free(Ho);
    free(NNZrow);
    freeHmat(H);

    return 0;
}
