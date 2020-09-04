#include <time.h>
#include "Hamiltonian.h"



int NaiveSetup(unsigned int Npar, unsigned int Morb, int L)
{

/** setup the hashing table of the multiconf. space with constrained
    momentum 'L' in a naive way, which first assemble the hashing
    table of the unconstrained (arbitrary 'L') problem. The total
    number of orbitals is Morb = 2*lmax+1 where 'lmax' is the max.
    orbital momentum, assuming their ordering -lmax , ... ,  lmax
 
    RETURN : size of constrained multiconfig. space with momentum 'L' **/

    int
        i,
        n,
        nc,
        count,
        totalMom;

    Iarray
        gen_config,
        HTindexes,
        ht;

    nc = NC(Npar,Morb); // size of entire config. space
    // Array to mark the config. indexes that have the momentum demanded
    HTindexes = iarrDef(nc);
    // general configurations without restrictions
    gen_config = iarrDef(Morb);

    count = 0;
    for (n = 0; n < nc; n++)
    {
        // scan the entire configurational space without restriction
        // selecting configurations that has the momentum demanded L
        indexToConfig(n,Npar,Morb,gen_config);
        totalMom = 0;
        for (i = 0; i < Morb; i++)
        {
            totalMom = totalMom + (i - Morb/2)*gen_config[i];
        }
        // Update the indexes of Fock states which have momentum L
        if (totalMom == L)
        {
            // record index whose Fock state has momentum L
            HTindexes[count] = n;
            count++;
        }
    }

    // Setup config.  with required momentum L
    // previously marked by index in HTindexes
    ht = iarrDef(count*Morb);
    for (n = 0; n < count; n++)
    {
        indexToConfig(HTindexes[n],Npar,Morb,&ht[n*Morb]);
    }

    free(HTindexes); free(ht); free(gen_config);
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
        Npar,   // number of particles
        lmax,   // number of orbitals
        totalL, // Momentum demanded for the config. space
        mcSize; // Size of config. space constrained

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
    printf("\nWorking on improved set up of L = %d space ",totalL);
    printf(" ."); printf("."); printf(".");
    start = clock(); // trigger to measure time
    mcSize = BFixedMom_mcsize(Npar,lmax,totalL);
    printf("..");
    printf("..");
    ht = BAssembleHT(Npar,lmax,totalL,mcSize);
    end = clock();   // finish time measure
    time_direct = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf(" Done !");

    // COMPUTE IN A NAIVE WAY WITHOUT RESTRICTIONS
    nc = NC(Npar,2*lmax+1); // Size of config. without restrictions
    if (nc*sizeof(int) < MEMORY_TOL)
    {
        printf("\nWorking on naive set up of L = %d space ",totalL);
        printf(" ...");
        start = clock(); // trigger to measure time
        k = NaiveSetup(Npar,2*lmax+1,totalL);
        end = clock();   // finish time measure
        time_naive = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf(".... Done !");

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
        printf("being evaluated because it exceeds memory tolerance ");
        printf("( > %.1lf )\n\n",((double) MEMORY_TOL) / 1E9);
    }



    // SHOW ON THE SCREEN SOME CONFIG. IF THERE ARE NOT TOO MANY
    if (mcSize < 500 && mcSize > 0)
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

    if (nnz > 0)
    {
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
            printf("\n\nFATAL ERROR : FUNCTION OF NONZERO ENTRIES\n\n");
            exit(EXIT_FAILURE);
        }
    }

    printf("\n\n\n");
    printf("\t******************************************************\n");
    printf("\t*                                                    *\n");
    printf("\t*        SOME MULTICONFIG. SPACE INFORMATIONS        *\n");
    printf("\t*                                                    *\n");
    printf("\t******************************************************\n");
    printf("\n\n\n");

    printf("There are %d Fock states with L = %d",mcSize,totalL);
    printf("\nthere are %d Fock states with arbitrary L",nc);
    if (nnz > 0)
    {
        printf("\nNumber of nonzero entries in H : %d | ",nnz);
        printf("Spartisty = %.5lf",1.0 - ((double) nnz)/mcSize/mcSize);
        printf("\nMaximum number of nonzero entries in a same row : %d",k);
    }
    else
    {
        printf("\nSparse matrix Hamiltonian would exceed allowed memory ");
        printf("%.1lf(GB)",MEMORY_TOL/1E9);
    }
    printf("\n\nMemory required for the:\n");
    l = ((double) mcSize*(2*lmax+1)*sizeof(int))/1E6;
    printf("\tmulticonfig. space constraining the momentum : %.1lf(Mb)\n",l);
    l = ((double) mcSize*(2*lmax+1) + nc)*sizeof(int)/1E6;
    printf("\tmulticonfig. space naively using all config. : %.1lf(Mb)\n",l);
    l = ((double) mcSize * sizeof(double complex))/1E6;
    printf("\tA many-body state in config. basis : %.1lf(Mb)\n",l);
    if (nnz > 0)
    {
        l = ((double)(nnz+mcSize)*sizeof(int)+nnz*sizeof(double complex))/1E6;
        printf("\tsparse 'Hamiltonian' matrix : %.1lf(Mb)\n",l);
    }

    // TIME REQUIRED IN EACH IMPORTANT DATA STRUCTURE SETUP
    printf("\n\nTime to setup multiconfig. constraining the momentum : ");
    if (time_direct < 0.1) printf("%.1lf(ms)",time_direct*1000);
    else printf("%.1lf(s)",time_direct);

    printf("\nTime to setup multiconfig. naively using all config. : ");
    if (time_naive < 0.1) printf("%.1lf(ms)",time_naive*1000);
    else printf("%.1lf(s)",time_naive);

    printf("\n\nDone.\n\n");

    for (i = 0; i < mcSize; i++) free(ht[i]);
    free(ht);
    free(NNZrow);

    return 0;
}
