#include <time.h>
#include "LFermiFockSpace.h"



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

    Farray
        occ,
        ht;

    Iarray
        HTindexes;

    count = 0;
    nc = NCF(Npar,Morb);
    occ = farrDef(Morb);
    HTindexes = iarrDef(nc);

    for (n = 0; n < nc; n++)
    {
        FindexToConfig(n,Npar,Morb,occ);
        totalMom = 0;
        for (i = 0; i < Morb; i++)
        {
            totalMom = totalMom + (i - Morb/2)*occ[i];
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

    ht = farrDef(count*Morb);
    // Copy the configurations that satisfy the momentum required
    for (n = 0; n < count; n++)
    {
        FindexToConfig(HTindexes[n],Npar,Morb,occ);
        for (i = 0; i < Morb; i++) ht[i+n*Morb] = occ[i];
    }

    free(HTindexes); free(ht); free(occ);
    return count;
}



int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        nc,
        Npar,
        lmax,
        totalL,
        mcSize;

    double
        time_naive,
        time_direct;

    clock_t
        start,
        end;

    Farray
        * ht;

    if (argc != 4)
    {
        printf("\n\nERROR: Need three integer numbers from command line");
        printf("\n\t1. Number of particles");
        printf("\n\t2. max. IPS angular momentum");
        printf("\n\t3. Total angular momentum\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);     // Number of particles
    sscanf(argv[2],"%d",&lmax);     // Maximum ang. mom. of the IPS
    sscanf(argv[3],"%d",&totalL);   // Total ang. mom. of the Fock states

    if (2*lmax+1 < Npar)
    {
        printf("\n\nERROR: The number of fermions (%d) cannot exceed ",Npar);
        printf("the number of orbitals to allocate them (%d) ",2*lmax+1);
        printf("according to Pauli exclusion principle\n\n");
        exit(EXIT_FAILURE);
    }

    // COMPUTE USING THE METHOD TO DIRECTLY SETUP
    // CONFIG. WHICH SATISFY THE REQUIRED 'totalL'
    start = clock(); // trigger to measure time
    mcSize = FFixedMom_mcsize(Npar,lmax,totalL);
    ht = FAssembleHT(Npar,lmax,totalL,mcSize);
    end = clock();   // finish time measure
    time_direct = ((double) (end - start)) / CLOCKS_PER_SEC;

    // COMPUTE IN A NAIVE WAY FIRST ASSEMBLING THE
    // THE HASHING TABLE WITHOUT RESTRICTIOS
    nc = NCF(Npar,2*lmax+1); // Size of config. without restrictions
    if (nc*sizeof(int) < 1E9)
    {
        start = clock(); // trigger to measure time
        k = NaiveSetup(Npar,2*lmax+1,totalL);
        end = clock();   // finish time measure
        time_naive = ((double) (end - start)) / CLOCKS_PER_SEC;

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



    // IF THE MULTICONF. SPACE IS NOT TOO LARGE PRINT THE CONFIG.
    if (mcSize < 5000)
    {
        if (2*lmax+1 > 7)
        {
            printf("\n\nWARNING: MANY SINGLE PARTICLE STATES -> ");
            printf("COMPACT PRINTING\n\n");

            printf(" Conf_index    Occupations as binary");
            for (i = 0; i < mcSize; i++)
            {
                printf("\n%9d      [ ", i);

                for (j = 0; j < 2*lmax+1; j++) printf("%d",ht[i][j]);
                printf(" ]");

                k = FgetIndex(lmax,mcSize,ht,ht[i]);
                printf(" %5d",k);
                if (k != i)
                {
                    printf("\n\nERROR : The routine to search conf. ");
                    printf("in hashing table gave wrong index.\n\n");
                    exit(EXIT_FAILURE);
                }

                k = 0;
                for (j = 0; j < 2*lmax+1; j++) k = k + ht[i][j];
                if (k != Npar)
                {
                    // Self consistency check, if the Fock state has
                    // the total number of particles supplied
                    printf("\n\nERROR: Number of particles in the Fock ");
                    printf("state above is different from total number of ");
                    printf("particles given %d\n\n",Npar);
                    exit(EXIT_FAILURE);
                }
            }
        }

        else
        {
            printf("\n\n");
            printf("conf_index  | ");
            for (i = 0; i < 2*lmax; i++) printf("%4d  , ",i-lmax);
            printf("%4d  | get_index",i-lmax);
            printf("\n");
            for (i = 0; i < mcSize; i++)
            {
                printf("\n %6d     [ ",i);
                for (j = 0; j < 2*lmax; j++) printf("%4d  , ",ht[i][j]);
                printf("%4d  ]",ht[i][j]);

                k = FgetIndex(lmax,mcSize,ht,ht[i]);
                printf(" %5d",k);

                k = 0;
                for (j = 0; j < 2*lmax+1; j++) k = k + ht[i][j];
                if (k != Npar)
                {
                    // Self consistency check, if the Fock state has
                    // the total number of particles supplied
                    printf("\n\nERROR: Number of particles in the Fock ");
                    printf("state above is different from total number of ");
                    printf("particles given %d\n\n",Npar);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    printf("\n\nThere are %d Fock states with L = %d",mcSize,totalL);
    printf("\nthere are %d Fock states with arbitrary L.",nc);

    printf("\nTime to setup multiconfig. constraining the momentum : ");
    if (time_direct < 0.1) printf("%.1lf(ms)",time_direct*1000);
    else printf("%.1lf(s)",time_direct);

    printf("\nTime to setup multiconfig. naively using all config. : ");
    if (time_naive < 0.1) printf("%.1lf(ms)",time_naive*1000);
    else printf(" %.1lf(s)",time_naive);

    printf("\nMemory required constraining the momentum : %.1lf(Mb)",
            ((double) mcSize*(2*lmax+1)*sizeof(int))/1E6);
    printf("\nMemory required naively using all config. : %.1lf(Mb)",
            ((double) mcSize*(2*lmax+1) + nc)*sizeof(int)/1E6);

    for (i = 0; i < mcSize; i++) free(ht[i]);
    free(ht);

    printf("\n\nDone\n\n");

    return 0;
}
