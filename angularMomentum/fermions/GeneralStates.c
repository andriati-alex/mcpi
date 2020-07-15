
/*****  AUTHOR INFORMATION
 
NAME : Alex Valerio Andriati
AFFILIATION : University of São Paulo - Brazil

Last update : July/01/2020



 *****  PROGRAM TO DEMONSTRATE HOW THE MAPPING OF CONFIGURATIONS WORKS

Demonstrate the basics functions useful for computation of many-body
quantities in Fock state basis.

**  COMPILATION :

icc demonstrateFockMap.c -o exe (if intel compiler is available)
gcc demonstrateFockMap.c -o exe

**  EXECUTION :

./exe Nparticles Norbitals

First command line argument is the number of particles and  the second
is the number of orbitals.  It prints on the screen all configurations
enumerated, that is, the hashing table. For each printed configuration
it perform two self-consistency checks using the functions that relate
the configurations and their index.

IMPORTANT NOTE : this is suppose to be a simple demonstration, then do
not use large number of particles or orbitals because it would mess up
the output on the screen

*****/



#include "FermiFockSpace.h"

#define PRINT_TOL 10000



int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        s,
        q,
        l,
        g,
        h,
        nc,
        Npar, // number of particles
        Morb; // number of orbitals

    Iarray
        IFmat;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("first the number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);

    if (Morb < Npar)
    {
        printf("\n\nERROR: The number of fermions (%d) cannot exceed ",Npar);
        printf("the number of orbitals to allocate them (%d) due to ",Morb);
        printf("Pauli exclusion principle\n\n");
        exit(EXIT_FAILURE);
    }

    nc = NC(Npar,Morb);

    printf("\n\nNumber of particles : %3d", Npar);
    printf(  "\nNumber of orbitals  : %3d", Morb);
    printf(  "\nNumber of configurations : %d", nc);

    IFmat = setupFocks(Npar,Morb);
    printf("\n\n\n");

    if (nc < PRINT_TOL)
    {

        if (Morb > 7)
        {
            printf("\n\nWARNING: lARGE SYSTEM CAN MESS UP THE OUTPUT, ");
            printf("COMPACT PRINTING\n\n");

            printf(" Configuration    Occupations as binary");
            for (i = 0; i < nc; i++)
            {
                printf("\n%9d         [ ", i);

                for (j = 0; j < Morb; j++) printf("%d",IFmat[j+i*Morb]);
                printf(" ]");

                k = FockToIndex(Npar,Morb,&IFmat[Morb*i]);
                if (k != i)
                {
                    // Self consistency check,  if  the  Fock  state in the
                    // hashing table that was constructed using IndexToFock
                    // function gives the correc index when converted  back
                    printf("\n\nERROR: Wrong map from FockToIndex\n\n");
                    exit(EXIT_FAILURE);
                }

                k = 0;
                for (j = 0; j < Morb; j++) k = k + IFmat[j+i*Morb];
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
            free(IFmat);
            printf("\n\nDone.\n\n");
            return 0;
        }

        printf("Configuration  [");

        for (i = 0; i < Morb - 1; i++) printf(" Orb%d ,",i + 1);
        printf(" Orb%d ]\n", Morb);

        printf("=============================================");
        printf("================================");

        for (i = 0; i < nc; i++)
        {

            printf("\n%8d       [", i);

            for (j = 0; j < Morb - 1; j++) printf(" %3d  ,",IFmat[j+i*Morb]);
            printf(" %3d  ]", IFmat[Morb - 1 + i*Morb]);

            k = FockToIndex(Npar,Morb,&IFmat[Morb*i]);
            if (k != i)
            {
                // Self consistency check, if the Fock state in the
                // hashing table that was constructed using IndexToFock
                // function gives the correc index when converted back
                printf("\n\nERROR: Wrong map from FockToIndex\n\n");
            }

            k = 0;
            for (j = 0; j < Morb; j++) k = k + IFmat[j+i*Morb];
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
        printf("\n\nWARNING : Not printing, conf. space too large with ");
        printf("total number of configurations larger than %d\n",PRINT_TOL);
    }



    free(IFmat);

    printf("\n\nDone.\n\n");
    return 0;
}
