
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of São Paulo - Brazil

 Last update : December/18/2019

---------------------------------------------------------------------------

 ****  PROGRAM TO DEMONSTRATE HOW THE MAPPING OF CONFIGURATIONS WORKS

 * Demonstrate the basics functions useful for computation of many-body
 * quantities in Fock state basis.
 *
 * COMPILATION :
 *
 * icc demonstrateFockMap.c -lm -o exe (if intel compiler is available)
 * gcc demonstrateFockMap.c -lm -o exe
 *
 * EXECUTION :
 *
 * ./exe Nparticles Norbitals
 *
 * First command line argument is the number of particles and  the second
 * is the number of orbitals. It prints on the screen all the Fock states
 * enumerated, that is, the hashing table. Moreover it shows the mappings
 * for some cases
 *
 * IMPORTANT NOTE : this is suppose to be a simple demonstration, then do
 * not use large number of particles or orbitals because it would mess up
 * the output on the screen
 *
--------------------------------------------------------------------------- **/



#include "configurationsMap.h"

int HaveMomentum(unsigned int Morb, int * conf, int lz)
{

/** Verify if the total momentum of the configuration is equal to 'lz'.
  * Here 'lz' denotes the momentum associate quantum number, which for
  * the physical problem can be either the angular(ring) or linear mom.
    Ascending order is assumed for the momentum quantum number occupat
    lz = [ -(Morb-1)/2 , ..., 0, ... (Morb-1)/2 ]     for  'Morb'  odd
    lz = [ -(Morb/2) , ..., 0, ... Morb/2-1 ]         for  'Morb' even **/

    int
        i,
        totalMom;

    totalMom = 0;
    for (i = 0; i < Morb; i++) totalMom = totalMom + (i - Morb/2)*conf[i];

    if (totalMom == lz) return 1;
    return 0;
}



int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        lz,
        nc,
        Npar, // number of particles
        Morb, // number of orbitals
        mustPrint;

    Iarray
        NCmat,
        IFmat;

    if (argc != 4)
    {
        printf("\n\nERROR: Need three integer numbers from command line ");
        printf("first the number of particles, second the number of states");
        printf("and third the total momentum quantum number.\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    sscanf(argv[3],"%d",&lz);
    nc = NC(Npar,Morb);

    printf("\n\nNumber of particles : %3d", Npar);
    printf(  "\nNumber of orbitals  : %3d", Morb);
    printf(  "\nNumber of configurations : %d", nc);

    if (nc > 5000 || Morb > 7)
    {
        printf("\n\nWARNING: lARGE SYSTEM CAN MESS UP THE OUTPUT\n\n");
    }

    NCmat = setupNCmat(Npar,Morb);
    IFmat = setupFocks(Npar,Morb);

    printf("\n\n\n");

    if (nc < 10000)
    {

        printf("Configuration  [");

        for (i = 0; i < Morb - 1; i++) printf(" Orb%d ,",i + 1);
        printf(" Orb%d ]", Morb);

        printf("\n=============================================");
        printf("================================");

        for (i = 0; i < nc; i++)
        {
            mustPrint = HaveMomentum(Morb,&IFmat[Morb*i],lz);
            if (mustPrint)
            {
                printf("\n%8d       [", i);
                for (j = 0; j < Morb - 1; j++)
                {
                    printf(" %3d  ,",IFmat[j+Morb*i]);
                }
                printf(" %3d  ]", IFmat[Morb - 1 + Morb*i]);

                k = FockToIndex(Npar,Morb,NCmat,&IFmat[Morb*i]);
                if (k != i)
                {
                    // Self consistency check, if the Fock state in the
                    // hashing table that was constructed using IndexToFock
                    // function gives the correc index when converted back
                    printf("\n\nERROR: Wrong map from FockToIndex\n\n");
                }
            }
        }
    }

    else
    {
        printf("\n\nWARNING : Not printing, too large system with ");
        printf("total number of configurations higher than 10000\n");
    }



    free(IFmat);
    free(NCmat);

    printf("\n\nDone.\n\n");
    return 0;
}
