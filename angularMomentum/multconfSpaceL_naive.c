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
        nL,
        Npar, // number of particles
        Morb, // number of orbitals
        mustPrint;

    Iarray
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
        printf("\n\nWARNING: LARGE SYSTEM PRODUCES BAD OUTPUT FORMAT\n\n");
    }

    IFmat = setupFocks(Npar,Morb);

    printf("\n\n\n");

    nL = 0;

    if (nc < 7500)
    {

        printf("Configurations with L = %d\n\n",lz);

        printf("Momentum IPS   [");

        for (i = 0; i < Morb - 1; i++) printf(" %3d  ,",i-Morb/2);
        printf(" %3d  ]",(Morb-1)-Morb/2);

        printf("\n=============================================");
        printf("================================");

        for (i = 0; i < nc; i++)
        {
            mustPrint = HaveMomentum(Morb,&IFmat[Morb*i],lz);
            if (mustPrint)
            {
                nL = nL + 1;
                printf("\n%8d       [", i);
                for (j = 0; j < Morb - 1; j++)
                {
                    printf(" %3d  ,",IFmat[j+Morb*i]);
                }
                printf(" %3d  ]", IFmat[Morb - 1 + Morb*i]);
            }
        }
    }

    else
    {
        printf("\n\nWARNING : Not printing, too large system or number ");
        printf("of configurations\n\n");
        for (i = 0; i < nc; i++)
        {
            mustPrint = HaveMomentum(Morb,&IFmat[Morb*i],lz);
            if (mustPrint) nL = nL + 1;
        }
    }



    free(IFmat);
    free(NCmat);

    printf("\n\nTotal number of states with L = %d : %d\n\n",lz,nL);

    printf("\n\nDone.\n\n");
    return 0;
}
