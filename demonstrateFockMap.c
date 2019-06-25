#include "confMap.h"



int main(int argc, char * argv[])
{

    int
        i,
        j,
        nc,
        Npar, // number of particles
        Morb; // number of orbitals

    Iarray
        occ,
        NCmat,
        IFmat;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    nc = NC(Npar,Morb);

    if (Npar > 99 || Morb > 99)
    {
        printf("\n\nWARNING: Bad output print format.");
        printf(" Number of particles or orbitals too large\n\n");
    }

    NCmat = MountNCmat(Npar,Morb);
    IFmat = MountFocks(Npar,Morb);
    occ = iarrDef(Morb);



    printf("\n\nNumber of particles : %3d", Npar);
    printf(  "\nNumber of orbitals  : %3d", Morb);
    printf(  "\nNumber of configurations : %d", nc);

    printf("\n\n=============================================\n\n");



    printf("All configurations :\n");
    for (i = 0; i < nc; i++)
    {
        printf("\n%8d  [", i);

        for (j = 0; j < Morb - 1; j++) printf(" %3d |", IFmat[i + nc*j]);
        printf(" %3d ]", IFmat[i + nc*(Morb-1)]);

        for (j = 0; j < Morb; j++) occ[j] = IFmat[i + nc*j];

        printf(" %d", FockToIndex(Npar,Morb,NCmat,occ));
    }

    free(IFmat);
    free(NCmat);
    free(occ);

    printf("\n\n");
    return 0;
}
