#include "confMap.h"



int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        nc,
        Npar, // number of particles
        Morb; // number of orbitals

    Iarray
        NCmat,
        IFmat,
        Map;

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

    if (Npar > 99 || Morb > 9)
    {
        printf("\n\nWARNING: Bad output print format.");
        printf(" Number of particles or orbitals too large\n\n");
    }

    NCmat = MountNCmat(Npar,Morb);
    IFmat = MountFocks(Npar,Morb);

    Map = JumpMapping(Npar,Morb,NCmat,IFmat);



    printf("\n\nNumber of particles : %3d", Npar);
    printf(  "\nNumber of orbitals  : %3d", Morb);
    printf(  "\nNumber of configurations : %d", nc);

    printf("\n\n=============================================\n\n");



    printf("Configuration Number  |  [ occupations ]  |  FockToIndex ");
    printf(" | Transitions\n");
    for (i = 0; i < nc; i++)
    {
        printf("\n%8d      [", i);

        for (j = 0; j < Morb - 1; j++) printf(" %3d |", IFmat[j + Morb*i]);
        printf(" %3d ]", IFmat[Morb - 1 + Morb*i]);

        k = FockToIndex(Npar,Morb,NCmat,&IFmat[Morb*i]);
        if (k != i)
        {
            printf("\n\nERROR: Wrong map from FockToIndex\n\n");
        }

        printf("      %d", k);

        printf("    %d -> %d", 1, 3);
        printf(" goes to index %d",Map[i+1*nc+3*nc*Morb]);
    }

    free(IFmat);
    free(NCmat);
    free(Map);

    printf("\n\n");
    return 0;
}
