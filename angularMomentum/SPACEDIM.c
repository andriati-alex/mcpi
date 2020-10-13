#include "LBoseBoseFockSpace.h"
#include "LBoseFermiFockSpace.h"

int main(int argc, char * argv[])
{

    char
        specType;

    int
        Npar_A,
        Npar_B,
        lmax_A,
        lmax_B,
        totalL,
        mcSize;

    CompoundSpace
        bb_space;

    BFCompoundSpace
        bf_space;

    if (argc < 5)
    {
        printf("\n\nINPUT ERROR: Need at least 4 command line parameters. ");
        printf("%d given.\n\n",argc-1);
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%c",&specType);

    switch (specType)
    {
        case '1':
            sscanf(argv[2],"%d",&Npar_A);
            sscanf(argv[3],"%d",&lmax_A);
            sscanf(argv[4],"%d",&totalL);
            mcSize = BFixedMom_mcsize(Npar_A,lmax_A,totalL);
            break;
        case 'b':
            if (argc < 7)
            {
                printf("\n\nINPUT ERROR: For mixtures, need at least 6 ");
                printf("command line parameters. %d given.\n\n",argc-1);
                exit(EXIT_FAILURE);
            }
            sscanf(argv[2],"%d",&Npar_A);
            sscanf(argv[3],"%d",&lmax_A);
            sscanf(argv[4],"%d",&Npar_B);
            sscanf(argv[5],"%d",&lmax_B);
            sscanf(argv[6],"%d",&totalL);
            bb_space = AllocCompBasis(Npar_A,Npar_B,lmax_A,lmax_B,totalL);
            mcSize = bb_space->size;
            break;
        case 'f':
            if (argc < 7)
            {
                printf("\n\nINPUT ERROR: For mixtures, need at least 6 ");
                printf("command line parameters. %d given.\n\n",argc-1);
                exit(EXIT_FAILURE);
            }
            sscanf(argv[2],"%d",&Npar_A);
            sscanf(argv[3],"%d",&lmax_A);
            sscanf(argv[4],"%d",&Npar_B);
            sscanf(argv[5],"%d",&lmax_B);
            sscanf(argv[6],"%d",&totalL); 
            bf_space = BoseFermiBasis(Npar_A,Npar_B,lmax_A,lmax_B,totalL);
            mcSize = bf_space->size;
            break;
    }

    printf("\nMulticonfig. space dimension :  %d\n\n",mcSize);
    return 0;
}
