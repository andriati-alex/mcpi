#include "DensityCorrelation.h"

int main(int argc, char * argv[])
{

    if (argv[1][0] == '1') SCANNING_ANALYSIS(argv[2]);
    if (argv[1][0] == 'b') SCANNING_MIXTURE_ANALYSIS(argv[2]);
    if (argv[1][0] == 'f') SCANNING_BOSEFERMI_ANALYSIS(argv[2]);
    printf("\n\nDone !\n\n");

    return 0;
}
