#include <time.h>
#include "LBoseFermiFockSpace.h"

int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        nc,
        Npar_B,
        Npar_F,
        lmax_B,
        lmax_F,
        totalL,
        mcSize;

    double
        l,
        time_used,
        mem_reqGB;

    clock_t
        start,
        end;

    BFCompoundSpace
        S;

    if (argc != 6)
    {
        printf("\n\nERROR: Need three integer numbers from command line");
        printf("\n\t1. Number of particles of species A");
        printf("\n\t2. Number of particles of species B");
        printf("\n\t3. max. IPS angular momentum species A");
        printf("\n\t4. max. IPS angular momentum species B");
        printf("\n\t5. Total angular momentum\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar_B);
    sscanf(argv[2],"%d",&Npar_F);
    sscanf(argv[3],"%d",&lmax_B);
    sscanf(argv[4],"%d",&lmax_F);
    sscanf(argv[5],"%d",&totalL);

    printf("\nAssembling the Bose-Fermi configurational space ");
    printf("...");

    start = clock(); // trigger to measure time
    S = BoseFermiBasis(Npar_B,Npar_F,lmax_B,lmax_F,totalL);
    end = clock();   // finish time measure
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nConfigurational space successfully set up in ");
    if (time_used > 0.1) printf("%.1lf(s).",time_used);
    else                 printf("%.0lf(ms).",time_used*1000);
    printf("\nThere are %d states with L = %d | ",S->size,totalL);
    mem_reqGB = ((double) EstimateMemoryBF(S))/1E9;
    printf("Estimated memory needed %.3lf(GB)",mem_reqGB);
    mem_reqGB = ((double) S->size * sizeof(double complex))/1E9;
    printf("\nEach Vector State in this basis cost %.3lf(GB)",mem_reqGB);

    if (S->size < 5000) PrintConfigBF(S);

    freeBoseFermiSpace(S);

    printf("\n\nDone\n\n");
    return 0;
}
