#include <time.h>
#include "GroundStateDrivers.h"

#define PI 3.141592653589793

/* icc -o exe groundStateLieb.c \
            -L${MKLROOT}/lib/intel64 \
            -lmkl_intel_lp64 \
            -lmkl_gnu_thread \
            -lmkl_core \
            -qopenmp \
            -O3 \
*/



int main(int argc, char * argv[])
{

    int
        i,
        Npar_A, // number of particles
        Npar_B, // number of particles
        lmax_A, // number of orbitals
        lmax_B, // number of orbitals
        totalL,
        lan_it,
        mcSize,
        coefMemory,
        nthreads;

    double
        l,
        E0,
        sum,
        real,
        imag,
        time_used;

    double
        g[3];

    Carray
        C,
        HoA,
        HoB;

    CompoundSpace
        MixSpace;

    nthreads = omp_get_max_threads() / 2;
    omp_set_num_threads(nthreads);

    if (argc != 7)
    {
        printf("\n\nERROR: Need three integer numbers from command line");
        printf("\n\t1. Number of particles of species A");
        printf("\n\t2. Number of particles of species B");
        printf("\n\t3. max. IPS angular momentum species A");
        printf("\n\t4. max. IPS angular momentum species B");
        printf("\n\t5. Total angular momentum");
        printf("\n\t6. Interaction strength\n\n");
        exit(EXIT_FAILURE);
    }

    // READ COMMAND LINE ARGUMENTS TO SETUP THE  CONFIGURATIONAL
    // SPACE WITH FIXED TOTAL ANGULAR  MOMENTUM  'totalL'  GIVEN
    // NUMBER OF PARTICLES AND THE MAXIMUM  INDIVIDUAL  MOMENTUM
    // 'lmax' OF ORBITALS. FINALLY, GET THE INTERACTION STRENGTH
    // PARAMETER 'g'
    sscanf(argv[1],"%d",&Npar_A);
    sscanf(argv[2],"%d",&Npar_B);
    sscanf(argv[3],"%d",&lmax_A);
    sscanf(argv[4],"%d",&lmax_B);
    sscanf(argv[5],"%d",&totalL);
    sscanf(argv[6],"%lf",&g[0]);

    printf("\nAssembling the configurational space ");
    printf("...");

    MixSpace = AllocCompBasis(Npar_A,Npar_B,lmax_A,lmax_B,totalL);
    printf("\nConfigurational space successfully set up");
    mcSize = MixSpace->size;

    HoA = carrDef(2*lmax_A+1);
    for (i = 0; i < 2*lmax_A+1; i++)
    {
        l = 2 * PI * (i-lmax_A);
        HoA[i] = 0.5*l*l;
    }

    HoB = carrDef(2*lmax_B+1);
    for (i = 0; i < 2*lmax_B+1; i++)
    {
        l = 2 * PI * (i-lmax_B);
        HoB[i] = 0.5*l*l;
    }

    g[1] = g[0];
    g[2] = 0;

    // define an initial guess for Lanczos tridiagonal decomposition
    C = carrDef(mcSize);
    sum = 0;
    for (i = 0; i < mcSize; i++)
    {
        real = (i % 3) - 1 + 4 * (i % 5);
        imag = 3 - (i % 7) + 2 * (i % 3);
        C[i] = real + I * imag;
        sum = sum + real*real + imag*imag;
    }
    // normalize to 1
    for (i = 0; i < mcSize; i++) C[i] = C[i] / sqrt(sum);

    if (mcSize < 201)
    {
        lan_it = mcSize - 1;
        time_used = omp_get_wtime(); // trigger to measure time
        E0 =  BOSEBOSE_GS(lan_it,MixSpace,C,HoA,HoB,g);
        time_used = omp_get_wtime() - time_used; // finish time measure
    }
    else
    {
        // Perform between 200 and 500 Lanczos  iterations
        // depending on the config. space size.  For  very
        // large spaces restart iterations to moderate the
        // memory usage.
        coefMemory = mcSize * sizeof(double complex);
        if (mcSize < 1000) lan_it = 200;
        else
        {
            lan_it = 10;
            while (lan_it*coefMemory < 2*MEMORY_TOL && lan_it < 300)
            {
                lan_it += 10;
            }
        }

        time_used = omp_get_wtime(); // trigger to measure time
        E0 =  BOSEBOSE_GS(lan_it,MixSpace,C,HoA,HoB,g);
        time_used = omp_get_wtime() - time_used; // finish time measure
    }

    printf("\n\nGround state energy per particle : %.10lf ",E0/(Npar_A+Npar_B));
    printf("(with %d Lanczos iterations)",lan_it);
    printf("\nTime to compute ground state (with %d threads) : ",nthreads);
    printf("%.1lf(s)",time_used);

    printf("\n\nDone.\n\n");

    free(C);
    free(HoA);
    free(HoB);
    freeCompSpace(MixSpace);

    return 0;
}
