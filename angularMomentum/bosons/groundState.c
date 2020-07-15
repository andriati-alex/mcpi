#include <time.h>
#include "HEigenState.h"

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
        nnz,
        Npar, // number of particles
        lmax, // number of orbitals
        totalL,
        mcSize,
        lan_it,
        coefMemory,
        nthreads;

    double
        l,
        g,
        E0,
        sum,
        real,
        imag,
        time_used;

    Iarray
        momentum;

    Carray
        C,
        Ho;

    nthreads = omp_get_max_threads() / 2;
    omp_set_num_threads(nthreads);

    if (argc != 5)
    {
        printf("\n\nERROR: Need four command line arguments");
        printf("\n\t1. Number of particles");
        printf("\n\t2. max. IPS angular momentum");
        printf("\n\t3. Total angular momentum");
        printf("\n\t4. Contact interaction strength\n\n");
        exit(EXIT_FAILURE);
    }

    // READ COMMAND LINE ARGUMENTS TO SETUP THE  CONFIGURATIONAL
    // SPACE WITH FIXED TOTAL ANGULAR  MOMENTUM  'totalL'  GIVEN
    // NUMBER OF PARTICLES AND THE MAXIMUM  INDIVIDUAL  MOMENTUM
    // 'lmax' OF ORBITALS. FINALLY, GET THE INTERACTION STRENGTH
    // PARAMETER 'g'
    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&lmax);
    sscanf(argv[3],"%d",&totalL);
    sscanf(argv[4],"%lf",&g);

    momentum = iarrDef(totalL+1);
    for (i = 0; i <= totalL; i++)
    {
        momentum[i] = i;
    }

    mcSize = BFixedMom_mcsize(Npar,lmax,totalL);

    Ho = carrDef(2*lmax+1);
    for (i = 0; i < 2*lmax+1; i++)
    {
        l = 2 * PI * (i-lmax);
        Ho[i] = 0.5*l*l;
    }

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
        // E0 = ground(lan_it,Npar,lmax,totalL,C,Ho,g);
        groundScanning(lan_it,Npar,lmax,momentum,totalL+1,Ho,g,"test.dat");
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
            while (lan_it*coefMemory < 2*MEMORY_TOL && lan_it < 500)
            {
                lan_it += 10;
            }
        }

        time_used = omp_get_wtime(); // trigger to measure time
        // E0 = ground(lan_it,Npar,lmax,totalL,C,Ho,g);
        groundScanning(lan_it,Npar,lmax,momentum,totalL+1,Ho,g,"test.dat");
        time_used = omp_get_wtime() - time_used; // finish time measure
    }

    printf("\n\nGround state energy per particle : %.10lf ",E0/Npar);
    printf("(with %d Lanczos iterations)",lan_it);
    printf("\nTime to compute ground state (with %d threads) : ",nthreads);
    printf("%.1lf(s)",time_used);

    printf("\n\nDone.\n\n");

    free(C);
    free(Ho);

    return 0;
}
