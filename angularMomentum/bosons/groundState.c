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
        j,
        k,
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

    clock_t
        start,
        end;

    Iarray
        NNZrow,
        * ht;

    Carray
        C,
        Ho;

    HConfMat
        H;

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



    // CONFIGURE MULTICONFIGURATIONAL SPACE - HASHING TABLE
    printf("\nConfiguring multiconf. space structures ...\n");

    start = clock(); // trigger to measure time
    mcSize = BFixedMom_mcsize(Npar,lmax,totalL);
    ht = BAssembleHT(Npar,lmax,totalL,mcSize);
    end = clock();   // stop time measure
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("\nConfigurational basis successfully setup, ");
    if (time_used < 0.1)
    {
        printf("Time to setup Conf. space : %.1lf(ms)",time_used*1000);
    }
    else
    {
        printf("Time to setup Conf. space : %.1lf(s)",time_used);
    }

    printf("\n\nConfiguring Hamiltonian Matrix ...\n");

    // Compute the total Number of NonZero elements in Hamiltonian matrix
    // 'nnz' and also record the Number of NonZero per row in NNZrow
    NNZrow = iarrDef(mcSize);
    nnz = NNZ_PerRow(Npar,lmax,mcSize,ht,NNZrow);

    k = 0; // maximum NNZ in a same row of the matrix
    j = 0; // Self-Consistency check with the returned value 'nnz'
    for (i = 0; i < mcSize; i++)
    {
        if (k < NNZrow[i]) k = NNZrow[i];
        j = j + NNZrow[i];
    }
    if (j != nnz)
    {
        printf("\n\nFATAL ERROR : FUNCTION TO COMPUTE ");
        printf("NONZERO ELEMENTS OF SPARSE MATRIX IS WRONG\n\n");
        exit(EXIT_FAILURE);
    }



    // ASSEMBLE THE HAMILTONIAN (SPARSE) MATRIX STRUCTURE



    // Configure the single-particle hamiltonian matrix.
    // Since it is diagonal  it is store in a  1D array.
    Ho = carrDef(2*lmax+1);
    for (i = 0; i < 2*lmax+1; i++)
    {
        l = 2 * PI * (i-lmax);
        Ho[i] = 0.5*l*l;
    }

    start = clock(); // trigger to measure time
    H = assembleH(Npar,lmax,mcSize,ht,Ho,g);
    end = clock();   // stop time measure
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    if (time_used < 0.1)
    {
        printf("\nTime to setup multiconf. Hamiltonian matrix : ");
        printf("%.1lf(ms)",time_used*1000);
    }
    else
    {
        printf("\nTime to setup multiconf. Hamiltonian matrix : ");
        printf("%.1lf(s)",time_used);
    }



    // DISPLAY INFO ABOUT THE CONFIGURATIONAL SPACE ASSEMBLED



    printf("\n\n----------------------------\n");
    printf("Multiconf. space information\n");
    printf("----------------------------\n");
    printf("Total number of states (with L = %2d subjected to ",totalL);
    printf("lmax = %d) is %d",lmax,mcSize);
    printf("\nNumber of nonzero entries in H : %d ",nnz);
    printf("Sparsity = %.5lf",1.0 - ((double) nnz)/mcSize/mcSize);
    printf("\nMaximum number of nonzero entries in a same row : %d",k);
    printf("\n=====================================");
    printf("=====================================\n");



    // COMPUTE THE GROUND STATE ENERGY AND STATE



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
        E0 = ground(lan_it,mcSize,H,C);
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
            while (lan_it*coefMemory < 2E9 && lan_it < 500)
            {
                lan_it += 10;
            }
        }

        time_used = omp_get_wtime(); // trigger to measure time
        E0 = ground(lan_it,mcSize,H,C);
        time_used = omp_get_wtime() - time_used; // finish time measure
    }

    printf("\n\nGround state energy per particle : %.10lf ",E0/Npar);
    printf("(with %d Lanczos iterations)",lan_it);
    printf("\nTime to compute ground state (with %d threads) : ",nthreads);
    printf("%.1lf(s)",time_used);

    printf("\n\nDone.\n\n");

    for (i = 0; i < mcSize; i++) free(ht[i]);
    free(C);
    free(ht);
    free(Ho);
    free(NNZrow);
    freeHmat(H);

    return 0;
}
