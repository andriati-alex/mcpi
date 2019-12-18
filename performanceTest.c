
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : November/02/2019

 -------------------------------------------------------------------------

 ****  TEST ROUTINES TO COMPUTE PHYSICAL OPERATORS
 *
 * COMPILE :
 *
 * icc performanceTest.c -lm -qopenmp -o exe (if intel compiler is available)
 * gcc performanceText.c -lm -fopenmp -o exe
 *
 * HOW TO EXECUTE :
 *
 * ./exe Nparticles Norbitals
 *
 * where Nparticles and Morbitals are command line arguments for the
 * number of particles and individual particle states respectively.
 * Execute each routine to compute density matrices and apply  the
 * Hamiltonian using different improvements and display the result
 * on screen.
 *
 * ----------------------------------------------------------------------- */

#include "onebodyMatrix.h"
#include "twobodyMatrix.h"
#include "hamiltonianMatrix.h"
#include "outTextFile.h"
#include <time.h>



int main(int argc, char * argv[])
{

    /*** NUMBER OF THREADS USED - CHOOSE ACCORDINGLY TO CPU LIMITS ***/
    omp_set_num_threads(2);

    int
        i,
        j,
        q,
        l,
        nc,
        Npar,
        Morb,
        nthreads;

    double
        sum,
        end_omp,
        start_omp,
        time_used;

    double complex
        z;

    clock_t
        start,
        end;

    Iarray
        Map,
        MapOT,
        MapTT,
        IFmat,
        NCmat,
        strideOT,
        strideTT;

    Carray
        C,
        out,
        out_X,
        out_XX,
        rho2,
        rho2_X,
        rho2_XM,
        Hint;

    Cmatrix
        rho1_XM,
        rho1_X,
        Ho;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    nthreads = omp_get_max_threads();

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    nc = NC(Npar,Morb);

    NCmat = setupNCmat(Npar,Morb);
    IFmat = setupFocks(Npar,Morb);

    printf("\nNumber of particles : %3d", Npar);
    printf("\nNumber of orbitals  : %3d", Morb);
    printf("\nNumber of configurations : %d", nc);

    printf("\n\n======================================\n\n");

    printf("MEMORY CONSUMPTION (in Mb)");

    printf("\n\nMemory for coefficients : %.1lf",
            ((double) nc*sizeof(double complex)) / 1E6);

    printf("\nMemory for Hashing table : %.1lf",
            ((double) nc*Morb*sizeof(int)) / 1E6);

    printf("\nMemory for single jump from one orbital(1J1O) : %.1lf",
            ((double) nc*Morb*Morb*sizeof(int)) / 1E6);



    strideTT = iarrDef(nc);
    strideOT = iarrDef(nc);
    Map = OneOneMap(Npar,Morb,NCmat,IFmat);
    MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
    MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

    printf("\nMemory for double jump from one orbital(2J1O) : %.1lf",
            ((double) strideOT[nc-1]*sizeof(int))/1E6);
    printf("\nMemory for double jump from two orbitals(2J2O) : %.1lf",
            ((double) strideTT[nc-1]*sizeof(int))/1E6);




    rho1_X = (double complex **) malloc(Morb * sizeof(double complex *));
    for (i = 0; i < Morb; i++)
    {
        rho1_X[i] = (double complex *) malloc(Morb * sizeof(double complex));
    }

    rho1_XM = (double complex **) malloc(Morb * sizeof(double complex *));
    for (i = 0; i < Morb; i++)
    {
        rho1_XM[i] = (double complex *) malloc(Morb * sizeof(double complex));
    }

    Ho = (double complex **) malloc(Morb * sizeof(double complex *));
    for (i = 0; i < Morb; i++)
    {
        Ho[i] = (double complex *) malloc(Morb * sizeof(double complex));
    }



    C = carrDef(nc);
    out = carrDef(nc);
    out_X = carrDef(nc);
    out_XX = carrDef(nc);
    rho2 = carrDef(Morb * Morb * Morb * Morb);
    rho2_X = carrDef(Morb * Morb * Morb * Morb);
    rho2_XM = carrDef(Morb * Morb * Morb * Morb);
    Hint = carrDef(Morb*Morb*Morb*Morb);

    sum = 0.0;
    for (i = 0; i < nc; i++)
    {
        C[i] = sin(20*((double)i)/nc) * (i%13) + ((i%8)-(i%3)) * I;
        sum = sum + creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
    }

    // normalize to 1
    for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);



    for (i = 0; i < Morb; i++)
    {
        Ho[i][i] = (i % 4) - 1;
        for (j = i + 1; j < Morb; j++)
        {
            Ho[i][j] = i*(j%3) - (i%4) + 5*(j%2) - I*(4.123*i/(j + 1));
            Ho[j][i] = conj(Ho[i][j]);
        }
    }


    for (i = 0; i < Morb*Morb*Morb*Morb; i++) Hint[i] = 1.234;

    for (i = 0; i < Morb; i++)
    {
        for (j = i + 1; j < Morb; j++)
        {
            for (q = 0; q < Morb; q++)
            {
                if (q == i || q == j) continue;
                for (l = q + 1; l < Morb; l++)
                {
                    if (l == i || l == j) continue;
                    // real part
                    z = i - 2 + 10 * (j % (i+1)) - q * l;
                    // imag part
                    z = z + I * ((double) i * q - j * l) / Morb;
                    Hint[i+j*Morb+q*Morb*Morb+l*Morb*Morb*Morb] = z;
                    Hint[i+j*Morb+l*Morb*Morb+q*Morb*Morb*Morb] = z;
                    Hint[j+i*Morb+l*Morb*Morb+q*Morb*Morb*Morb] = z;
                    Hint[j+i*Morb+q*Morb*Morb+l*Morb*Morb*Morb] = z;
                    Hint[q+l*Morb+i*Morb*Morb+j*Morb*Morb*Morb] = conj(z);
                    Hint[q+l*Morb+j*Morb*Morb+i*Morb*Morb*Morb] = conj(z);
                    Hint[l+q*Morb+i*Morb*Morb+j*Morb*Morb*Morb] = conj(z);
                    Hint[l+q*Morb+j*Morb*Morb+i*Morb*Morb*Morb] = conj(z);
                }
            }
        }
    }



    printf("\n\n======================================\n\n");

    printf("TIME DEMANDED");



    start = clock();
    for (i = 0; i < 50; i++) OBrho_X(Npar,Morb,NCmat,IFmat,C,rho1_X);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 50;
    printf("\n\nTime to setup rho1 with hashing table");
    printf(" : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 50; i++) OBrho_XM(Npar,Morb,Map,NCmat,IFmat,C,rho1_XM);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 50;
    printf("\n\nTime to setup rho1 using 1J1O mapping");
    printf(" : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++) TBrho(Npar,Morb,NCmat,IFmat,C,rho2);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to setup rho2 with hashing table");
    printf(" : %.3lfms",time_used*1000);

    start = clock();
    for (i = 0; i < 5; i++)
    {
        TBrho_X(Npar,Morb,Map,MapOT,strideOT,NCmat,IFmat,C,rho2_X);
    }
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to setup rho2 with 1J1O/2J1O mapping : ");
    printf("%.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++)
    {
        TBrho_XX(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                NCmat,IFmat,C,rho2_XM);
    }
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to setup rho2 with all mappings : ");
    printf("%.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++) applyHconf(Npar,Morb,NCmat,IFmat,C,Ho,Hint,out);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to apply H with Hashing Table");
    printf(" : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++)
    {
        applyHconf_X(Npar,Morb,Map,MapOT,strideOT,NCmat,IFmat,C,Ho,Hint,out_X);
    }
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to apply H with 1J1O/2J1O mapping : ");
    printf("%.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++)
    {
        applyHconf_XX(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                      IFmat,C,Ho,Hint,out_XX);
    }
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to apply H with all mappings : ");
    printf("%.3lfms", time_used * 1000);


    // Parallelized routine
    time_used = 0;
    for (i = 0; i < 5; i++)
    {
        start_omp = omp_get_wtime();
        applyHconf_omp(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                       IFmat,C,Ho,Hint,out);
        end_omp = omp_get_wtime();
        time_used = time_used + ((double) (end_omp - start_omp));
    }
    time_used = time_used / 5;
    printf("\n\nTime to apply H with all mappings(%d threads) : ",nthreads);
    printf("%.3lfms", time_used * 1000);



    free(C);
    free(out);
    free(out_X);
    free(out_XX);
    free(Map);

    for(i = 0; i < Morb; i++) free(rho1_X[i]);
    free(rho1_X);

    for(i = 0; i < Morb; i++) free(rho1_XM[i]);
    free(rho1_XM);

    for(i = 0; i < Morb; i++) free(Ho[i]);
    free(Ho);

    free(Hint);

    free(IFmat);
    free(NCmat);

    free(strideOT);
    free(strideTT);
    free(MapOT);
    free(MapTT);
    free(rho2);
    free(rho2_X);

    printf("\n\nDone.\n\n");
    return 0;
}
