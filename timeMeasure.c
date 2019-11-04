/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of São Paulo - Brazil

 Last update : November/02/2019

---------------------------------------------------------------------------

 ****  EXECUTABLE TO MEASURE TIME FOR DIFFENRENT NUMBER OF PARTICLES ****
 *
 * compilation :
 * -------------
 *
 * icc timeMeasure.c -lm -qopenmp -o exe (if available)
 * gcc timeMeasure.c -lm -fopenmp -o exe
 *
 * comments :
 * ----------
 *
 * On execution, call routines to compute density matrices and apply the
 * hamiltonian matrix in the Fock space, exploring various improvements,
 * varying the size of the configurational space through the  number of
 * particles or number of individual particle states. The time demanded
 * by the different implementations are recorded 'time_used.dat'  file,
 * measured im miliseconds. The output file data is organized as follows
 *
 * column 1 : Number of particles
 * column 2 : Number of orbitals
 * ---------------------------------- TIME ------------------------------
 * column 3 : one-body density matrix both routines
 * column 4 : one-body density matrix hashing table
 * column 5 : one-body density matrix mappings
 * column 6 : two-body density matrix hashing table
 * column 7 : two-body density matrix jumps from one orbital mappings
 * column 8 : two-body density matrix jumps from two orbitals mappings
 * column 9 : Apply Hamiltonian with hashing table
 * column 10 : Apply Hamiltonian with jumps from one orbital mappings
 * column 11 : Apply Hamiltonian with jumps from two orbitals mappings
 *
 * The time for each of these routines is the average of 10 runs and
 * additionaly a file with the standart deviation is supplied
 *
 * OBS.: Time is measured using clock() function from time.h library
 *
 * ----------------------------------------------------------------------- */

#include "onebodyMatrix.h"
#include "twobodyMatrix.h"
#include "hamiltonianMatrix.h"
#include "outTextFile.h"
#include <time.h>



int main(int argc, char * argv[])
{

    int
        i,
        j,
        q,
        l,
        nc,
        Npar,
        Morb;

    double
        t[10],
        sum,
        time_used,
        time_std;

    FILE
        * times,
        * times_std;

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
        rho1,
        Ho;

    // record time and fluctuations by standart deviation
    times = fopen("time_used.dat", "w");
    times_std = fopen("time_std.dat", "w");

/** COLLECT TIME VARYING CONFIGURATIONAL PARAMETER, NUMBER OF  PARTICLES
  * OR NUMBER OF INDIVIDUAL PARTICLE STATES. TO SWITCH WICH PARAMETER IS
  * FIXED WITH THE ONE IS BEING VARIED, JUST SWAP THE PLACES  OF  'Npar'
  * AND 'Morb' VARIABLES IN THE NEXT TWO LINES                       **/

    Morb = 3;

    for (Npar = 30; Npar < 105; Npar += 50)
    {

        printf("\n\n");

        nc = NC(Npar,Morb);

        NCmat = setupNCmat(Npar,Morb);
        IFmat = setupFocks(Npar,Morb);

        printf("\nNumber of particles : %3d", Npar);
        printf("\nNumber of orbitals  : %3d", Morb);
        printf("\nNumber of configurations : %d", nc);

        printf("\n\n======================================\n\n");

        strideTT = iarrDef(nc);
        strideOT = iarrDef(nc);
        Map = OneOneMap(Npar,Morb,NCmat,IFmat);
        MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
        MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

        rho1 = (double complex **) malloc(Morb * sizeof(double complex *));
        for (i = 0; i < Morb; i++)
        {
            rho1[i] = (double complex *) malloc(Morb * sizeof(double complex));
        }

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

        printf("TIME DEMANDED");

        fprintf(times,"%d  ",Npar);
        fprintf(times,"%d  ",Morb);
        fprintf(times_std,"%d  ",Npar);
        fprintf(times_std,"%d  ",Morb);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            OBrho(Npar,Morb,NCmat,C,rho1);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to setup rho1 : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            OBrho_X(Npar,Morb,NCmat,IFmat,C,rho1_X);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to setup rho1 with Fock states");
        printf(" : (%.3lf +- %.4lf)ms", time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            OBrho_XM(Npar,Morb,Map,NCmat,IFmat,C,rho1_XM);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to setup rho1 using Maps and Fock states");
        printf(" : (%.3lf +- %.4lf)ms", time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            TBrho(Npar,Morb,NCmat,IFmat,C,rho2);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to setup rho2 : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            TBrho_X(Npar,Morb,Map,MapOT,strideOT,NCmat,IFmat,C,rho2_X);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to setup rho2 with one-map : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);


        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            TBrho_XX(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                    NCmat,IFmat,C,rho2_XM);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to setup rho2 with two-map : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);


        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            applyHconf(Npar,Morb,NCmat,IFmat,C,Ho,Hint,out);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to apply H with no map : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);


        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            applyHconf_X(Npar,Morb,Map,MapOT,strideOT,NCmat,IFmat,C,Ho,Hint,out_X);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to apply H with one-map : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = clock();
            applyHconf_XX(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                    IFmat,C,Ho,Hint,out_XX);
            end = clock();
            t[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("\n\nTime to apply H with two-map : (%.3lf +- %.4lf)ms",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        free(C);
        free(out);
        free(out_X);
        free(out_XX);
        free(Map);

        for(i = 0; i < Morb; i++) free(rho1[i]);
        free(rho1);

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
        free(rho2_XM);

        fprintf(times,"\n");
        fprintf(times_std,"\n");

    }

    fclose(times);
    fclose(times_std);

    printf("\n\nDone.\n\n");
    return 0;
}