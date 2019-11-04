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
 * icc timeHomp.c -qopenmp -lm -o exe (if available)
 * gcc timeHomp.c -fopenmp -lm -o exe
 *
 * comments :
 * ----------
 *
 * On execution, call parallelized routine to apply the  hamiltonian on
 * the coefficients of a state in Fock basis, varying the configuration
 * space. The result is recorded in 'time_used_omp.dat'  file where the
 * first column is the number of particles, the second is the number of
 * orbitals and the third the average time demanded by the routine in
 * 10 runs in miliseconds. Additionally another file with the  standart
 * deviation of the runs is recorded as well
 *
 * OBS.: Time is measured using 'omp_get_wtime()' routine
 *
 * ----------------------------------------------------------------------- */

#include "hamiltonianMatrix.h"
#include "outTextFile.h"



int main(int argc, char * argv[])
{

    /*** NUMBER OF THREADS USED ***/
    omp_set_num_threads(2);

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

    double
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
        Hint;

    Cmatrix
        Ho;

    // record time and fluctuations by standart deviation
    times = fopen("time_used_omp.dat", "w");
    times_std = fopen("time_std_omp.dat", "w");


/** COLLECT TIME VARYING CONFIGURATIONAL PARAMETER, NUMBER OF  PARTICLES
  * OR NUMBER OF INDIVIDUAL PARTICLE STATES. TO SWITCH WICH PARAMETER IS
  * FIXED WITH THE ONE IS BEING VARIED, JUST SWAP THE PLACES  OF  'Npar'
  * AND 'Morb' VARIABLES IN THE NEXT TWO LINES                       **/

    Npar = 5;

    for (Morb = 4; Morb < 12; Morb++)
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

        Ho = (double complex **) malloc(Morb * sizeof(double complex *));
        for (i = 0; i < Morb; i++)
        {
            Ho[i] = (double complex *) malloc(Morb * sizeof(double complex));
        }



        C = carrDef(nc);
        out = carrDef(nc);
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

        fprintf(times,"%d  ",Npar);
        fprintf(times,"%d  ",Morb);
        fprintf(times_std,"%d  ",Npar);
        fprintf(times_std,"%d  ",Morb);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = omp_get_wtime();
            applyHconf_omp(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                    IFmat,C,Ho,Hint,out);
            end = omp_get_wtime();
            t[i] = ((double) (end - start));
            time_used += t[i];
        }
        time_used = time_used / 10;

        time_std = 0;
        for (i = 0; i < 10; i++)
        {
            time_std += (time_used - t[i]) * (time_used - t[i]);
        }
        time_std = sqrt(time_std / 10);

        printf("Time demanded : (%.3lf +- %.4lf)ms \n\n",
                time_used * 1000, time_std * 1000);

        fprintf(times, "%.6lf  ", time_used * 1000);
        fprintf(times_std, "%.6lf  ", time_std * 1000);

        free(C);
        free(out);
        free(Map);

        for(i = 0; i < Morb; i++) free(Ho[i]);
        free(Ho);

        free(Hint);

        free(IFmat);
        free(NCmat);

        free(strideOT);
        free(strideTT);
        free(MapOT);
        free(MapTT);

        fprintf(times,"\n");
        fprintf(times_std,"\n");

    }

    fclose(times);
    fclose(times_std);

    printf("\n\nDone.\n\n");
    return 0;
}