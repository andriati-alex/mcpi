/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of São Paulo - Brazil

 Last update : 08/13/2019

---------------------------------------------------------------------------

 ****  EXECUTABLE TO MEASURE TIME FOR DIFFENRENT NUMBER OF PARTICLES ****
 *
 * compilation :
 * -------------
 *
 * icc timeMeasure.c -lm -o exe (if available)
 * gcc timeMeasure.c -lm -o exe
 *
 * comments :
 * ----------
 *
 *  This executable record time demanded in each routine varying the number
 *  of particles in the system. It can be also used to sweep the number of
 *  single particle states with simple changes in the first for loop by
 *  allowing 'Morb' to vary and set Npar fixed.
 *
 * ----------------------------------------------------------------------- */

#include "hamiltonianMatrix.h"
#include "outTextFile.h"



int main(int argc, char * argv[])
{

    omp_set_num_threads(8);

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
    times = fopen("time_used_omp_orb.dat", "w");
    times_std = fopen("time_std_omp_orb.dat", "w");

    Npar = 5;

    for (Morb = 4; Morb < 26; Morb++)
    {

        printf("\n\n");

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

        printf("\nMemory for Fock states : %.1lf",
                ((double) nc*Morb*sizeof(int)) / 1E6);

        printf("\nMemory for one to one Map : %.1lf",
                ((double) nc*Morb*Morb*sizeof(int)) / 1E6);



        strideTT = iarrDef(nc);
        strideOT = iarrDef(nc);
        Map = OneOneMap(Npar,Morb,NCmat,IFmat);
        MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
        MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

        printf("\nMemory for one-two Map : %.1lf",
                ((double) strideOT[nc-1]*sizeof(int))/1E6);
        printf("\nMemory for two-two Map : %.1lf",
                ((double) strideTT[nc-1]*sizeof(int))/1E6);

        printf("\nMemory for two-body matrix elements : %.1lf",
                ((double) Morb*Morb*Morb*Morb*sizeof(double complex))/1E6);



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



        printf("\n\n======================================\n\n");

        printf("TIME DEMANDED");

        fprintf(times,"%d  ",Npar);
        fprintf(times,"%d  ",Morb);
        fprintf(times_std,"%d  ",Npar);
        fprintf(times_std,"%d  ",Morb);

        time_used = 0;
        for (i = 0; i < 10; i++)
        {
            start = omp_get_wtime();
            applyHconf_XX(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
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

        printf("\n\nTime to apply H with two-map : (%.3lf +- %.4lf)ms",
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
