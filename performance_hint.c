#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include <omp.h>

#define PI 3.14159265359





int fac(int n)
{   // return n !
    int n_fac = 1;
    for (int i = 1; i < n; i++) n_fac = n_fac * (i + 1);
    return n_fac;
}





int NC(int N, int M)
{   // Return # of possible configurations of N particles in M orbitals
    int i, n = 1;
    if  (M > N)
    {
        for (i = N + M - 1; i > M - 1; i --) n = n * i;
        return n / fac(N);
    }
    else
    {
        for (i = N + M - 1; i > N; i --) n = n * i;
        return n / fac(M - 1);
    }
}





void IndexToFock(int k, int N, int M, int ** NCmat, int * v)
{
    int x;

    int  i,
         m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    /* -------------------------------------------------------------
     * Put the particles in orbitals while has combinations to spend
    ----------------------------------------------------------------- */
    while ( k > 0 )
    {
        while ( k - NCmat[N][m] < 0 ) m = m - 1;
        x = k - NCmat[N][m];
        while ( x >= 0 )
        {
            v[m] = v[m] + 1; // One more particle in orbital m
            N = N - 1;       // Less one particle to setup
            k = x;
            x = x - NCmat[N][m];
        }
    }

    /* -------------------------------------------------------------
     * Put the particles in orbitals while has combinations to spend
    ----------------------------------------------------------------- */
    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}





int FockToIndex(int N, int M, int ** NCmat, int * v)
{
    int i, n;
    int k = 0;

    /* ---------------------------------------------------
     * Empty one by one orbital starting from the last one
    ------------------------------------------------------ */
    for (i = M - 1; i > 0; i--)
    {
        n = v[i]; // Number of particles in the orbital
        while (n > 0)
        {
            k = k + NCmat[N][i]; // number of combinations needed
            N = N - 1;           // decrease the number of particles
            n = n - 1;
        }
    }

    return k;
}



Carray carrDef(int n)
{
    double complex * ptr;

    ptr = (double complex * ) malloc( n * sizeof(double complex) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

Cmatrix cmatDef(int m, int n)
{

/** Complex matrix of m rows and n columns **/

    int i;

    double complex ** ptr;

    ptr = (double complex ** ) malloc( m * sizeof(double complex *) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for (complex *)\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = carrDef(n);

    return ptr;
}



void cprint(double complex z)
{

/** printf complex number in coordinates with 2 decimal digits**/

    printf("(%9.2E,%9.2E )", creal(z), cimag(z));
}





void carr_txt(char fname [], int M, Carray v)
{

/** Record a array of complex elements in a text file in a
  * suitable format to import as numpy array with python.
**/

    int
        j;

    double
        real,
        imag;

    FILE
        * data_file = fopen(fname, "w");



    if (data_file == NULL)
    {
        printf("\n\n\n\tERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < M - 1; j ++)
    {

        real = creal(v[j]);
        imag = cimag(v[j]);

        if (imag >= 0) fprintf(data_file, "(%.15E+%.15Ej)", real, imag);
        else           fprintf(data_file, "(%.15E%.15Ej)", real, imag);

        fprintf(data_file, "\n");
    }

    real = creal(v[M-1]);
    imag = cimag(v[M-1]);

    if (imag >= 0) fprintf(data_file, "(%.15E+%.15Ej)", real, imag);
    else           fprintf(data_file, "(%.15E%.15Ej)", real, imag);

    fclose(data_file);
}





void cmat_txt (char fname [], int m, int n, Cmatrix A)
{

    int
        i,
        j;

    double
        real,
        imag;

    FILE
        * f = fopen(fname, "w");



    if (f == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\n\tERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }



    for (i = 0; i < m - 1; i++)
    {

        for (j = 0; j < n; j++)
        {

            real = creal(A[i][j]);
            imag = cimag(A[i][j]);

            if (imag >= 0) fprintf(f, "(%.15E+%.15Ej) ", real, imag);
            else           fprintf(f, "(%.15E%.15Ej) ", real, imag);
        }

        fprintf(f, "\n");
    }

    for (j = 0; j < n; j++)
    {

        real = creal(A[m-1][j]);
        imag = cimag(A[m-1][j]);

        if (imag >= 0) fprintf(f, "(%.15E+%.15Ej) ", real, imag);
        else           fprintf(f, "(%.15E%.15Ej) ", real, imag);
    }

    fclose(f);
}



double complex Csimps(int n, Carray f, double h)
{

    int
        i;

    double complex
        sum;

    sum = 0;

    if (n < 3)
    {
        printf("\n\n\tERROR : less than 3 point to integrate by simps !\n\n");
        exit(EXIT_FAILURE);
    }

    if (n % 2 == 0)
    {

    //  Case the number of points is even then must integrate the last
    //  chunk using simpson's 3/8 rule to maintain accuracy

        for (i = 0; i < (n - 4); i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }

        sum = sum * h / 3; // End 3-point simpsons intervals
        sum = sum + (f[n-4] + 3 * (f[n-3] + f[n-2]) + f[n-1]) * 3 * h / 8;

    }

    else
    {

        for (i = 0; i < n - 2; i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }

        sum = sum * h / 3; // End 3-point simpsons intervals

    }

    return sum;

}





int ** MountNCmat(int N, int M)
{   // Matrix of all possible configurations with
    // # particles < N and M < # of orbitals
    int i,
        j;

    int ** NCmat = (int ** ) malloc((N + 1) * sizeof(int * ));

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i] = (int * ) malloc((M + 1) * sizeof(int));
    }

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i][0] = 0;
        NCmat[i][1] = 1;
        for (j = 2; j < M + 1; j++) NCmat[i][j] = NC( i , j );
    }

    return NCmat;
}





int ** MountFocks(int N, int M, int ** NCmat)
{   // All possible occupation vectors organized in rows
    int k;

    int ** ItoFock = (int **) malloc(NC(N, M) * sizeof(int *));

    for (k = 0; k < NC(N, M); k++)
    {
        ItoFock[k] = (int * ) malloc(M * sizeof(int));
        IndexToFock(k, N, M, NCmat, ItoFock[k]);
    }

    return ItoFock;
}





void SetupHint_X (int Morb, int Mpos, Cmatrix Omat, double dx, double g,
     Carray Hint)
{

/** Configure two-body hamiltonian matrix elements in a chosen orbital basis
  * for contact interactions
  *
  * Output parameter : Hint
  *
  *************************************************************************/

    int i,
        k,
        s,
        q,
        l,
        M,
        M2,
        M3,
        ik,
        iq,
        is,
        il;

    double complex
        Integral;

    Carray
        orb,
        toInt;

    M  = Morb;
    M2 = M * M;
    M3 = M * M2;

#pragma omp parallel private(i,k,s,q,l,ik,iq,is,il,Integral,orb,toInt) \
    firstprivate(g,dx)

    {

    toInt = carrDef(Mpos);

    orb = carrDef(Morb*Mpos);

    for (k = 0; k < Morb; k++)
    {
        for (i = 0; i < Mpos; i++) orb[k * Mpos + i] = Omat[k][i];
    }

#pragma omp for schedule(static)
    for (k = 0; k < Morb; k++)
    {

        ik = k * Mpos;

        for (i = 0; i < Mpos; i++)
        {
            toInt[i] = conj(orb[ik+i]*orb[ik+i]) * orb[ik+i]*orb[ik+i];
        }

        Hint[k * (1 + M + M2 + M3)] = g * Csimps(Mpos,toInt,dx);

        for (s = k + 1; s < Morb; s++)
        {

            is = s * Mpos;

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(orb[ik+i]*orb[is+i]) * orb[ik+i]*orb[ik+i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[k + s * M + k * M2 + k * M3] = Integral;
            Hint[s + k * M + k * M2 + k * M3] = Integral;
            Hint[k + k * M + k * M2 + s * M3] = conj(Integral);
            Hint[k + k * M + s * M2 + k * M3] = conj(Integral);

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(orb[is+i]*orb[ik+i]) * orb[is+i]*orb[is+i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[s + k * M + s * M2 + s * M3] = Integral;
            Hint[k + s * M + s * M2 + s * M3] = Integral;
            Hint[s + s * M + s * M2 + k * M3] = conj(Integral);
            Hint[s + s * M + k * M2 + s * M3] = conj(Integral);

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(orb[ik+i]*orb[is+i]) * orb[is+i]*orb[ik+i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[k + s * M + s * M2 + k * M3] = Integral;
            Hint[s + k * M + s * M2 + k * M3] = Integral;
            Hint[s + k * M + k * M2 + s * M3] = Integral;
            Hint[k + s * M + k * M2 + s * M3] = Integral;

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(orb[ik+i]*orb[ik+i]) * orb[is+i]*orb[is+i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[k + k * M + s * M2 + s * M3] = Integral;
            Hint[s + s * M + k * M2 + k * M3] = conj(Integral);

            for (q = s + 1; q < Morb; q++)
            {

                iq = q * Mpos;

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(orb[ik+i]*orb[is+i]) * \
                               orb[iq+i]*orb[ik+i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[k + s * M + q * M2 + k * M3] = Integral;
                Hint[k + s * M + k * M2 + q * M3] = Integral;
                Hint[s + k * M + k * M2 + q * M3] = Integral;
                Hint[s + k * M + q * M2 + k * M3] = Integral;

                Hint[k + q * M + s * M2 + k * M3] = conj(Integral);
                Hint[k + q * M + k * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + k * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + s * M2 + k * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(orb[is+i]*orb[ik+i]) * \
                               orb[iq+i]*orb[is+i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[s + k * M + q * M2 + s * M3] = Integral;
                Hint[k + s * M + q * M2 + s * M3] = Integral;
                Hint[k + s * M + s * M2 + q * M3] = Integral;
                Hint[s + k * M + s * M2 + q * M3] = Integral;

                Hint[s + q * M + k * M2 + s * M3] = conj(Integral);
                Hint[s + q * M + s * M2 + k * M3] = conj(Integral);
                Hint[q + s * M + s * M2 + k * M3] = conj(Integral);
                Hint[q + s * M + k * M2 + s * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(orb[iq+i]*orb[is+i]) * \
                               orb[ik+i]*orb[iq+i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[q + s * M + k * M2 + q * M3] = Integral;
                Hint[q + s * M + q * M2 + k * M3] = Integral;
                Hint[s + q * M + q * M2 + k * M3] = Integral;
                Hint[s + q * M + k * M2 + q * M3] = Integral;

                Hint[k + q * M + s * M2 + q * M3] = conj(Integral);
                Hint[k + q * M + q * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + s * M2 + q * M3] = conj(Integral);
                Hint[q + k * M + q * M2 + s * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(orb[ik+i]*orb[ik+i]) * \
                               orb[iq+i]*orb[is+i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[k + k * M + q * M2 + s * M3] = Integral;
                Hint[k + k * M + s * M2 + q * M3] = Integral;
                Hint[q + s * M + k * M2 + k * M3] = conj(Integral);
                Hint[s + q * M + k * M2 + k * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(orb[is+i]*orb[is+i]) * \
                               orb[ik+i]*orb[iq+i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[s + s * M + k * M2 + q * M3] = Integral;
                Hint[s + s * M + q * M2 + k * M3] = Integral;
                Hint[k + q * M + s * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + s * M2 + s * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(orb[iq+i]*orb[iq+i]) * \
                               orb[ik+i]*orb[is+i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[q + q * M + k * M2 + s * M3] = Integral;
                Hint[q + q * M + s * M2 + k * M3] = Integral;
                Hint[k + s * M + q * M2 + q * M3] = conj(Integral);
                Hint[s + k * M + q * M2 + q * M3] = conj(Integral);

                for (l = q + 1; l < Morb; l++)
                {

                    il = l * Mpos;

                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(orb[ik+i] * orb[is+i]) * \
                                   orb[iq+i] * orb[il+i];
                    }

                    Integral = g * Csimps(Mpos,toInt,dx);

                    Hint[k + s * M + q * M2 + l * M3] = Integral;
                    Hint[k + s * M + l * M2 + q * M3] = Integral;
                    Hint[s + k * M + q * M2 + l * M3] = Integral;
                    Hint[s + k * M + l * M2 + q * M3] = Integral;

                    Hint[q + l * M + k * M2 + s * M3] = conj(Integral);
                    Hint[l + q * M + k * M2 + s * M3] = conj(Integral);
                    Hint[l + q * M + s * M2 + k * M3] = conj(Integral);
                    Hint[q + l * M + s * M2 + k * M3] = conj(Integral);

                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(orb[ik+i] * orb[iq+i]) * \
                                   orb[is+i] * orb[il+i];
                    }

                    Integral = g * Csimps(Mpos,toInt,dx);

                    Hint[k + q * M + s * M2 + l * M3] = Integral;
                    Hint[k + q * M + l * M2 + s * M3] = Integral;
                    Hint[q + k * M + s * M2 + l * M3] = Integral;
                    Hint[q + k * M + l * M2 + s * M3] = Integral;

                    Hint[s + l * M + k * M2 + q * M3] = conj(Integral);
                    Hint[s + l * M + q * M2 + k * M3] = conj(Integral);
                    Hint[l + s * M + q * M2 + k * M3] = conj(Integral);
                    Hint[l + s * M + k * M2 + q * M3] = conj(Integral);

                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(orb[ik+i] * orb[il+i]) * \
                                   orb[is+i] * orb[iq+i];
                    }

                    Integral = g * Csimps(Mpos,toInt,dx);

                    Hint[k + l * M + s * M2 + q * M3] = Integral;
                    Hint[k + l * M + q * M2 + s * M3] = Integral;
                    Hint[l + k * M + s * M2 + q * M3] = Integral;
                    Hint[l + k * M + q * M2 + s * M3] = Integral;

                    Hint[s + q * M + k * M2 + l * M3] = conj(Integral);
                    Hint[s + q * M + l * M2 + k * M3] = conj(Integral);
                    Hint[q + s * M + l * M2 + k * M3] = conj(Integral);
                    Hint[q + s * M + k * M2 + l * M3] = conj(Integral);

                }
            }
        }
    }

    free(toInt);

    free(orb);

    }

}





void SetupHint (int Morb, int Mpos, Cmatrix Omat, double dx, double g,
     Carray Hint)
{

/** Configure two-body hamiltonian matrix elements in a chosen orbital basis
  * for contact interactions
  *
  * Output parameter : Hint
  *
  *************************************************************************/

    int i,
        k,
        s,
        q,
        l,
        M,
        M2,
        M3;

    double complex
        Integral;

    Carray
        toInt;

    M  = Morb;
    M2 = M * M;
    M3 = M * M2;

    toInt = carrDef(Mpos);

    for (k = 0; k < Morb; k++)
    {

        for (i = 0; i < Mpos; i++)
        {
            toInt[i] = conj(Omat[k][i]*Omat[k][i]) * Omat[k][i]*Omat[k][i];
        }

        Hint[k * (1 + M + M2 + M3)] = g * Csimps(Mpos,toInt,dx);

        for (s = k + 1; s < Morb; s++)
        {

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(Omat[k][i]*Omat[s][i]) * Omat[k][i]*Omat[k][i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[k + s * M + k * M2 + k * M3] = Integral;
            Hint[s + k * M + k * M2 + k * M3] = Integral;
            Hint[k + k * M + k * M2 + s * M3] = conj(Integral);
            Hint[k + k * M + s * M2 + k * M3] = conj(Integral);

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(Omat[s][i]*Omat[k][i]) * Omat[s][i]*Omat[s][i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[s + k * M + s * M2 + s * M3] = Integral;
            Hint[k + s * M + s * M2 + s * M3] = Integral;
            Hint[s + s * M + s * M2 + k * M3] = conj(Integral);
            Hint[s + s * M + k * M2 + s * M3] = conj(Integral);

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(Omat[k][i]*Omat[s][i]) * Omat[s][i]*Omat[k][i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[k + s * M + s * M2 + k * M3] = Integral;
            Hint[s + k * M + s * M2 + k * M3] = Integral;
            Hint[s + k * M + k * M2 + s * M3] = Integral;
            Hint[k + s * M + k * M2 + s * M3] = Integral;

            for (i = 0; i < Mpos; i++)
            {
                toInt[i] = conj(Omat[k][i]*Omat[k][i]) * Omat[s][i]*Omat[s][i];
            }

            Integral = g * Csimps(Mpos,toInt,dx);

            Hint[k + k * M + s * M2 + s * M3] = Integral;
            Hint[s + s * M + k * M2 + k * M3] = conj(Integral);

            for (q = s + 1; q < Morb; q++)
            {

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(Omat[k][i]*Omat[s][i]) * \
                               Omat[q][i]*Omat[k][i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[k + s * M + q * M2 + k * M3] = Integral;
                Hint[k + s * M + k * M2 + q * M3] = Integral;
                Hint[s + k * M + k * M2 + q * M3] = Integral;
                Hint[s + k * M + q * M2 + k * M3] = Integral;

                Hint[k + q * M + s * M2 + k * M3] = conj(Integral);
                Hint[k + q * M + k * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + k * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + s * M2 + k * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(Omat[s][i]*Omat[k][i]) * \
                               Omat[q][i]*Omat[s][i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[s + k * M + q * M2 + s * M3] = Integral;
                Hint[k + s * M + q * M2 + s * M3] = Integral;
                Hint[k + s * M + s * M2 + q * M3] = Integral;
                Hint[s + k * M + s * M2 + q * M3] = Integral;

                Hint[s + q * M + k * M2 + s * M3] = conj(Integral);
                Hint[s + q * M + s * M2 + k * M3] = conj(Integral);
                Hint[q + s * M + s * M2 + k * M3] = conj(Integral);
                Hint[q + s * M + k * M2 + s * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(Omat[q][i]*Omat[s][i]) * \
                               Omat[k][i]*Omat[q][i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[q + s * M + k * M2 + q * M3] = Integral;
                Hint[q + s * M + q * M2 + k * M3] = Integral;
                Hint[s + q * M + q * M2 + k * M3] = Integral;
                Hint[s + q * M + k * M2 + q * M3] = Integral;

                Hint[k + q * M + s * M2 + q * M3] = conj(Integral);
                Hint[k + q * M + q * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + s * M2 + q * M3] = conj(Integral);
                Hint[q + k * M + q * M2 + s * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(Omat[k][i]*Omat[k][i]) * \
                               Omat[q][i]*Omat[s][i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[k + k * M + q * M2 + s * M3] = Integral;
                Hint[k + k * M + s * M2 + q * M3] = Integral;
                Hint[q + s * M + k * M2 + k * M3] = conj(Integral);
                Hint[s + q * M + k * M2 + k * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(Omat[s][i]*Omat[s][i]) * \
                               Omat[k][i]*Omat[q][i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[s + s * M + k * M2 + q * M3] = Integral;
                Hint[s + s * M + q * M2 + k * M3] = Integral;
                Hint[k + q * M + s * M2 + s * M3] = conj(Integral);
                Hint[q + k * M + s * M2 + s * M3] = conj(Integral);

                for (i = 0; i < Mpos; i++)
                {
                    toInt[i] = conj(Omat[q][i]*Omat[q][i]) * \
                               Omat[k][i]*Omat[s][i];
                }

                Integral = g * Csimps(Mpos,toInt,dx);

                Hint[q + q * M + k * M2 + s * M3] = Integral;
                Hint[q + q * M + s * M2 + k * M3] = Integral;
                Hint[k + s * M + q * M2 + q * M3] = conj(Integral);
                Hint[s + k * M + q * M2 + q * M3] = conj(Integral);

                for (l = q + 1; l < Morb; l++)
                {

                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(Omat[k][i] * Omat[s][i]) * \
                                   Omat[q][i] * Omat[l][i];
                    }

                    Integral = g * Csimps(Mpos,toInt,dx);

                    Hint[k + s * M + q * M2 + l * M3] = Integral;
                    Hint[k + s * M + l * M2 + q * M3] = Integral;
                    Hint[s + k * M + q * M2 + l * M3] = Integral;
                    Hint[s + k * M + l * M2 + q * M3] = Integral;

                    Hint[q + l * M + k * M2 + s * M3] = conj(Integral);
                    Hint[l + q * M + k * M2 + s * M3] = conj(Integral);
                    Hint[l + q * M + s * M2 + k * M3] = conj(Integral);
                    Hint[q + l * M + s * M2 + k * M3] = conj(Integral);

                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(Omat[k][i] * Omat[q][i]) * \
                                   Omat[s][i] * Omat[l][i];
                    }

                    Integral = g * Csimps(Mpos,toInt,dx);

                    Hint[k + q * M + s * M2 + l * M3] = Integral;
                    Hint[k + q * M + l * M2 + s * M3] = Integral;
                    Hint[q + k * M + s * M2 + l * M3] = Integral;
                    Hint[q + k * M + l * M2 + s * M3] = Integral;

                    Hint[s + l * M + k * M2 + q * M3] = conj(Integral);
                    Hint[s + l * M + q * M2 + k * M3] = conj(Integral);
                    Hint[l + s * M + q * M2 + k * M3] = conj(Integral);
                    Hint[l + s * M + k * M2 + q * M3] = conj(Integral);

                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(Omat[k][i] * Omat[l][i]) * \
                                   Omat[s][i] * Omat[q][i];
                    }

                    Integral = g * Csimps(Mpos,toInt,dx);

                    Hint[k + l * M + s * M2 + q * M3] = Integral;
                    Hint[k + l * M + q * M2 + s * M3] = Integral;
                    Hint[l + k * M + s * M2 + q * M3] = Integral;
                    Hint[l + k * M + q * M2 + s * M3] = Integral;

                    Hint[s + q * M + k * M2 + l * M3] = conj(Integral);
                    Hint[s + q * M + l * M2 + k * M3] = conj(Integral);
                    Hint[q + s * M + l * M2 + k * M3] = conj(Integral);
                    Hint[q + s * M + k * M2 + l * M3] = conj(Integral);

                }
            }
        }
    }

    free(toInt);

}





void TBrho(int N, int M, int ** NCmat, int ** IF, Carray C, Carray rho)
{

    int i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        t,
        nc,
        M2,
        M3,
        * v;
    
    double
        mod2,   // |Cj| ^ 2
        sqrtOf; // Factors from the action of creation/annihilation

    double complex
        RHO;





    // Auxiliar to memory access
    M2 = M * M;
    M3 = M * M * M;

    nc = NCmat[N][M];



    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {

        RHO = 0;

        #pragma omp parallel for private(i,mod2) reduction(+:RHO)
        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO  = RHO + mod2 * IF[i][k] * (IF[i][k] - 1);
        }

        rho[k + M * k + M2 * k + M3 * k] = RHO;
    }



    /* ---------------------------------------------
     * Rule 2: Creation on k s / Annihilation on k s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {

            RHO = 0;

            #pragma omp parallel for private(i,mod2) reduction(+:RHO)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[i][k] * IF[i][s];
            }

            // commutation of bosonic operators is used
            // to fill elements by exchange  of indexes
            rho[k + s * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + k * M3] = RHO;
            rho[k + s * M + s * M2 + k * M3] = RHO;
        }
    }



    /* ---------------------------------------------
     * Rule 3: Creation on k k / Annihilation on q q
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = sqrt((double)(v[k]-1)*v[k]*(v[q]+1)*(v[q]+2));
                    v[k] -= 2;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            // Use 2-index-'hermiticity'
            rho[k + k * M + q * M2 + q * M3] = RHO;
            rho[q + q * M + k * M2 + k * M3] = conj(RHO);
        }
    }



    /* ---------------------------------------------
     * Rule 4: Creation on k k / Annihilation on k l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            rho[k + k * M + k * M2 + l * M3] = RHO;
            rho[k + k * M + l * M2 + k * M3] = RHO;
            rho[l + k * M + k * M2 + k * M3] = conj(RHO);
            rho[k + l * M + k * M2 + k * M3] = conj(RHO);
        }
    }



    /* ---------------------------------------------
     * Rule 5: Creation on k s / Annihilation on s s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = v[s] * sqrt((double)v[k]*(v[s]+1));
                    v[k] -= 1;
                    v[s] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            rho[k + s * M + s * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + s * M3] = RHO;
            rho[s + s * M + s * M2 + k * M3] = conj(RHO);
            rho[s + s * M + k * M2 + s * M3] = conj(RHO);
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (k = q + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (l = q + 1; l < M; l++)
        {
            for (k = l + 1; k < M; k++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.0: Creation on k s / Annihilation on s l (s < k < l)
    ------------------------------------------------------------------- */
    for (s = 0; s < M; s++)
    {
        for (k = s + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.1: Creation on k s / Annihilation on s l (k < s < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (l = s + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.2: Creation on k s / Annihilation on s l (k < l < s)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            for (s = l + 1; s < M; s++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k]*(v[l]+1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* ---------------------------------------------
     * Rule 8: Creation on k s / Annihilation on q l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = 0; s < M; s++)
        {
            if (s == k) continue;

            for (q = 0; q < M; q++)
            {
                if (q == s || q == k) continue;

                for (l = 0; l < M; l ++)
                {

                    RHO = 0;

                    if (l == k || l == s || l == q) continue;

                    #pragma omp parallel private(i,j,t,sqrtOf,v) \
                    reduction(+:RHO)
                    {
                        v = (int *) malloc(M * sizeof(int));
                        #pragma omp for
                        for (i = 0; i < nc; i++)
                        {
                            if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                            for (t = 0; t < M; t++) v[t] = IF[i][t];
                            sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                            v[k] -= 1;
                            v[s] -= 1;
                            v[q] += 1;
                            v[l] += 1;
                            j = FockToIndex(N, M, NCmat, v);
                            RHO += conj(C[i]) * C[j] * sqrtOf;
                        }
                        free(v);
                    }

                    rho[k + s * M + q * M2 + l * M3] = RHO;
                }   // Finish l loop
            }       // Finish q loop
        }           // Finish s loop
    }               // Finish k loop


    /*       ------------------- END OF ROUTINE -------------------       */
}





doublec nonlinear (int M, int k, int n, double g, Cmatrix Orb,
        Cmatrix Rinv, Carray R2, Cmatrix Ho, Carray Hint )
{

/** For a orbital 'k' computed at discretized position 'n' calculate
  * the right-hand-side part of MCTDHB orbital's equation of  motion
  * that is nonlinear, part because of projections that made the eq.
  * an integral-differential equation, and other part due to contact
  * interactions. Assume that Rinv, R2 are  defined  by  the  set of
  * configuration-state coefficients as the inverse of  one-body and
  * two-body density matrices respectively. Ho and Hint are  assumed
  * to be defined accoding to 'Orb' variable as well.
  *
**/



    int a,
        j,
        s,
        q,
        l,
        M2,
        M3,
        ind;



    double complex
        G,
        X;



    X = 0;
    M2 = M * M;
    M3 = M * M * M;



    for (s = 0; s < M; s++)
    {
        // Subtract one-body projection
        X = X - Ho[s][k] * Orb[s][n];

        for (a = 0; a < M; a++)
        {

            for (q = 0; q < M; q++)
            {
                // Particular case with the two last indices equals
                // to take advantage of the symmetry afterwards

                G = Rinv[k][a] * R2[a + M*s + M2*q + M3*q];

                // Sum interacting part contribution
                X = X + g * G * conj(Orb[s][n]) * Orb[q][n] * Orb[q][n];

                // Subtract interacting projection
                for (j = 0; j < M; j++)
                {
                    ind = j + s * M + q * M2 + q * M3;
                    X = X - G * Orb[j][n] * Hint[ind];
                }

                for (l = q + 1; l < M; l++)
                {
                    G = 2 * Rinv[k][a] * R2[a + M*s + M2*q + M3*l];

                    // Sum interacting part
                    X = X + g * G * conj(Orb[s][n]) * Orb[l][n] * Orb[q][n];

                    // Subtract interacting projection
                    for (j = 0; j < M; j++)
                    {
                        ind = j + s * M + l * M2 + q * M3;
                        X = X - G * Orb[j][n] * Hint[ind];
                    }
                }
            }
        }
    }

    return X;
}





doublec nonlinear_X (int M, int k, int n, double g, Cmatrix Orb,
        Cmatrix Rinv, Carray R2, Cmatrix Ho, Carray Hint )
{

/** For a orbital 'k' computed at discretized position 'n' calculate
  * the right-hand-side part of MCTDHB orbital's equation of  motion
  * that is nonlinear, part because of projections that made the eq.
  * an integral-differential equation, and other part due to contact
  * interactions. Assume that Rinv, R2 are  defined  by  the  set of
  * configuration-state coefficients as the inverse of  one-body and
  * two-body density matrices respectively. Ho and Hint are  assumed
  * to be defined accoding to 'Orb' variable as well.
  */



    int a,
        j,
        s,
        q,
        l,
        M2,
        M3,
        ind;



    double complex
        G,
        Ginv,
        X;



    X = 0;
    M2 = M * M;
    M3 = M * M * M;



    for (a = 0; a < M; a++)
    {
        // Subtract one-body projection
        X = X - Ho[a][k] * Orb[a][n];

        for (q = 0; q < M; q++)
        {
            // Particular case with the two last indices equals
            // to take advantage of the symmetry afterwards

            G = Rinv[k][a] * R2[a + M*a + M2*q + M3*q];

            // Sum interacting part contribution
            X = X + g * G * conj(Orb[a][n]) * Orb[q][n] * Orb[q][n];

            // Subtract interacting projection
            for (j = 0; j < M; j++)
            {
                ind = j + a * M + q * M2 + q * M3;
                X = X - G * Orb[j][n] * Hint[ind];
            }

            for (l = q + 1; l < M; l++)
            {
                G = 2 * Rinv[k][a] * R2[a + M*a + M2*q + M3*l];

                // Sum interacting part
                X = X + g * G * conj(Orb[a][n]) * Orb[l][n] * Orb[q][n];

                // Subtract interacting projection
                for (j = 0; j < M; j++)
                {
                    ind = j + a * M + l * M2 + q * M3;
                    X = X - G * Orb[j][n] * Hint[ind];
                }
            }
        }

        for (s = a + 1; s < M; s++)
        {

            for (q = 0; q < M; q++)
            {
                // Particular case with the two last indices equals
                // to take advantage of the symmetry afterwards

                G = Rinv[k][a] * R2[a + M*s + M2*q + M3*q];
                Ginv = Rinv[k][s] * R2[a + M*s + M2*q + M3*q];

                // Sum interacting part contribution
                X = X + g * (G*conj(Orb[s][n]) + Ginv*conj(Orb[a][n])) * \
                    Orb[q][n]*Orb[q][n];

                // Subtract interacting projection
                for (j = 0; j < M; j++)
                {
                    ind = j + s * M + q * M2 + q * M3;
                    X = X - G * Orb[j][n] * Hint[ind];
                    ind = j + a * M + q * M2 + q * M3;
                    X = X - Ginv * Orb[j][n] * Hint[ind];
                }

                for (l = q + 1; l < M; l++)
                {
                    G = 2 * Rinv[k][a] * R2[a + M*s + M2*q + M3*l];
                    Ginv = 2 * Rinv[k][s] * R2[a + M*s + M2*q + M3*l];

                    // Sum interacting part
                    X = X + g * (G*conj(Orb[s][n]) + Ginv*conj(Orb[a][n])) * \
                            Orb[l][n]*Orb[q][n];

                    // Subtract interacting projection
                    for (j = 0; j < M; j++)
                    {
                        ind = j + s * M + l * M2 + q * M3;
                        X = X - G * Orb[j][n] * Hint[ind];
                        ind = j + a * M + l * M2 + q * M3;
                        X = X - Ginv * Orb[j][n] * Hint[ind];
                    }
                }
            }
        }
    }

    return X;
}










int main(int argc, char * argv[])
{

    omp_set_num_threads(omp_get_max_threads() / 2);

    int
        i,
        j,
        k,
        nc,
        Morb,
        Mpos,
        Npar,
        ** NCmat,
        ** IF;

    double
        sum,
        start,
        time_used;

    Carray
        C,
        Hint_X,
        rho2,
        Hint;

    Cmatrix
        Ho,
        dOdt,
        rho_inv,
        orb;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of grid points and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Mpos);
    sscanf(argv[2],"%d",&Morb);

    Npar = 4;

    Hint = (doublec * ) malloc(Morb*Morb*Morb*Morb*sizeof(doublec));
    Hint_X = (doublec * ) malloc(Morb*Morb*Morb*Morb*sizeof(doublec));
    rho2 = (doublec * ) malloc(Morb*Morb*Morb*Morb*sizeof(doublec));

    nc = NC(Npar,Morb);
    C = (doublec *) malloc(nc * sizeof(doublec));
    sum = 0.0;
    for (i = 0; i < nc; i++)
    {
        C[i] = sin( 20 * ((double) i) / nc) * (i % 13) - (i % 8) * I;
        sum = sum + creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
    }
    for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);

    NCmat = MountNCmat(Npar,Morb);
    IF = MountFocks(Npar,Morb,NCmat);

    TBrho(Npar,Morb,NCmat,IF,C,rho2);

    orb = cmatDef(Morb,Mpos);
    dOdt = cmatDef(Morb,Mpos);
    Ho = cmatDef(Morb,Morb);
    rho_inv = cmatDef(Morb,Morb);

    for (i = 0; i < Morb; i++)
    {
        for (j = 0; j < Mpos; j++)
        {
            orb[i][j] = ((i * j) % 8) * sin( (PI / Mpos) * j) - \
                        I * ((i + j) % 5) * cos(6 * (PI / Mpos) * j);
        }

        for (j = i + 1; j < Morb; j++)
        {
            Ho[i][j] = (i - j) + I * (i % (j + 1)) * 3;
            Ho[j][i] = conj(Ho[i][j]);
            rho_inv[i][j] = - 5 * sin(i * j / 3.0) + I * cos(1.0 * j - i);
            rho_inv[j][i] = conj(rho_inv[i][j]);
        }

        Ho[i][i] = 5.0 / (i + 1);
        rho_inv[i][i] =  -7.3456 / (i + 1);
    }


    start = omp_get_wtime();
    for (i = 0; i < 50; i++) SetupHint_X(Morb,Mpos,orb,0.02,13.764539,Hint_X);
    time_used = (double) (omp_get_wtime() - start) / 50;
    printf("\n\nTime to setup Hint_X : %.1lfms", time_used * 1000);

    start = omp_get_wtime();
    for (i = 0; i < 50; i++) SetupHint(Morb,Mpos,orb,0.02,13.764539,Hint);
    time_used = (double) (omp_get_wtime() - start) / 50;
    printf("\n\nTime to setup Hint : %.1lfms", time_used * 1000);

    start = omp_get_wtime();
    #pragma omp parallel for private(k, j) schedule(static)
    for (k = 0; k < Morb; k++)
    {
        for (j = 0; j < Mpos; j++)
            dOdt[k][j] = - I * \
            nonlinear(Morb, k, j, 12.4534, orb, rho_inv, rho2, Ho, Hint);
    }
    time_used = (double) (omp_get_wtime() - start);
    printf("\n\nTime to compute nonlinear : %.1lfms", time_used * 1000);

    cmat_txt ("dOdt.dat", Morb, Mpos, dOdt);



    start = omp_get_wtime();
    #pragma omp parallel for private(k, j) schedule(static)
    for (k = 0; k < Morb; k++)
    {
        for (j = 0; j < Mpos; j++)
            dOdt[k][j] = - I * \
            nonlinear_X(Morb, k, j, 12.4534, orb, rho_inv, rho2, Ho, Hint);
    }
    time_used = (double) (omp_get_wtime() - start);
    printf("\n\nTime to compute nonlinear_X : %.1lfms", time_used * 1000);

    cmat_txt ("dOdt_X.dat", Morb, Mpos, dOdt);



    carr_txt("Hint_X.dat",Morb*Morb*Morb*Morb,Hint_X);
    carr_txt("Hint.dat",Morb*Morb*Morb*Morb,Hint);
    

    free(Hint);
    free(Hint_X);

    for(i = 0; i < Morb; i++) free(orb[i]);
    free(orb);

    printf("\n\n");
    return 0;
}
