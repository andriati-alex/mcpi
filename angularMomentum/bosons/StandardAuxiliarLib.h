#ifndef _StandardAuxiliarLib_h
#define _StandardAuxiliarLib_h

#include <math.h>
#include "DataTypesDefinition.h"





void sepline()
{

/** print in current screen a separation line **/

    printf("\n=======================================");
    printf("=======================================\n");
}



void TimePrint(double t)
{
    
/** format and print time in days / hours / minutes given 't' in seconds **/

    int
        tt = (int) t,
        days  = 0,
        hours = 0,
        mins  = 0;

    if ( tt / 86400 > 0 )
    { days  = tt / 86400; tt = tt % 86400; }

    if ( tt / 3600  > 0 )
    { hours = tt / 3600;  tt = tt % 3600;  }

    if ( tt / 60    > 0 )
    { mins  = tt / 60;    tt = tt % 60;    }

    if (days > 0)
    {
        printf("%d day(s) %d hour(s) %d minute(s)", days, hours, mins);
    }
    else
    {
        if (hours > 0) printf("%d hour(s) %d minute(s)",hours,mins);
        else printf("%d minute(s)",mins);
    }
}



Rarray rarrDef(unsigned int n)
{
    double * ptr;

    ptr = (double * ) malloc(n*sizeof(double));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for double");
        printf(" Size requested : %d double elements\n\n",n);
        exit(EXIT_FAILURE);
    }
    return ptr;
}



Carray carrDef(unsigned int n)
{
    double complex * ptr;

    ptr = (double complex * ) malloc(n * sizeof(double complex));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for complex.");
        printf(" Size requested : %d double complex elements\n\n",n);
        exit(EXIT_FAILURE);
    }
    return ptr;
}



Iarray iarrDef(unsigned int n)
{
    int * ptr;

    ptr = (int * ) malloc(n * sizeof(int));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for integer.");
        printf(" Size requested : %d integer elements\n\n",n);
        exit(EXIT_FAILURE);
    }
    return ptr;
}



Farray farrDef(unsigned int n)
{
    char * ptr;

    ptr = (char * ) malloc(n*sizeof(char));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for array of char.");
        printf(" Size requested : %ld\n\n",n*sizeof(char));
        exit(EXIT_FAILURE);
    }
    return ptr;
}



Cmatrix cmatDef(unsigned int m, unsigned int n)
{

    int i;

    double complex ** ptr;

    ptr = (double complex ** ) malloc(m * sizeof(double complex *));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for (complex *)");
        printf(" Requested %d pointers to complex\n\n",m);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < m; i++) ptr[i] = carrDef(n);
    return ptr;
}



void initGuess(unsigned int nc, Carray C)
{

/** SET UP A NORMALIZED VECTOR WITH RANDOM VALUES **/

    int
        i;

    double
        sum,
        real,
        imag;

    sum = 0;
    for (i = 0; i < nc; i++)
    {
        real = (i % 3) - 1 + 4 * (i % 5);
        imag = 3 - (i % 7) + 2 * (i % 3);
        // real = ((i % 4) * (3.1864 - 7 % (i+1)));
        // imag = ((i * i) % 13 - 15.34576 * (i % 3));
        sum = sum + real*real + imag*imag;
        C[i] = real + I * imag;
    }
    for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);
}



double complex carrDot(int n, Carray v1, Carray v2)
{

/** Convetional L2 scalar product for complex vectors **/

    int i;

    double complex z = 0;
    for (i = 0; i < n; i++) z = z + conj(v1[i]) * v2[i];
    return z;
}



double carrNorm(int n, Carray v)
{

/** Conventional L2 norm of complex vectors **/

    int i;

    double mod = 0;
    for (i = 0; i < n; i++)
    {
        mod = mod + creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
    return sqrt(mod);
}



int assert_BoseNocc(int Norbs, int Npar, Iarray occ)
{

/** Auxiliar function to assess if the occupation numbers are valid, i.e
    no negative entries and conserving the total number of particles **/

    int
        i,
        N;

    N = 0;
    for (i = 0; i < Norbs; i++)
    {
        if (occ[i] < 0) return 0;
        N = N + occ[i];
    }
    if (N != Npar) return 0;
    return 1; // Everything fine
}



int assert_FermiNocc(int Norbs, int Npar, Farray occ)
{

/** Auxiliar function to assess if the occupation numbers are valid, i.e
    no negative entries and conserving the total number of particles **/

    int
        i,
        N;

    N = 0;
    for (i = 0; i < Norbs; i++)
    {
        if (occ[i] < 0 || occ[i] > 1) return 0;
        N = N + occ[i];
    }
    if (N != Npar) return 0;
    return 1; // Everything fine
}



FILE * openFileWrite(char fname [])
{
    FILE
        * arq;

    arq = fopen(fname,"w");
    if (arq == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n",fname);
        exit(EXIT_FAILURE);
    }
    return arq;
}



FILE * openFileRead(char fname [])
{
    FILE
        * arq;

    arq = fopen(fname,"r");
    if (arq == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n",fname);
        exit(EXIT_FAILURE);
    }
    return arq;
}



void ReachNewLine(FILE * f)
{

/** Read until get new line in a opened file. **/

    char
        sentinel;

    while (1)
    {
        fscanf(f,"%c",&sentinel);
        if (sentinel == '\n' || sentinel == EOF) return;
    }
}



void carrAppend(FILE * f, int M, Carray v)
{

/** Given a opened file write complex array in the succeeding lines **/

    int
        j;

    double
        real,
        imag;

    if (f == NULL)
    {
        printf("\n\nERROR: NULL file in carrAppend routine\n\n");
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < M; j ++)
    {
        real = creal(v[j]);
        imag = cimag(v[j]);

        if (imag >= 0) fprintf(f, "\n(%.15E+%.15Ej)", real, imag);
        else           fprintf(f, "\n(%.15E%.15Ej)", real, imag);
    }
}



unsigned int NumberOfLines(char fname [])
{

    int
        i;

    char
        c;

    FILE
        * arq;

    i = 0;
    arq = openFileRead(fname);
    // jump comment line which are initiated by #
    while ((c = getc(arq)) != EOF)
    {
        if (c == '\n') i++;
    }
    fclose(arq);
    return i;
}



#endif
