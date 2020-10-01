#ifndef _StandardAuxiliarLib_h
#define _StandardAuxiliarLib_h

#include <math.h>
#include "DataTypesDefinition.h"





void sepline()
{

/** print in the screen a separation line **/

    printf("\n=======================================");
    printf("=======================================\n");
}



void LAPACK_PROBLEM(int k, char funcName [])
{
    printf("\n\nERROR IN LAPACK CALL IN FUNCTION %s\n\n",funcName);
    if (k < 0) printf("Illegal value in LAPACK_dstev parameter %d\n\n",-k);
    else       printf("LAPACK_dstev algorithm failed to converge\n\n");
    exit(EXIT_FAILURE);
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



Rmatrix rmatDef(unsigned int m, unsigned int n)
{

    int i;

    double ** ptr;

    ptr = (double ** ) malloc(m * sizeof(double *));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for real matrix\n\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < m; i++) ptr[i] = rarrDef(n);
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



void cmatFree(unsigned int m, Cmatrix M)
{
    int
        i;
    for (i = 0; i < m; i++) free(M[i]);
    free(M);
}



void rmatFree(unsigned int m, Rmatrix M)
{
    int
        i;
    for (i = 0; i < m; i++) free(M[i]);
    free(M);
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



void carrNormalize(unsigned int n, Carray v)
{

/** Normalize vector v according to L2 norm **/

    int
        i;
    double
        norm;
    norm = carrNorm(n,v);
    for (i = 0; i < n; i++) v[i] = v[i] / norm;
}



double rarrNorm(unsigned int n, Rarray x)
{

/** Conventional L2 norm of real vectors **/

    int
        i;
    double
        sum;

    sum = 0;
    for (i = 0; i < n; i++) sum = sum + x[i]*x[i];
    return sqrt(sum);
}



void rarrNormalize(unsigned int n, Rarray x)
{

/** Normalize vector x accorsing to L2 norm **/

    int
        i;
    double
        norm;
    norm = rarrNorm(n,x);
    for (i = 0; i < n; i++) x[i] = x[i] / norm;
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
    while ((c = getc(arq)) != EOF)
    {
        if (c == '\n') i++;
    }
    fclose(arq);
    return i;
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



void carr_txt(char fname [], int M, Carray v)
{

/** Record a array of complex elements in a text file **/

    int
        j;
    double
        real,
        imag;
    FILE *
        data_file;

    data_file = openFileWrite(fname);

    for (j = 0; j < M; j ++)
    {
        real = creal(v[j]);
        imag = cimag(v[j]);
        if (imag >= 0) fprintf(data_file, "(%.15E+%.15Ej)", real, imag);
        else           fprintf(data_file, "(%.15E%.15Ej)", real, imag);
        fprintf(data_file, "\n");
    }
}



void carr_input_txt(char fname [], int M, Carray v)
{
    int
        i,
        j;
    double
        real,
        imag;
    FILE
        * data_file;

    j = NumberOfLines(fname);
    if (j != M)
    {
        printf("\n\n\nINPUT ERROR : The number of lines in the input ");
        printf("file %s(%d lines) is not the required of ",fname,j);
        printf("%d lines\n\n\n",M);
        exit(EXIT_FAILURE);
    }
    data_file = openFileRead(fname);
    for (j = 0; j < M; j++)
    {
        i = fscanf(data_file,"(%lf%lfj)\n",&real,&imag);
        v[j] = real + I * imag;
    }
    fclose(data_file);
}



void carr_inline(FILE * f, int M, Carray v)
{

/** Given a opened file write complex array in current line of buffer **/

    int
        j;
    double
        real,
        imag;

    if (f == NULL)
    {
        printf("\n\n\nERROR: NULL file in carr_inline routine\n\n");
        exit(EXIT_FAILURE);
    }
    for (j = 0; j < M; j ++)
    {
        real = creal(v[j]);
        imag = cimag(v[j]);
        if (imag >= 0) fprintf(f,"(%.15E+%.15Ej) ",real,imag);
        else           fprintf(f,"(%.15E%.15Ej) ",real,imag);
    }
    fprintf(f,"\n");
}



void parLine(char fname [], int line, int * N, int * lmax, int * L,
             double * v, double * g)
{

/** READ A LINE OF INPUT FILE TO SET THE PARAMETERS FOR 1 SPECIES **/

    int
        i;

    FILE
        * in_file;

    in_file = openFileRead(fname);

    // jump lines to get to requested 'line'
    for (i = 0; i < line; i++) ReachNewLine(in_file);

    // read parameters
    fscanf(in_file,"%d",N);     // Number of particles
    fscanf(in_file,"%d",lmax);  // max. individual ang. momentum
    fscanf(in_file,"%d",L);     // total angular momentum constraint
    fscanf(in_file,"%lf",v);    // frame velocity
    fscanf(in_file,"%lf",g);    // contact interaction strength parameter

    fclose(in_file);
}



void mixParLine(char fname [], int line, int * NA, int * lmaxA, int * NB,
                int * lmaxB, int * L, double * v, double * mi, double g [])
{

/** READ A LINE OF INPUT FILE TO SET THE PARAMETERS FOR 2 SPECIES **/

    int
        i;

    FILE
        * in_file;

    in_file = openFileRead(fname);

    // jump lines to get to requested 'line'
    for (i = 0; i < line; i++) ReachNewLine(in_file);

    // read parameters
    fscanf(in_file,"%d",NA);    // Num. of particles A (aways bosons)
    fscanf(in_file,"%d",lmaxA); // max. individual ang. momemtum
    fscanf(in_file,"%d",NB);    // Num. of particles B (bosons or fermions)
    fscanf(in_file,"%d",lmaxB); // max. individual ang. momentum
    fscanf(in_file,"%d",L);     // Total ang. momentum constraint
    fscanf(in_file,"%lf",v);     // frame velocity
    fscanf(in_file,"%lf",mi);    // mass imbalance
    // CONTACT INTERACTION STRENGTH PARAMETERS
    fscanf(in_file,"%lf",&g[0]);    // of particles A
    fscanf(in_file,"%lf",&g[1]);    // of particles B (ignored for fermions)
    fscanf(in_file,"%lf",&g[2]);    // interspecies interaction

    fclose(in_file);
}



void parLine_time(char fname [], int line, int * N, int * lmax, int * L,
                  double * v, double * g, int * Nsteps, double * dt)
{

/** READ A LINE OF INPUT FILE TO SET THE PARAMETERS FOR 1 SPECIES **/

    int
        i;

    FILE
        * in_file;

    in_file = openFileRead(fname);

    // jump lines to get to requested 'line'
    for (i = 0; i < line; i++) ReachNewLine(in_file);

    // read parameters
    fscanf(in_file,"%d",N);     // Number of particles
    fscanf(in_file,"%d",lmax);  // max. individual ang. momentum
    fscanf(in_file,"%d",L);     // total angular momentum constraint
    fscanf(in_file,"%lf",v);    // frame velocity
    fscanf(in_file,"%lf",g);    // contact interaction strength parameter
    fscanf(in_file,"%d",Nsteps);
    fscanf(in_file,"%lf",dt);

    fclose(in_file);
}



void mixParLine_time(char fname [], int line, int * NA, int * lmaxA,
                     int * NB, int * lmaxB, int * L, double * v,
                     double * mi, double g [], int * Nsteps, double * dt)
{

/** READ A LINE OF INPUT FILE TO SET THE PARAMETERS FOR 2 SPECIES **/

    int
        i;

    FILE
        * in_file;

    in_file = openFileRead(fname);

    // jump lines to get to requested 'line'
    for (i = 0; i < line; i++) ReachNewLine(in_file);

    // read parameters
    fscanf(in_file,"%d",NA);    // Num. of particles A (aways bosons)
    fscanf(in_file,"%d",lmaxA); // max. individual ang. momemtum
    fscanf(in_file,"%d",NB);    // Num. of particles B (bosons or fermions)
    fscanf(in_file,"%d",lmaxB); // max. individual ang. momentum
    fscanf(in_file,"%d",L);     // Total ang. momentum constraint
    fscanf(in_file,"%lf",v);     // frame velocity
    fscanf(in_file,"%lf",mi);    // mass imbalance
    // CONTACT INTERACTION STRENGTH PARAMETERS
    fscanf(in_file,"%lf",&g[0]);    // of particles A
    fscanf(in_file,"%lf",&g[1]);    // of particles B (ignored for fermions)
    fscanf(in_file,"%lf",&g[2]);    // interspecies interaction
    fscanf(in_file,"%d",Nsteps);
    fscanf(in_file,"%lf",dt);

    fclose(in_file);
}



void rmatBlockMult(int n, int init, Rmatrix A, Rmatrix B, Rmatrix AB)
{

/** Perform the matrix multiplication of A(left) with B(right),
    both square,  beginning at index 'init' (lower right block)  **/

    int
        i,
        j,
        l;

    double
        x;

    for (i = init; i < n; i++)
    {
        for (j = init; j < n; j++)
        {
            x = 0;
            for (l = init; l < n; l++) x = x + A[i][l]*B[l][j];
            AB[i][j] = x;
        }
    }
}



void rmatTranspose(int n, Rmatrix A, Rmatrix AT)
{
    int
        i,
        j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++) AT[i][j] = A[j][i];
    }
}



void rmatCopy(int n, int init, Rmatrix inp, Rmatrix out)
{
    int
        i,
        j;

    for (i = init; i < n; i++)
    {
        for (j = init; j < n; j++) out[i][j] = inp[i][j];
    }
}



unsigned int numericalIdentityMatrix(int n, Rmatrix M)
{
/** Return boolean 1 if the matrix M is the Identity matrix in numerical
    precision and 0 otherwise. Default tolerance is 10^(-12) for the avg
    of the absolute valued of non diagonal elements                  **/
    int
        i,
        j;
    double
        offdiag;

    offdiag = 0;
    for (i = 0; i < n; i++)
    {
        for (j = i+1; j < n; j++)
        {
            offdiag = offdiag + (fabs(M[i][j]) + fabs(M[j][i]))/n;
        }
    }
    if (offdiag > 1E-12) return 0;
    return 1;
}



void QR_rec(int dim, int init, Rmatrix A, Rmatrix Qstep, Rmatrix Q, Rmatrix aux)
{

/** RECURSIVE IMPLEMENTATION USING HOUSEHOLDER REFLECTIONS FOLLOWING
    WIKIPEDIA CONVENTIONS. IT IS WORTH CONSULTING STOER AND BULIRSCH
    SEC. 3.7 (MORE SPECIFICALLY PAGES 227-228)                   **/

    int
        i,
        j,
        n;

    double
        alpha;

    Rarray
        x;

    if (init == dim-1) return;

    n = dim-init;   // dimension of bulk matrix A'
    x = rarrDef(n); // vector for new Q computation
    for (j = init; j < dim; j++) x[j-init] = A[j][init];
    if (x[0] < 0) alpha = rarrNorm(n,x);
    else          alpha = - rarrNorm(n,x);
    x[0] = x[0] - alpha; // update to get 'u'
    rarrNormalize(n,x);  // normalize it to get 'v'
    // SETUP NEW Q in 'Qstep'
    for (i = 0; i < init; i++)
    {
        Qstep[i][i] = 1;
        for (j = i+1; j < dim; j++)
        {
            Qstep[i][j] = 0;
            Qstep[j][i] = 0;
        }
    }
    for (i = init; i < dim; i++)
    {
        Qstep[i][i] = 1.0 - 2 * x[i-init]*x[i-init];
        for (j = i+1; j < dim; j++)
        {
            Qstep[i][j] = - 2 * x[i-init]*x[j-init];
            Qstep[j][i] = - 2 * x[j-init]*x[i-init];
        }
    }
    rmatBlockMult(dim,0,Qstep,Q,aux);   // one new multiplication
    rmatCopy(dim,0,aux,Q);              // update final Q
    rmatBlockMult(dim,init,Qstep,A,aux);
    rmatCopy(dim,init,aux,A);
    free(x);
    QR_rec(dim,init+1,A,Qstep,Q,aux);
}



void QRdecomp(unsigned int dim, Rmatrix A, Rmatrix Q)
{

/** RECURSIVE IMPLEMENTATION USING HOUSEHOLDER REFLECTIONS FOLLOWING
    WIKIPEDIA CONVENTIONS. IT IS WORTH CONSULTING STOER AND BULIRSCH
    SEC. 3.7 (MORE SPECIFICALLY PAGES 227-228)
    ================================================================
    Input  : Real square matrix 'A' and its dimension 'dim'
    Output : Real square matrix 'Q' which have orthonormal columns
    In the computations the input matrix 'A' is destroyed since it
    is used as workspace memory                                  **/

    int
        i,
        j;

    double
        alpha;

    Rarray
        x;

    Rmatrix
        Qstep,
        aux;

    if (dim < 2) return;

    aux = rmatDef(dim,dim);
    Qstep = rmatDef(dim,dim);

    x = rarrDef(dim);       // vector for Q computation
    for (j = 0; j < dim; j++) x[j] = A[j][0]; // copy the first column
    if (x[0] < 0) alpha = rarrNorm(dim,x);
    else          alpha = - rarrNorm(dim,x);
    x[0] = x[0] - alpha;    // update to get 'u'
    rarrNormalize(dim,x);   // normalize it to get 'v'
    // SETUP FIRST Q TO INITIATE RECURSION
    for (i = 0; i < dim; i++)
    {
        Q[i][i] = 1.0 - 2 * x[i]*x[i];
        for (j = i+1; j < dim; j++)
        {
            Q[i][j] = - 2 * x[i]*x[j];
            Q[j][i] = - 2 * x[j]*x[i];
        }
    }
    free(x);
    rmatBlockMult(dim,0,Q,A,aux);
    rmatCopy(dim,0,aux,A);
    QR_rec(dim,1,A,Qstep,Q,aux); // Finish recursive part
    // Final step - transpose
    rmatCopy(dim,0,Q,aux);
    rmatTranspose(dim,aux,Q);
    // Free memory
    rmatFree(dim,aux);
    rmatFree(dim,Qstep);
}



#endif
