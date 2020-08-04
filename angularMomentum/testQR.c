#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef double * Rarray;
typedef double ** Rmatrix;



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



Rmatrix rmatDef(unsigned int m, unsigned int n)
{

    int i;

    double ** ptr;

    ptr = (double ** ) malloc(m * sizeof(double *));
    if (ptr == NULL)
    {
        printf("\n\nMEMORY ERROR : malloc fail for matrix\n\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < m; i++) ptr[i] = rarrDef(n);
    return ptr;
}



double rarrNorm(unsigned int n, Rarray x)
{
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
    int
        i;
    double
        norm;
    norm = rarrNorm(n,x);
    for (i = 0; i < n; i++) x[i] = x[i] / norm;
}



void SpecialMult(int n, int init, Rmatrix A, Rmatrix B, Rmatrix AB)
{
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



void Transpose(int n, Rmatrix A, Rmatrix AT)
{
    int
        i,
        j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++) AT[i][j] = A[j][i];
    }
}



void RmatCopy(int n, int init, Rmatrix inp, Rmatrix out)
{
    int
        i,
        j;

    for (i = init; i < n; i++)
    {
        for (j = init; j < n; j++) out[i][j] = inp[i][j];
    }
}



void QR_rec(int dim, int init, Rmatrix A, Rmatrix Qstep, Rmatrix Q, Rmatrix aux)
{
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
    // copy from the first column
    for (j = init; j < dim; j++)
    {
        x[j-init] = A[j][init];
    }
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
    SpecialMult(dim,0,Qstep,Q,aux); // one new multiplication
    RmatCopy(dim,0,aux,Q);          // update final Q
    SpecialMult(dim,init,Qstep,A,aux);
    RmatCopy(dim,init,aux,A);
    free(x);
    QR_rec(dim,init+1,A,Qstep,Q,aux);
}



void QRdecomp(unsigned int dim, Rmatrix A, Rmatrix Q)
{
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
    SpecialMult(dim,0,Q,A,aux);
    RmatCopy(dim,0,aux,A);
    QR_rec(dim,1,A,Qstep,Q,aux);
    // Final step - transpose
    RmatCopy(dim,0,Q,aux);
    Transpose(dim,aux,Q);
    for (i = 0; i < dim; i++)
    {
        free(aux[i]);
        free(Qstep[i]);
    }
    free(aux);
    free(Qstep);
}



int main()
{
    int
        i,
        j,
        N;

    Rmatrix
        A,
        Q;

    N = 500;

    A = rmatDef(N,N);
    Q = rmatDef(N,N);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i][j] = ((i+1) % 3) * ((j+1) % 7 - 3) * 0.17827;
        }
    }

    QRdecomp(N,A,Q);

    /*
    A[0][0] = 12;
    A[1][0] = 6;
    A[2][0] = -4;
    A[0][1] = -51;
    A[1][1] = 167;
    A[2][1] = 24;
    A[0][2] = 4;
    A[1][2] = -68;
    A[2][2] = -41;
    QRdecomp(3,A,Q);
    for (i = 0; i < 3; i++)
    {
        printf("\n");
        for (j = 0; j < 3; j++)
        {
            printf(" %7.4lf",Q[i][j]);
        }
    }
    */

    printf("\n\n");
    return 0;
}
