#ifndef _HEigenState_h
#define _HEigenState_h

#include <mkl.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "HMatrix.h"



double complex carrDot(int n, Carray v1, Carray v2)
{

/** Convetional L2 scalar product for complex vectors **/

    int i;

    double complex z = 0;

    for (i = 0; i < n; i++) z = z + conj(v1[i]) * v2[i];

    return z;
}



double carrMod(int n, Carray v)
{

/** Conventional L2 modulus of complex vectors **/

    int i;

    double mod = 0;

    for (i = 0; i < n; i++)
    {
        mod = mod + creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }

    return sqrt(mod);
}



void matmul(int n, Iarray rows, Iarray cols, Carray vals,
            Carray vin, Carray vout)
{

/** Parallelized routine to perform matrix-vector multiplication from
    a sparse matrix structure given by 'rows', 'cols' and 'vals'  **/

    int
        i,
        j,
        threadId,
        nthreads;

    double complex
        z;

#pragma omp parallel private(i,j,z,threadId,nthreads)
    {
        threadId = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        for (i = threadId; i < n; i += nthreads)
        {
            z = 0;
            for (j = rows[i]; j < rows[i+1]; j++)
            {
                z = z + vals[j] * vin[cols[j]];
            }
            vout[i] = z;
        }
    }
}



int lanczos(int lm, int nc, HConfMat H, Carray diag, Carray offdiag,
            Cmatrix lvec)
{

/** Improved lanczos iterations with reorthogonalization for the vector
    of coefficients of the many-body state represented in configuration
    basis. Lanczos procedure give a (real)tridiagonal whose the spectra
    approximate the original one, as better as more iterations are done.
    It work, however, just with a routine to apply the original  matrix
    to any vector, what can save memory.

    OUTPUT PARAMETERS :
        lvec - Lanczos vectors used to convert eigenvectors of the
               resulting tridiagonal system back to the original one
        diag - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix

    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred) **/

    int
        i,
        j,
        k,
        threadId,
        nthreads;

    double
        tol,
        maxCheck;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

    printf("\n\nCALLING LANCZOS ITERATIONS\n\n");
    printf("-- Progress");

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    //applyHconf_omp(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,IF,
    //               lvec[0],Ho,Hint,HC);
    matmul(nc,H->rows,H->cols,H->vals,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrMod(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol) return i;

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        // apply Hamiltonian in a vector from configurational basis
        matmul(nc,H->rows,H->cols,H->vals,lvec[i+1],HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - offdiag[i] * lvec[i][j];
        }

        diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - diag[i+1]*lvec[i+1][j];
        }

        // Additional re-orthogonalization procedure. The main process
        // is parallelized because depending on the number  of Lanczos
        // iterations it may be the most demanding time,  beating even
        // the time to apply the Hamiltonian
        for (j = 0; j < i + 2; j++) ortho[j] = carrDot(nc, lvec[j], HC);

        #pragma omp parallel private(j,k,threadId,nthreads,hc_update)
        {
            threadId = omp_get_thread_num();
            nthreads = omp_get_num_threads();

            for (j = threadId; j < nc; j += nthreads)
            {
                hc_update = HC[j];
                for (k = 0; k < i + 2; k++)
                {
                    hc_update = hc_update - lvec[k][j] * ortho[k];
                }
                HC[j] = hc_update;
            }
        }

        printf("\n  %3d/%d",i+1,lm);
    }

    printf("\n===========   FINISHED LANCZOS ITERATIONS\n");

    free(ortho);
    free(HC);

    return lm;
}



double ground(int Niter, int nc, HConfMat H, Carray C)
{

/** Find the lowest Eigenvalue using Lanczos tridiagonal decomposition
    for the hamiltonian in configurational space, with orbitals fixed.
    Use up to Niter(unless the a breakdown occur) in Lanczos method to
    obtain a basis-fixed ground state approximation  of  the truncated
    configuration space.

    INPUT/OUTPUT : C (end up as eigenvector approximation)

    RETURN : Lowest eigenvalue found **/



    int
        i,
        k,
        j,
        predictedIter;

    double
        sentinel;

    Rarray
        d,
        e,
        eigvec;

    // Elements of tridiagonal lanczos matrix
    Carray
        diag,
        offdiag;

    Cmatrix
        lvec;

    // variables to call LAPACK routine. eigvec matrix is stored in
    // a vector in row major order
    d = rarrDef(Niter);
    e = rarrDef(Niter);
    eigvec = rarrDef(Niter * Niter);

    // output of lanczos iterative method
    diag = carrDef(Niter);
    offdiag = carrDef(Niter);

    // Lanczos Vectors (organize by rows of the following matrix)
    lvec = cmatDef(Niter,nc);

    // initiate date to call lanczos
    offdiag[Niter-1] = 0;
    for (i = 0; i < nc; i++)  lvec[0][i] = C[i];

    // Call Lanczos to setup tridiagonal matrix and lanczos vectors
    predictedIter = Niter;

    Niter = lanczos(Niter,nc,H,diag,offdiag,lvec);

    if (Niter < predictedIter)
    {
        printf("\n\nWARNING : lanczos iterations exit before expected.");
        printf(" Supposed to do %d but did %d\n\n",predictedIter,Niter);
    }

    // Transfer data to use lapack routine
    for (k = 0; k < Niter; k++)
    {
        d[k] = creal(diag[k]);    // Supposed to be real
        e[k] = creal(offdiag[k]); // Supposed to be real
        for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
    }

    k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
    if (k != 0)
    {
        printf("\n\nERROR IN DIAGONALIZATION\n\n");
        if (k < 0)
        {
            printf("Illegal value in parameter %d\n\n",-k);
        }
        else
        {
            printf("Algorithm failed to converge\n\n");
        }
        exit(EXIT_FAILURE);
    }

    sentinel = 1E10;
    // Get Index of smallest eigenvalue, keep it on j
    for (k = 0; k < Niter; k++)
    {
        if (sentinel > d[k]) { sentinel = d[k];   j = k; }
    }

    // Update C with the coefficients of ground state
    for (i = 0; i < nc; i++)
    {
        C[i] = 0;
        for (k = 0; k < Niter; k++) C[i] += lvec[k][i] * eigvec[k * Niter + j];
    }

    free(d);
    free(e);
    free(eigvec);
    free(diag);
    free(offdiag);
    for (i = 0; i < predictedIter; i++) free(lvec[i]);
    free(lvec);

    return sentinel;
}

#endif
