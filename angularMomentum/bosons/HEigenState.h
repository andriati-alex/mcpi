#ifndef _HEigenState_h
#define _HEigenState_h

#include <mkl.h>
#include "HMatrix.h"



void sepline()
{

/** print in current screen a separation line **/

    printf("\n=======================================");
    printf("=======================================\n");
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
        printf("\n\nERROR: NULL file in carr_inline routine\n\n");
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < M; j ++)
    {
        real = creal(v[j]);
        imag = cimag(v[j]);

        if (imag >= 0) fprintf(f, "(%.15E+%.15Ej) ", real, imag);
        else           fprintf(f, "(%.15E%.15Ej) ", real, imag);
    }

    fprintf(f,"\n");
}



void initGuess(unsigned int nc, Carray C)
{
    int
        i;

    double
        sum,
        real,
        imag;

    sum = 0;
    for (i = 0; i < nc; i++)
    {
        real = ((i % 4) * (3.1864 - 7 % (i+1)));
        imag = ((i * i) % 13 - 15.34576 * (i % 3));
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



int lanczos_mat(int lm, int nc, HConfMat H, Carray diag, Carray offdiag,
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



int lanczos_func(int lm, int nc, int lmax, Iarray * ht, Carray Ho, double g,
                 Carray diag, Carray offdiag, Cmatrix lvec)
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
    actH(lmax,nc,ht,Ho,g,lvec[0],HC);
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
        // matmul(nc,H->rows,H->cols,H->vals,lvec[i+1],HC);
        actH(lmax,nc,ht,Ho,g,lvec[i+1],HC);

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



double ground(int Niter, int Npar, int lmax, int total_mom, Carray C,
              Carray Ho, double g)
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
        nc,
        predictedIter;

    double
        GSenergy;

    Iarray
        * ht;

    Rarray
        d,
        e,
        eigvec;

    Carray
        diag,
        offdiag;

    Cmatrix
        lvec;

    HConfMat
        H;

    printf("\n\n * EVALUATING LOWEST ENERGY STATE WITH LANCZOS\n");
    // Configurational basis setup
    nc = BFixedMom_mcsize(Npar,lmax,total_mom);
    ht = BAssembleHT(Npar,lmax,total_mom,nc);
    H = assembleH(Npar,lmax,nc,ht,Ho,g);

    // variables to call LAPACK routine. eigvec matrix is stored in
    // a vector in row major order
    d = rarrDef(Niter);
    e = rarrDef(Niter);
    eigvec = rarrDef(Niter * Niter);
    // output of lanczos iterative method - Tridiagonal decomposition
    diag = carrDef(Niter);
    offdiag = carrDef(Niter);
    // Lanczos Vectors (organized by rows of the following matrix)
    lvec = cmatDef(Niter,nc);
    // initiate date to call lanczos. The first vector is the input guess
    offdiag[Niter-1] = 0;
    for (i = 0; i < nc; i++)  lvec[0][i] = C[i];



    /***   CALL LANCZOS   ***/
    predictedIter = Niter;
    if (H == NULL)
    {
        printf("   Without the matrix set up");
        Niter = lanczos_func(Niter,nc,lmax,ht,Ho,g,diag,offdiag,lvec);
    }
    else
    {
        printf("   With the sparse matrix set up");
        Niter = lanczos_mat(Niter,nc,H,diag,offdiag,lvec);
    }



    // Assert if some breakdown occur evaluating Lanczos method
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



    /***   CALL LAPACK FOR TRIDIAGONAL LANCZOS OUTPUT   ***/
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



    GSenergy = d[0];
    j = 0;
    // Update C with the coefficients of ground state
    for (i = 0; i < nc; i++)
    {
        C[i] = 0;
        for (k = 0; k < Niter; k++)
        {
            C[i] = C[i] + lvec[k][i] * eigvec[k * Niter + j];
        }
    }

    free(d);
    free(e);
    free(eigvec);
    free(diag);
    free(offdiag);
    // free matrices
    for (i = 0; i < predictedIter; i++) free(lvec[i]);
    free(lvec);
    for (i = 0; i < nc; i++) free(ht[i]);
    free(ht);
    // If the Hamiltonian matrix could be allocated free it
    if (H != NULL) freeHmat(H);

    return GSenergy;
}



void groundScanning(int Niter, int Npar, int lmax, Iarray total_mom,
                    int n_cases, Carray Ho, double g, char fname [])
{
    int
        i,
        nc;

    double
        E0;

    Carray
        C;

    FILE
        * out_file;

    out_file = fopen(fname,"w");
    if (out_file == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n",fname);
        exit(EXIT_FAILURE);
    }

    printf("\n\n\n");
    printf("\t*******************************************************\n");
    printf("\t*                                                     *\n");
    printf("\t*  ANGULAR MOMENTUM SCANNING OF LOWEST ENERGY STATES  *\n");
    printf("\t*                                                     *\n");
    printf("\t*******************************************************\n");

    fprintf(out_file,"# momentum x ground_energy x Many-Body-State\n");

    for (i = 0; i < n_cases; i++)
    {
        printf("\n\n\n\nCOMPUTING GROUND STATE %d/%d ",i+1,n_cases);
        printf("WITH L = %d",total_mom[i]);
        sepline();

        nc = BFixedMom_mcsize(Npar,lmax,total_mom[i]);
        C = carrDef(nc);
        initGuess(nc,C);
        E0 = ground(Niter,Npar,lmax,total_mom[i],C,Ho,g);
        fprintf(out_file,"(%.1E+0.0j) ",(double)total_mom[i]);
        fprintf(out_file,"(%.5E+0.0j) ",E0/Npar);
        carr_inline(out_file,nc,C);
        free(C);
    }

    fclose(out_file);
}

#endif
