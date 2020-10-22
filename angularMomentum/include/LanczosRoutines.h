#ifndef _LanczosRoutines_h
#define _LanczosRoutines_h

#include <mkl.h>
#include "Hamiltonian.h"



/** DIFFERENT 'FLAVOURS' OF LANCZOS ALGORITHM
    =========================================
    The routines here perform the Lanczos iterations to  obtain  the
    tridiagonal decomposition, which shall be furter serves as input
    to LAPACK routine for diagonalization of  symmetric  tridiagonal
    system. Thus the main output parameters from Lanczos method  are
    the diagonal ('diag') and off-diagonal ('offdiag') arrays  which
    contain the values of the output tridiagonal matrix.

    The routines follow a stop condition criteria 'stopLanczos'. The
    condition depends on the relative variation of the lowest eigval
    produced after some iterations. The relative variation parameter
    is defined in 'DataTypesDefinition.h' as LNCZS_STOP_TOL

    In all cases,  the called 'Lanczos vectors' (lvec) are stored to
    further compute the many-body ground state through coeff. in the
    configurational basis.  Without  them  it  would  be possible to
    obtain the energy  only  and  without  the  re-orthogonalization
    procedure, which improve accuracy due to numerical fluctuations.  **/



int stopLanczos(int Niter, Carray diag, Carray offdiag, double * previousE,
                int supp)
{

/** IN BETWEEN LANCZOS INTERACTIONS, PROVIDE A CONVERGENCE CRITERION TO STOP
    Compute the ground state energy with the outcome of Lanczos iterations
    at intermediate iterations to analyze if it is worth to continue or if
    a relevant accuracy was already achieved                           **/

    int
        j,
        k,
        shallStop;

    double
        err;

    Rarray
        d,
        e,
        eigvec;

    shallStop = 0; // boolean with information to stop iterations or not

    // variables to call LAPACK routine. 'eigvec' matrix is  stored
    // in row major order, where the eigenvectors are along columns
    d = rarrDef(Niter);
    e = rarrDef(Niter);
    eigvec = rarrDef(Niter*Niter);

    // Transfer data to use lapack routine
    for (k = 0; k < Niter; k++)
    {
        d[k] = creal(diag[k]);    // Supposed to be real
        e[k] = creal(offdiag[k]); // Supposed to be real
        for (j = 0; j < Niter; j++) eigvec[k*Niter+j] = 0;
    }
    e[Niter-1] = 0;

    /***   CALL LAPACK FOR TRIDIAGONAL LANCZOS OUTPUT   ***/
    k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
    if (k != 0)  LAPACK_PROBLEM(k,"stopLanczos");

    // compute relative error of the current lowest eigenvalue approx.
    // If the eigenvalue is close enough to zero (up to pre-defined
    // tolerance by LNCZS_STOP_TOL) compute just the absolute error
    // to avoid zero-division
    if (fabs(d[0]) > LNCZS_STOP_TOL)
    {
        err = fabs(creal(offdiag[Niter-2])*eigvec[(Niter-1)*Niter]/d[0]);
    }
    else
    {
        err = fabs(creal(offdiag[Niter-2])*eigvec[(Niter-1)*Niter]);
    }

    if (!supp) printf("\nRelativeError : %.7lf",err);

    // if (1.0 - d[0]/(*previousE) < LNCZS_STOP_TOL) shallStop = 1;
    if (err < LNCZS_STOP_TOL) shallStop = 1;
    *previousE = d[0];
    free(d);
    free(e);
    free(eigvec);
    return shallStop;
}



void updateLvec(int nc, int kp, Cmatrix lvec, Rmatrix Q)
{

/** UPDATE LANCZOS VECTORS ACCORDING TO IMPLICIT RESTART TECHNIQUE **/

    int
        m,
        i,
        j;
    double complex
        sum;
    Carray
        aux;

    aux = carrDef(kp);

    for (m = 0; m < nc; m++)
    {
        for (i = 0; i < kp; i++)
        {
            sum = 0;
            for (j = 0; j < kp; j++)
            {
                sum = sum + lvec[j][m] * Q[j][i];
            }
            aux[i] = sum;
        }
        for (i = 0; i < kp; i++) lvec[i][m] = aux[i];
    }
    free(aux);
}



void updateTrid(int kp, Carray diag, Carray offdiag, Rmatrix Q)
{

/** UPDATE TRIDIAGONAL MATRIX TO RESTART LANCZOS **/

    int
        i,
        j,
        k;
    Rmatrix
        A,
        QT,
        aux;

    // p = Niter/2; k = Niter - p;
    k = kp - kp/2;
    A = rmatDef(kp,kp);
    QT = rmatDef(kp,kp);
    aux = rmatDef(kp,kp);

    // Set a full matrix to compute the matrix mult.
    for (i = 0; i < kp; i++) { for (j = 0; j < kp; j++) A[i][j] = 0; }
    for (i = 0; i < kp - 1; i++)
    {
        A[i][i] = creal(diag[i]);
        A[i][i+1] = creal(offdiag[i]);
        A[i+1][i] = creal(offdiag[i]);
    }
    A[i][i] = creal(diag[i]);

    rmatTranspose(kp,Q,QT);
    rmatBlockMult(kp,0,A,Q,aux);
    rmatBlockMult(kp,0,QT,aux,A);

    // reconstruction of tridiagonal part up to iteration  'k'  from
    // the submatrix block of dimension 'k' in matrix A note that an
    // extra element is set in 'offdiag' that is outside  the  block
    // This elements is required from the algorithm proposed
    for (i = 0; i < k; i++)
    {
        diag[i] = A[i][i];
        offdiag[i] = A[i][i+1];
    }

/*  UNCOMMENT THIS SECTION TO PRINT NEW TRIDIAGONAL BLOCK MATRIX TO RESTART
    printf("\n\nTransformed tridiagonal Matrix\n\n");
    for (i = 0; i < kp; i++)
    {
        printf("\n");
        for (j = 0; j < kp; j++) printf("%10.1E",A[i][j]);
    }
    printf("\n\n");
*/

    rmatFree(kp,A);
    rmatFree(kp,QT);
    rmatFree(kp,aux);
}



void impRestart(int nc, int kp, Cmatrix lvec, Carray diag, Carray offdiag,
                double beta, Carray f)
{

/** UPDATE TRIDIAGONAL DECOMPOSITION AND LANCZOS VECTORS OF THE FIRST  'k'
    ITERATIONS GIVEN 'kp' ITERATIONS PERFORMED TO RESTART
    ======================================================================
    [1]
    "An implicitly restarted Lanczos method for large symmetric eigenvalue
    problems", D. Calvetti, L. Reichel & D.C. Sorensen. Vol 2, 1994.
    ----------------------------------------------------------------------
    [2]
    "Lecture Notes on Solving Large Scale Eigenvalue Problems", P. Arbens.
    Computer Science Department, ETH Zurich, 2016.
    Specially see algorithm 11.3 for the recursive QR decompositions using
    the shifts
    ----------------------------------------------------------------------
    Here is selected p = kp / 2 shifts to restart  the  Lanczos  algorithm
    from iteration k = kp - p (in this way, p <= k). From  the  references
    we must perform a QR-decomposition RECURSIVELY using  the  shifts  and
    the tridiagonal matrix from the 'kp' iterations. For each matrix Q get
    from the shifts we must multiply by all the previously  'Q's,  process
    that ends up in the final matrix Q to restart the Lanczos algorithm.
    INPUT PARAMETERS:
        nc - size of Lanczos vectors(multiconfigurational space dimention)
        kp - number of Lanczos iterations done (kp = k + p)
        lvec - Lanczos vectors along rows
        diag and offdiag - Lanczos tridiagonal decompostion
        beta - the 'beta' coefficient associated to (kp+1) Lanczos vector
        f - the (kp+1) Lanczos vector
    OUTPUT PARAMETERS:
        diag and offdiag - First 'k' elements updated (restarted)
        lvec - First 'k' vectors updated (restarted)
        'f' - New initial vector to continue from iteration 'k'        **/

    int
        i,
        j,
        k,
        p;

    Rarray
        d,
        e,
        eigvec,
        shifts;

    Rmatrix
        A,
        Q,
        Qk,
        aux;

    p = kp / 2;             // number of shifts used (highest eigenvalues)
    d = rarrDef(kp);
    e = rarrDef(kp);
    eigvec = rarrDef(kp*kp);
    shifts = rarrDef(p);    // shifts to be applied
    A = rmatDef(kp,kp);     // full tridiagonal matrix
    Q = rmatDef(kp,kp);     // Matrix to transform into the restarted problem
    Qk = rmatDef(kp,kp);    // Matrix for each QR-decomposition
    aux = rmatDef(kp,kp);   // keep A matrix values

    // Initialize the matrix to update(restart) Lanczos vectors
    // As it is the multiplication of each  QR-decompostion  of
    // each shift, start with the identity matrix
    for (i = 0; i < kp; i++)
    {
        for (j = 0; j < kp; j++) Q[i][j] = 0;
        Q[i][i] = 1.0;
    }

    // Transfer data to use lapack routine
    for (k = 0; k < kp; k++)
    {
        d[k] = creal(diag[k]);    // Supposed to be real
        e[k] = creal(offdiag[k]); // Supposed to be real
        for (j = 0; j < kp; j++) eigvec[k*kp+j] = 0;
    }
    e[kp-1] = 0;

    /***   CALL LAPACK FOR TRIDIAGONAL LANCZOS OUTPUT   ***/
    k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',kp,d,e,eigvec,kp);
    if (k != 0)  LAPACK_PROBLEM(k,"impRestart");

    /*** DEFINE THE SHIFTS BASED ON LARGEST EIGENVALUES ***/
    k = kp - p; // Iteration to continue after the restart
    for (i = 0; i < p; i++) shifts[i] = d[i+k];

    // Step 1. initialize tridiagonal matrix for decomp
    for (i = 0; i < kp; i++) { for (j = 0; j < kp; j++) A[i][j] = 0; }
    for (i = 0; i < kp - 1; i++)
    {
        A[i][i] = creal(diag[i]);
        A[i][i+1] = creal(offdiag[i]);
        A[i+1][i] = creal(offdiag[i]);
    }
    A[i][i] = creal(diag[i]);
    for (k = 0; k < p; k++)
    {
        rmatCopy(kp,0,A,aux);           // safety copy
        for (i = 0; i < kp; i++) A[i][i] = A[i][i] - shifts[k];
        // Step 2. Perform QR decomposition for the current shift 'k'
        QRdecomp(kp,A,Qk);
        // Step 3. Matrix multiplication by the very new Qk obtained
        rmatBlockMult(kp,0,Q,Qk,A);     // uses A as auxiliar output
        rmatCopy(kp,0,A,Q);             // Q holds all matrix multiplications
        // step 4. Recursively update A for next QR-decomposition. Ref [2]
        rmatBlockMult(kp,0,aux,Qk,A);   // 'aux' has the old 'A'
        rmatTranspose(kp,Qk,aux);
        rmatBlockMult(kp,0,aux,A,Qk);   // use Qk as auxiliar output
        rmatCopy(kp,0,Qk,A);
    }

    // assert the matrix Q matrix from all QR-decomposition is unitary
    rmatTranspose(kp,Q,Qk);
    rmatBlockMult(kp,0,Qk,Q,aux);
    if (!numericalIdentityMatrix(kp,aux))
    {
        printf("\n\nERROR : Matrix Q in Lanczos restart is not unitary ");
        printf("according to numerical precision\n\n");
        exit(EXIT_FAILURE);
    }

    // UPDATE/RESTART tridiagonal system and Lanczos vectors
    updateLvec(nc,kp,lvec,Q);
    updateTrid(kp,diag,offdiag,Q);

    // update initial vector to continue from iteration 'k'
    // equivalent to f^(plus) after Eq. (2.12)  in Ref [1]
    k = kp - p;
    for (i = 0; i < nc; i++)
    {
        f[i] = offdiag[k-1] * lvec[k][i] + f[i] * beta * Q[kp-1][k-1];
    }

/*  UNCOMMENT THIS SECTION TO PRINT Q MATRIX
    printf("\n\nResult of all QR-decomp. multiplications\n\n");
    for (i = 0; i < kp; i++)
    {
        printf("\n");
        for (j = 0; j < kp; j++) printf("%10.1E",Q[i][j]);
    }
    printf("\n\n");
*/

    free(d);
    free(e);
    free(eigvec);
    free(shifts);
    rmatFree(kp,A);
    rmatFree(kp,Q);
    rmatFree(kp,Qk);
    rmatFree(kp,aux);
}



int LNCZS_HMAT(int lm, int Npar, int nc, HConfMat H, Carray diag,
               Carray offdiag, Cmatrix lvec, int supp)
{

/** IMPROVED LANCZOS ITERATIONS WITH REORTHOGONALIZATION USING H. MATRIX
    OUTPUT PARAMETERS :
        lvec    - Lanczos vectors used to convert eigenvectors  of  the
                  resulting tridiagonal system back to the original one
        diag    - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix
    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred)
    OTHER PARAMETERS ARE REQUIRED TO APPLY HAMILTONIAN IN CONFIG. BASIS **/

    int
        i,
        j,
        k,
        Niter,
        threadId,
        nthreads;
    double
        tol,
        maxCheck,
        energy,
        extraB;
    double complex
        hc_update;
    Carray
        HC,
        ortho;

    if (!supp)
    {
        printf("\nLANCZOS ITERATIONS\n");
        printf(" * Progress");
    }

    Niter = 0;

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    matmul(nc,H->rows,H->cols,H->vals,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    energy = creal(diag[0]);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
            printf("WHICH USES HAMILTONIAN MATRIX\n\n");
            return i;
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        // apply Hamiltonian in a vector from configurational basis
        matmul(nc,H->rows,H->cols,H->vals,lvec[i+1],HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - offdiag[i] * lvec[i][j];
        }

        diag[i + 1] = carrDot(nc,lvec[i+1],HC);

        // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
        if ((i+1) % 10 == 0)
        {
            if (stopLanczos(i+1,diag,offdiag,&energy,supp))
            {
                if (!supp)
                {
                    printf(" | Energy per particle : %.7lf",energy/Npar);
                    printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                }
                free(HC);
                free(ortho);
                return i+1;
            }
            if (!supp) printf(" | Energy per particle : %.7lf",energy/Npar);
        }
        // -------------------------------------------------------

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

        Niter++;
        if (!supp) printf("\n  %3d/%d",i+1,lm);
    }

    /************************************************************************
    ************   RESTART ITERATIONS TO IMPROVE LANCZOS SPACE   ************
    *************************************************************************/

    while (!stopLanczos(lm,diag,offdiag,&energy,1))
    {
        if (!supp)
        {
            printf(" | Energy per particle : %.7lf",energy/Npar);
            printf("\n------------- IMPLICIT RESTART");
        }

        // Implicit restart algorithm requires one extra Lanczos vector
        extraB = carrNorm(nc,HC);                           // beta k+p
        for (j = 0; j < nc; j++) HC[j] = HC[j] / extraB;    // Lvector : k+p+1

        // ===================================================================
        // Restart iterations in k = lm - lm/2
        impRestart(nc,lm,lvec,diag,offdiag,extraB,HC);
        for (i = lm - lm/2 - 1; i < lm - 1; i++)
        {
            offdiag[i] = carrNorm(nc,HC);

            if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);
            // If method break return number of iterations achieved
            if (creal(offdiag[i]) / maxCheck < tol)
            {
                free(HC);
                free(ortho);
                printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
                printf("FOR TWO SPECIES BOSONIC MIXTURE\n\n");
                return i;
            }

            for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

            // apply Hamiltonian in a vector from configurational basis
            matmul(nc,H->rows,H->cols,H->vals,lvec[i+1],HC);

            for (j = 0; j < nc; j++)
            {
                HC[j] = HC[j] - offdiag[i] * lvec[i][j];
            }

            diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

            // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
            /*
            if ((i+1) % 10 == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy))
                {
                    printf(" | Energy per particle : %.7lf",energy/Npar);
                    printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                    free(HC);
                    free(ortho);
                    return i+1;
                }
                printf(" | Energy per particle : %.7lf",energy/Npar);
            }
            */
            // -------------------------------------------------------

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
            Niter++;
            if (!supp) printf("\n  %3d/%d | %d",i+1,lm,Niter);
        }
    }

    if (!supp) printf("\n============= FINISHED LANCZOS ITERATIONS\n");

    free(ortho);
    free(HC);

    return lm;
}



int LNCZS_HACT(int lm, int nc, int lmax, Iarray * ht, Carray Ho, double g,
                 Carray diag, Carray offdiag, Cmatrix lvec)
{

/** IMPROVED LANCZOS ITERATIONS WITH REORTHOGONALIZATION WITHOUT H. MATRIX
    OUTPUT PARAMETERS :
        lvec    - Lanczos vectors used to convert eigenvectors  of  the
                  resulting tridiagonal system back to the original one
        diag    - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix
    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred)
    OTHER PARAMETERS ARE REQUIRED TO APPLY HAMILTONIAN IN CONFIG. BASIS **/

    int
        i,
        j,
        k,
        Npar,
        Niter,
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy,
        extraB;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

    printf("\nLANCZOS ITERATIONS\n");
    printf(" * Progress");

    Niter = 0;

    // Compute Number of particles
    Npar = 0;
    for (i = 0; i < 2 * lmax + 1; i++) Npar = Npar + ht[0][i];

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    actH(lmax,nc,ht,Ho,g,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    energy = creal(diag[0]);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
            printf("WITHOUT USING HAMILTONIAN MATRIX FOR 1 SPECIE\n\n");
            return i;
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        // apply Hamiltonian in a vector from configurational basis
        actH(lmax,nc,ht,Ho,g,lvec[i+1],HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - offdiag[i] * lvec[i][j];
        }

        diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

        // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
        if ((i+1) % 10 == 0)
        {
            if (stopLanczos(i+1,diag,offdiag,&energy,0))
            {
                printf(" | Energy per particle : %.7lf",energy/Npar);
                free(HC);
                free(ortho);
                printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                return i+1;
            }
            printf(" | Energy per particle : %.7lf",energy/Npar);
        }
        // -------------------------------------------------------

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

        Niter++;
        printf("\n  %3d/%d",i+1,lm);
    }

    /************************************************************************
    ************   RESTART ITERATIONS TO IMPROVE LANCZOS SPACE   ************
    *************************************************************************/

    while (!stopLanczos(lm,diag,offdiag,&energy,0))
    {
        printf(" | Energy per particle : %.7lf",energy/Npar);
        printf("\n------------- IMPLICIT RESTART");

        // Implicit restart algorithm requires one extra Lanczos vector
        extraB = carrNorm(nc,HC);                           // beta k+p
        for (j = 0; j < nc; j++) HC[j] = HC[j] / extraB;    // Lvector : k+p+1

        // ===================================================================
        // Restart iterations in k = lm - lm/2
        impRestart(nc,lm,lvec,diag,offdiag,extraB,HC);
        for (i = lm - lm/2 - 1; i < lm - 1; i++)
        {
            offdiag[i] = carrNorm(nc,HC);

            if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);
            // If method break return number of iterations achieved
            if (creal(offdiag[i]) / maxCheck < tol)
            {
                free(HC);
                free(ortho);
                printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
                printf("FOR TWO SPECIES BOSONIC MIXTURE\n\n");
                return i;
            }

            for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

            actH(lmax,nc,ht,Ho,g,lvec[i+1],HC);

            for (j = 0; j < nc; j++)
            {
                HC[j] = HC[j] - offdiag[i] * lvec[i][j];
            }

            diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

            // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
            /*
            if ((i+1) % 10 == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy,0))
                {
                    printf(" | Energy per particle : %.7lf",energy/Npar);
                    free(HC);
                    free(ortho);
                    printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                    return i+1;
                }
                printf(" | Energy per particle : %.7lf",energy/Npar);
            }
            */
            // -------------------------------------------------------

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
            Niter++;
            printf("\n  %3d/%d | %d",i+1,lm,Niter);
        }
    }

    printf("\n============= FINISHED LANCZOS ITERATIONS\n");

    free(ortho);
    free(HC);

    return lm;
}



int LNCZS_BBMIX(int lm, CompoundSpace MixSpace, Carray HoA, Carray HoB,
                double g [], Carray diag, Carray offdiag, Cmatrix lvec)
{

/** IMPROVED LANCZOS ITERATIONS WITH REORTHOGONALIZATION FOR MIXTURE OF BOSONS
    OUTPUT PARAMETERS :
        lvec    - Lanczos vectors used to convert eigenvectors  of  the
                  resulting tridiagonal system back to the original one
        diag    - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix
    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred)
    OTHER PARAMETERS ARE REQUIRED TO APPLY HAMILTONIAN IN CONFIG. BASIS **/

    int
        i,
        j,
        k,
        nc,
        Npar,
        Niter,
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy,
        extraB;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

    Npar = MixSpace->Na + MixSpace->Nb;
    nc = MixSpace->size;
    Niter = 0;

    printf("\nLANCZOS ITERATIONS\n");
    printf(" * Progress");

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    mixture_actH(MixSpace,HoA,HoB,g,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    energy = creal(diag[0]);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
            printf("FOR TWO SPECIES BOSONIC MIXTURE\n\n");
            return i;
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        mixture_actH(MixSpace,HoA,HoB,g,lvec[i+1],HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - offdiag[i] * lvec[i][j];
        }

        diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

        // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
        if ((i+1) % 10 == 0)
        {
            if (stopLanczos(i+1,diag,offdiag,&energy,0))
            {
                printf(" | Energy per particle : %.7lf",energy/Npar);
                free(HC);
                free(ortho);
                printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                return i+1;
            }
            printf(" | Energy per particle : %.7lf",energy/Npar);
        }
        // -------------------------------------------------------

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

        Niter++;
        printf("\n  %3d/%d",i+1,lm);
    }

    /************************************************************************
    ************   RESTART ITERATIONS TO IMPROVE LANCZOS SPACE   ************
    *************************************************************************/

    while (!stopLanczos(lm,diag,offdiag,&energy,0))
    {
        printf(" | Energy per particle : %.7lf",energy/Npar);
        printf("\n------------- IMPLICIT RESTART");

        // Implicit restart algorithm requires one extra Lanczos vector
        extraB = carrNorm(nc,HC);                           // beta k+p
        for (j = 0; j < nc; j++) HC[j] = HC[j] / extraB;    // Lvector : k+p+1

        // ===================================================================
        // Restart iterations in k = lm - lm/2
        impRestart(nc,lm,lvec,diag,offdiag,extraB,HC);
        for (i = lm - lm/2 - 1; i < lm - 1; i++)
        {
            offdiag[i] = carrNorm(nc,HC);

            if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);
            // If method break return number of iterations achieved
            if (creal(offdiag[i]) / maxCheck < tol)
            {
                free(HC);
                free(ortho);
                printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
                printf("FOR TWO SPECIES BOSONIC MIXTURE\n\n");
                return i;
            }

            for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

            mixture_actH(MixSpace,HoA,HoB,g,lvec[i+1],HC);

            for (j = 0; j < nc; j++)
            {
                HC[j] = HC[j] - offdiag[i] * lvec[i][j];
            }

            diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

            // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
            /*
            if ((i+1) % 10 == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy,0))
                {
                    printf(" | Energy per particle : %.7lf",energy/Npar);
                    free(HC);
                    free(ortho);
                    printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                    return i+1;
                }
                printf(" | Energy per particle : %.7lf",energy/Npar);
            }
            */
            // -------------------------------------------------------

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
            Niter++;
            printf("\n  %3d/%d | %d",i+1,lm,Niter);
        }
    }

    printf("\n============= FINISHED LANCZOS ITERATIONS\n");

    free(ortho);
    free(HC);

    return lm;
}



int LNCZS_BFMIX(int lm, BFCompoundSpace MixSpace, Carray HoB, Carray HoF,
                double g [], Carray diag, Carray offdiag, Cmatrix lvec)
{

/** IMPROVED LANCZOS ITERATIONS WITH REORTHOGONALIZATION FOR BOSE-FERMI MIX
    OUTPUT PARAMETERS :
        lvec    - Lanczos vectors used to convert eigenvectors  of  the
                  resulting tridiagonal system back to the original one
        diag    - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix
    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred)
    OTHER PARAMETERS ARE REQUIRED TO APPLY HAMILTONIAN IN CONFIG. BASIS **/

    int
        i,
        j,
        k,
        nc,
        Npar,
        Niter,
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy,
        extraB;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

    Niter = 0;
    nc = MixSpace->size;
    Npar = MixSpace->Nb + MixSpace->Nf;

    printf("\nLANCZOS ITERATIONS\n");
    printf(" * Progress");

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    bosefermi_actH(MixSpace,HoB,HoF,g,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    energy = creal(diag[0]);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
            printf("FOR TWO SPECIES BOSE-FERMI MIXTURE\n\n");
            return i;
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        bosefermi_actH(MixSpace,HoB,HoF,g,lvec[i+1],HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - offdiag[i] * lvec[i][j];
        }

        diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

        // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
        if ((i+1) % 10 == 0)
        {
            if (stopLanczos(i+1,diag,offdiag,&energy,0))
            {
                printf(" | Energy per particle : %.7lf",energy/Npar);
                free(HC);
                free(ortho);
                printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                return i+1;
            }
            printf(" | Energy per particle : %.7lf",energy/Npar);
        }
        // -------------------------------------------------------

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

        Niter++;
        printf("\n  %3d/%d",i+1,lm);
    }

    /************************************************************************
    ************   RESTART ITERATIONS TO IMPROVE LANCZOS SPACE   ************
    *************************************************************************/

    while (!stopLanczos(lm,diag,offdiag,&energy,0))
    {
        printf(" | Energy per particle : %.7lf",energy/Npar);
        printf("\n------------- IMPLICIT RESTART");

        // Implicit restart algorithm requires one extra Lanczos vector
        extraB = carrNorm(nc,HC);                           // beta k+p
        for (j = 0; j < nc; j++) HC[j] = HC[j] / extraB;    // Lvector : k+p+1

        // ===================================================================
        // Restart iterations in k = lm - lm/2
        impRestart(nc,lm,lvec,diag,offdiag,extraB,HC);
        for (i = lm - lm/2 - 1; i < lm - 1; i++)
        {
            offdiag[i] = carrNorm(nc,HC);

            if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);
            // If method break return number of iterations achieved
            if (creal(offdiag[i]) / maxCheck < tol)
            {
                free(HC);
                free(ortho);
                printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
                printf("FOR TWO SPECIES BOSE-FERMI MIXTURE\n\n");
                return i;
            }

            for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

            bosefermi_actH(MixSpace,HoB,HoF,g,lvec[i+1],HC);

            for (j = 0; j < nc; j++)
            {
                HC[j] = HC[j] - offdiag[i] * lvec[i][j];
            }

            diag[i + 1] = carrDot(nc,lvec[i + 1],HC);

            // ASSESS STOP CONDITION BASED ON ENERGY VARIATION -------
            /*
            if ((i+1) % 10 == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy,0))
                {
                    printf(" | Energy per particle : %.7lf",energy/Npar);
                    free(HC);
                    free(ortho);
                    printf("\n============= FINISHED LANCZOS ITERATIONS\n");
                    return i+1;
                }
                printf(" | Energy per particle : %.7lf",energy/Npar);
            }
            */
            // -------------------------------------------------------

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
            Niter++;
            printf("\n  %3d/%d | %d",i+1,lm,Niter);
        }
    }

    printf("\n===========   FINISHED LANCZOS ITERATIONS\n");

    free(ortho);
    free(HC);

    return lm;
}



void LNCZS_HMAT_TIME(int nc, HConfMat H, Carray diag, Carray offdiag,
                     Cmatrix lvec)
{

/** SHORT LANCZOS ITERATIONS FOR TIME PROPARATION
    OUTPUT PARAMETERS :
        lvec    - Lanczos vectors used to convert eigenvectors  of  the
                  resulting tridiagonal system back to the original one
        diag    - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix
    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred)
    OTHER PARAMETERS ARE REQUIRED TO APPLY HAMILTONIAN IN CONFIG. BASIS **/

    int
        i,
        j,
        k,
        lm,
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

    lm = LNCZS_TIME_IT;

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    matmul(nc,H->rows,H->cols,H->vals,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nERROR : BREAKDOWN OCCURRED IN LANCZOS(MATRIX) ");
            printf("METHOD FOR TIME INTEGRATION\n\n");
            exit(EXIT_FAILURE);
        }

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
    }

    free(ortho);
    free(HC);
}



void LNCZS_TIME(int nc, int lmax, Iarray * ht, Carray Ho, double g,
               Carray diag, Carray offdiag, Cmatrix lvec)
{

/** SHORT LANCZOS ITERATIONS FOR TIME PROPARATION
    OUTPUT PARAMETERS :
        lvec    - Lanczos vectors used to convert eigenvectors  of  the
                  resulting tridiagonal system back to the original one
        diag    - diagonal elements of tridiagonal symmetric matrix
        offdiag - symmetric elements of tridiagonal matrix
    RETURN :
        number of itertion done (just 'lm' if none breakdown occurred)
    OTHER PARAMETERS ARE REQUIRED TO APPLY HAMILTONIAN IN CONFIG. BASIS **/

    int
        i,
        j,
        k,
        lm,
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

    lm = LNCZS_TIME_IT;

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    actH(lmax,nc,ht,Ho,g,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nERROR : BREAKDOWN OCCURRED IN LANCZOS(MATRIX) ");
            printf("METHOD FOR TIME INTEGRATION\n\n");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        // apply Hamiltonian in a vector from configurational basis
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
    }

    free(ortho);
    free(HC);
}



void LNCZS_BBMIX_TIME(CompoundSpace MixSpace, Carray HoA, Carray HoB,
                      double g [], Carray diag, Carray offdiag, Cmatrix lvec)
{

    int
        i,
        j,
        k,
        nc,
        lm,
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

    lm = LNCZS_TIME_IT;

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    nc = MixSpace->size;
    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Initiate the method
    mixture_actH(MixSpace,HoA,HoB,g,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
            printf("FOR TWO SPECIES BOSONIC MIXTURE IN TIME EVOLUTION\n\n");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        mixture_actH(MixSpace,HoA,HoB,g,lvec[i+1],HC);

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
        for (j = 0; j < i + 2; j++) ortho[j] = carrDot(nc,lvec[j],HC);

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
    }

    free(ortho);
    free(HC);
}



void LNCZS_BFMIX_TIME(BFCompoundSpace MixSpace, Carray HoB, Carray HoF,
                      double g [], Carray diag, Carray offdiag, Cmatrix lvec)
{

    int
        i,
        j,
        k,
        nc,
        lm,
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

    lm = LNCZS_TIME_IT; // DEFAULT NUMBER OF ITERATIONS

    nc = MixSpace->size;
    HC = carrDef(nc);
    ortho = carrDef(lm);

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    // Initiate the method
    bosefermi_actH(MixSpace,HoB,HoF,g,lvec[0],HC);
    diag[0] = carrDot(nc,lvec[0],HC);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrNorm(nc,HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol)
        {
            free(HC);
            free(ortho);
            printf("\n\nWARNING : BREAKDOWN OCCURRED IN LANCZOS METHOD FOR ");
            printf("TWO SPECIES BOSE-FERMI MIXTURE IN TIME EVOLUTION\n\n");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];

        bosefermi_actH(MixSpace,HoB,HoF,g,lvec[i+1],HC);

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
        for (j = 0; j < i + 2; j++) ortho[j] = carrDot(nc,lvec[j],HC);

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
    }

    free(ortho);
    free(HC);
}



#endif
