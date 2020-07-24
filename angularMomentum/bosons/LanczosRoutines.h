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



int stopLanczos(int Niter, Carray diag, Carray offdiag, double * previousE)
{

/** IN BETWEEN LANCZOS INTERACTIONS, PROVIDE A CONVERGENCE CRITERION TO SETOP
    Compute the ground state energy with the outcome of Lanczos iterations
    at intermediate iterations to analyze if it is worth to continue or if
    a relevant accuracy was already achieved                              **/

    int
        j,
        k,
        shallStop;

    Rarray
        d,
        e,
        eigvec;

    shallStop = 0; // boolean with information to stop iterations or not

    // variables to call LAPACK routine. 'eigvec' matrix is stored in
    // row major order, where the eigenvectors are along columns
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
    if (k != 0)
    {
        printf("\n\nERROR IN STOP CHECKING OF LANCZOS ALGORITHM\n\n");
        if (k < 0)
        {
            printf("Illegal value in LAPACK_dstev parameter %d\n\n",-k);
        }
        else
        {
            printf("LAPACK algorithm failed to converge\n\n");
        }
        exit(EXIT_FAILURE);
    }

    if (1.0 - d[0]/(*previousE) < LNCZS_STOP_TOL) shallStop = 1;
    *previousE = d[0];
    free(d);
    free(e);
    free(eigvec);
    return shallStop;
}



int LNCZS_HMAT(int lm, int nc, HConfMat H, Carray diag, Carray offdiag,
            Cmatrix lvec)
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
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

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
    //applyHconf_omp(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,IF,
    //               lvec[0],Ho,Hint,HC);
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

        if (lm > 50)
        {
            if ((i+1) % (lm / 10) == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy))
                {
                    free(HC);
                    free(ortho);
                    printf("\nACHIEVED CONVERGENCE FOR GROUND STATE ENERGY\n");
                    return i+1;
                }
            }
        }

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

    printf("\n===========   MAX. LANCZOS ITERATIONS REACHED\n");

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
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

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
            printf("WITHOUT USING HAMILTONIAN MATRIX\n\n");
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

        if (lm > 50)
        {
            if ((i+1) % (lm / 10) == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy))
                {
                    free(HC);
                    free(ortho);
                    printf("\nACHIEVED CONVERGENCE FOR GROUND STATE ENERGY\n");
                    return i+1;
                }
            }
        }

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

    printf("\n===========   MAX. LANCZOS ITERATIONS REACHED\n");

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
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

    printf("\nLANCZOS ITERATIONS\n");
    printf(" * Progress");

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
        if (lm > 50)
        {
            if ((i+1) % (lm / 10) == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy))
                {
                    free(HC);
                    free(ortho);
                    printf("\n * ACHIEVED CONVERGENCE FOR GROUND STATE ENERGY\n");
                    return i+1;
                }
            }
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

        printf("\n  %3d/%d",i+1,lm);
    }

    printf("\n===========   FINISHED LANCZOS ITERATIONS\n");

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
        threadId,
        nthreads;

    double
        tol,
        maxCheck,
        energy;

    double complex
        hc_update;

    Carray
        HC,
        ortho;

    printf("\nLANCZOS ITERATIONS\n");
    printf(" * Progress");

    // Check for a source of breakdown in the algorithm to do not
    // divide by zero. Instead of zero use  a  tolerance (tol) to
    // avoid numerical instability
    maxCheck = 0;
    tol = 1E-15;

    nc = MixSpace->size;

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
            printf("FOR TWO SPECIES BOSONIC MIXTURE\n\n");
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
        if (lm > 50)
        {
            if ((i+1) % (lm / 10) == 0)
            {
                if (stopLanczos(i+1,diag,offdiag,&energy))
                {
                    free(HC);
                    free(ortho);
                    printf("\n * ACHIEVED CONVERGENCE FOR GROUND STATE ENERGY\n");
                    return i+1;
                }
            }
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

        printf("\n  %3d/%d",i+1,lm);
    }

    printf("\n===========   FINISHED LANCZOS ITERATIONS\n");

    free(ortho);
    free(HC);

    return lm;
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
        maxCheck,
        energy;

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
            printf("\n\nERROR : BREAKDOWN OCCURRED IN LANCZOS METHOD ");
            printf("FOR TIME INTEGRATION\n\n");
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



#endif
