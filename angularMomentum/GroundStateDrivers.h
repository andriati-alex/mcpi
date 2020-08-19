#ifndef _GroundStateDrivers_h
#define _GroundStateDrivers_h

#include "LanczosRoutines.h"



/** FUNCTIONS CALLED IN THE EXECUTABLE PROGRAMS
    ===========================================
    These functions provide an interface to read parameters
    call Lanczos algorithm and record output data.      **/



void parLine(char fname [], int line, int * N, int * lmax, int * L, double * g)
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
    fscanf(in_file,"%d",N);
    fscanf(in_file,"%d",lmax);
    fscanf(in_file,"%d",L);
    fscanf(in_file,"%lf",g);

    fclose(in_file);
}



void mixParLine(char fname [], int line, int * NA, int * lmaxA,
                int * NB, int * lmaxB, int * L, double g [])
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
    fscanf(in_file,"%d",NA);
    fscanf(in_file,"%d",lmaxA);
    fscanf(in_file,"%d",NB);
    fscanf(in_file,"%d",lmaxB);
    fscanf(in_file,"%d",L);
    fscanf(in_file,"%lf",&g[0]);
    fscanf(in_file,"%lf",&g[1]);
    fscanf(in_file,"%lf",&g[2]);

    fclose(in_file);
}



double IMAGTIME(int Nsteps, double dt, int Npar, int lmax, int total_mom,
                Carray C, Carray Ho, double g)
{

/** SHORT ITERATIVE LANCZOS(SIL) INTEGRATOR IN IMAGINARY TIME
    =========================================================
    INPUT/OUTPUT : C - coefficients at time  t = dt*Nsteps
    RETURN : Energy of the state at time t = dt*Nsteps **/

    int
        i,
        k,
        j,
        nc,
        step,
        Niter,
        predictedIter;

    double
        sum,
        GSenergy,
        previousE;

    Iarray
        * ht;

    Rarray
        d,
        e,
        aux,
        eigvec,
        Clanczos;

    Carray
        diag,
        offdiag;

    Cmatrix
        lvec;

    Niter = LNCZS_TIME_IT;

    printf("\n\n * EVALUATING LOWEST ENERGY STATE WITH IMAGINARY TIME\n");
    printf("\n T-step/Nsteps    time          E / Npar");
    sepline();
    // Configurational basis setup
    nc = BFixedMom_mcsize(Npar,lmax,total_mom);
    ht = BAssembleHT(Npar,lmax,total_mom,nc);

    // variables to call LAPACK routine. eigvec matrix is stored in
    // a vector in row major order
    d = rarrDef(Niter);
    e = rarrDef(Niter);
    eigvec = rarrDef(Niter*Niter);
    // output of lanczos iterative method - Tridiagonal decomposition
    diag = carrDef(Niter);
    offdiag = carrDef(Niter);
    offdiag[Niter-1] = 0;
    // Lanczos Vectors (organized by rows of the following matrix)
    // The first row is sometimes used as auxiliar output from H. action
    lvec = cmatDef(Niter,nc);
    aux = rarrDef(Niter);
    Clanczos = rarrDef(Niter);



    // Compute initial energy (GUESS)
    actH(lmax,nc,ht,Ho,g,C,lvec[0]);
    previousE = creal(carrDot(nc,C,lvec[0]));

    // Initiate time evolution
    for (step = 0; step < Nsteps; step++)
    {
        /***   CALL LANCZOS FOR TRIDIAGONAL DECOMPOSITION   ***/
        for (i = 0; i < nc; i++)  lvec[0][i] = C[i];
        LNCZS_TIME(nc,lmax,ht,Ho,g,diag,offdiag,lvec);

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
            printf("\n\nERROR IN DIAGONALIZATION IN IMAG. TIME\n\n");
            if (k < 0) printf("Illegal value in parameter %d\n\n",-k);
            else       printf("Algorithm failed to converge\n\n");
            exit(EXIT_FAILURE);
        }

        /*** USE TIME EVOLUTION IN LANCZOS SPACE AND THEN TRANSFORM BACK ***/
        Clanczos[0] = 1.0;
        for (k = 1; k < Niter; k++) Clanczos[k] = 0;

        for (k = 0; k < Niter; k++)
        {   // Solve in diagonal basis and for this apply eigvec trasformation
            aux[k] = eigvec[0*Niter + k] * Clanczos[0];
            // aux[k] = aux[k] * cexp(- I * d[k] * dt); // REAL TIME CASE
            aux[k] = aux[k] * exp(- d[k] * dt);
        }

        for (k = 0; k < Niter; k++)
        {   // Backward transformation from diagonal representation
            Clanczos[k] = 0;
            for (j = 0; j < Niter; j++)
            {
                Clanczos[k] += eigvec[k*Niter + j] * aux[j];
            }
        }

        for (i = 0; i < nc; i++)
        {   // Return from Lanczos space to configurational space
            C[i] = 0;
            for (j = 0; j < Niter; j++) C[i] += lvec[j][i] * Clanczos[j];
        }

        // renormalization
        sum = 0;
        for (i = 0; i < nc; i++)
        {
            sum = sum + creal(C[i])*creal(C[i]) + cimag(C[i])*cimag(C[i]);
        }
        for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);

        // update energy
        actH(lmax,nc,ht,Ho,g,C,lvec[0]);
        GSenergy = creal(carrDot(nc,C,lvec[0]));

        if ((step+1) % (Nsteps / 1000) == 0)
        {
            printf("\n %6d/%d   ",step+1,Nsteps);
            printf(" %10.6lf",(step+1)*dt);
            printf(" %16.9lf",GSenergy/Npar);
        }

        if ((step+1) % (Nsteps / 20) == 0)
        {
            if (1.0 - GSenergy/previousE < 1E-9)
            {
                // ACHIEVED CONVERGENCE CRITERIA
                free(d);
                free(e);
                free(eigvec);
                free(aux);
                free(diag);
                free(offdiag);
                free(Clanczos);
                // free matrices
                for (i = 0; i < Niter; i++) free(lvec[i]);
                free(lvec);
                for (i = 0; i < nc; i++) free(ht[i]);
                free(ht);

                printf("\n\nACHIEVED CONVERGENCE CRITERION");

                return GSenergy;
            }
            else previousE = GSenergy;
        }
    }



    free(d);
    free(e);
    free(eigvec);
    free(aux);
    free(diag);
    free(offdiag);
    free(Clanczos);
    // free matrices
    for (i = 0; i < Niter; i++) free(lvec[i]);
    free(lvec);
    for (i = 0; i < nc; i++) free(ht[i]);
    free(ht);

    return GSenergy;
}



double MIXTURE_IMAGTIME(int Nsteps, double dt, CompoundSpace S, Carray C,
                        Carray HoA, Carray HoB, double g [])
{

    int
        i,
        k,
        j,
        nc,
        Npar,
        step,
        Niter,
        predictedIter;

    double
        sum,
        GSenergy,
        previousE;

    Rarray
        d,
        e,
        aux,
        eigvec,
        Clanczos;

    Carray
        diag,
        offdiag;

    Cmatrix
        lvec;

    Niter = LNCZS_TIME_IT;
    nc = S->size;
    Npar = S->Na + S->Nb;

    printf("\n\nIMAGINARY TIME FOR BOSONIC MIXTURES\n");
    printf("\n T-step/Nsteps    time          E / Npar");
    sepline();

    // variables to call LAPACK routine. eigvec matrix is stored in
    // a vector in row major order
    d = rarrDef(Niter);
    e = rarrDef(Niter);
    eigvec = rarrDef(Niter*Niter);
    // output of lanczos iterative method - Tridiagonal decomposition
    diag = carrDef(Niter);
    offdiag = carrDef(Niter);
    offdiag[Niter-1] = 0;
    // Lanczos Vectors (organized by rows of the following matrix)
    // The first row is sometimes used as auxiliar output from H. action
    lvec = cmatDef(Niter,nc);
    aux = rarrDef(Niter);
    Clanczos = rarrDef(Niter);

    // Compute initial energy (GUESS)
    mixture_actH(S,HoA,HoB,g,C,lvec[0]);
    previousE = creal(carrDot(nc,C,lvec[0]));

    // Initiate time evolution
    for (step = 0; step < Nsteps; step++)
    {
        /***   CALL LANCZOS FOR TRIDIAGONAL DECOMPOSITION   ***/
        for (i = 0; i < nc; i++)  lvec[0][i] = C[i];
        LNCZS_BBMIX_TIME(Niter,S,HoA,HoB,g,diag,offdiag,lvec);

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
            printf("\n\nERROR IN DIAGONALIZATION IN IMAG. TIME\n\n");
            if (k < 0) printf("Illegal value in parameter %d\n\n",-k);
            else       printf("Algorithm failed to converge\n\n");
            exit(EXIT_FAILURE);
        }

        /*** USE TIME EVOLUTION IN LANCZOS SPACE AND THEN TRANSFORM BACK ***/
        Clanczos[0] = 1.0;
        for (k = 1; k < Niter; k++) Clanczos[k] = 0;

        for (k = 0; k < Niter; k++)
        {   // Solve in diagonal basis and for this apply eigvec trasformation
            aux[k] = eigvec[0*Niter + k] * Clanczos[0];
            // aux[k] = aux[k] * cexp(- I * d[k] * dt); // REAL TIME CASE
            aux[k] = aux[k] * exp(- d[k] * dt);
        }
        for (k = 0; k < Niter; k++)
        {   // Backward transformation from diagonal representation
            Clanczos[k] = 0;
            for (j = 0; j < Niter; j++)
            {
                Clanczos[k] += eigvec[k*Niter + j] * aux[j];
            }
        }
        for (i = 0; i < nc; i++)
        {   // Return from Lanczos space to configurational space
            C[i] = 0;
            for (j = 0; j < Niter; j++) C[i] += lvec[j][i] * Clanczos[j];
        }

        // renormalization
        sum = 0;
        for (i = 0; i < nc; i++)
        {
            sum = sum + creal(C[i])*creal(C[i]) + cimag(C[i])*cimag(C[i]);
        }
        for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);

        // update energy
        mixture_actH(S,HoA,HoB,g,C,lvec[0]);
        GSenergy = creal(carrDot(nc,C,lvec[0]));

        if ((step+1) % (Nsteps / 1000) == 0)
        {
            printf("\n %6d/%d   ",step+1,Nsteps);
            printf(" %10.6lf",(step+1)*dt);
            printf(" %16.9lf",GSenergy/Npar);
        }

        if ((step+1) % (Nsteps / 20) == 0)
        {
            if (1.0 - GSenergy/previousE < 1E-9)
            {
                // ACHIEVED CONVERGENCE CRITERIA
                free(d);
                free(e);
                free(eigvec);
                free(aux);
                free(diag);
                free(offdiag);
                free(Clanczos);
                // free matrices
                for (i = 0; i < Niter; i++) free(lvec[i]);
                free(lvec);

                printf("\n\nACHIEVED CONVERGENCE CRITERION FOR IMAG TIME\n");

                return GSenergy;
            }
            else previousE = GSenergy;
        }
    }

    free(d);
    free(e);
    free(eigvec);
    free(aux);
    free(diag);
    free(offdiag);
    free(Clanczos);
    // free matrices
    for (i = 0; i < Niter; i++) free(lvec[i]);
    free(lvec);

    return GSenergy;
}



double GROUND_STATE(int Niter, int Npar, int lmax, int total_mom, Carray C,
                    Carray Ho, double g)
{

/** Find the lowest Eigenvalue using Lanczos tridiagonal decomposition
    for the Hamiltonian in configurational space, with orbitals fixed.
    Use up to Niter unless achieve convergence criterion or the method
    breakdown. The first case is preferable.

    INPUT/OUTPUT : C input guess / output ground state coefficients
    RETURN : Lowest eigenvalue found                            **/

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

    printf("\n * EVALUATING LOWEST ENERGY STATE WITH LANCZOS\n");

    // empty system
    if (Npar == 0) return 0;

    // Configurational basis setup
    nc = BFixedMom_mcsize(Npar,lmax,total_mom);
    // Assert the config. space is not void
    if (nc == 0)
    {
        printf("\n\nDIAGONALIZATION PROCESS ABORTED : The parameters ");
        printf("requested, Npar = %d, lmax = %d, L = %d",Npar,lmax,total_mom);
        printf(" generated empty config. space.\n\n");
        exit(EXIT_FAILURE);
    }

    ht = BAssembleHT(Npar,lmax,total_mom,nc);
    H = boseH(Npar,lmax,nc,ht,Ho,g);
    if (H == NULL) printf("   Without set up (sparse) Hamiltonian matrix\n");
    else           printf("   Using (sparse) Hamiltonian matrix\n");

    if (nc == 1)
    {
        // In this case there must be only one config. because there
        // is only one single particle state, the one with l = 0
        diag = carrDef(1);
        C[0] = 1;
        matmul(nc,H->rows,H->cols,H->vals,C,diag);
        GSenergy = creal(carrDot(nc,C,diag));
        free(diag);
        free(ht[0]);
        free(ht);
        freeHmat(H);
        return GSenergy;
    }

    // USE RESTARTS IF THE MEMORY FOR LANCZOS VECTORS WILL
    // EXCEED THE TOLERANCE DEFINED BEFORE COMPLETING  THE
    // DEFAULT NUMBER OF LANCZOS ITERATIONS 'MAX_LNCZS_IT'
    // RESTART EVERY 20 ITERATIONS ACCUMULATED
    if (nc > MAX_LNCZS_IT && Niter < MAX_LNCZS_IT) Niter = 20;

    // VARIABLES TO CALL LAPACK ROUTINE.
    d = rarrDef(Niter);             // diagonal elements
    e = rarrDef(Niter);             // off-diagonal elements
    eigvec = rarrDef(Niter*Niter);  // eigenvectors along columns
    
    // OUTPUT OF LANCZOS ITERATIVE METHOD - TRIDIAGONAL DECOMPOSITION
    diag = carrDef(Niter);      // diagonal elements
    offdiag = carrDef(Niter);   // off-diagonal elements
    offdiag[Niter-1] = 0;

    // LANCZOS VECTORS (ORGANIZED BY ROWS OF THE FOLLOWING MATRIX)
    lvec = cmatDef(Niter,nc);
    // THE FIRST VECTOR IS THE INPUT GUESS
    for (i = 0; i < nc; i++) lvec[0][i] = C[i];

    /***********************************
     ***   CALL LANCZOS ITERATIONS   ***
     ***********************************/

    predictedIter = Niter;
    if (H == NULL) Niter = LNCZS_HACT(Niter,nc,lmax,ht,Ho,g,diag,offdiag,lvec);
    else           Niter = LNCZS_HMAT(Niter,Npar,nc,H,diag,offdiag,lvec);

    // TRANSFER DATA TO USE LAPACK ROUTINE
    for (k = 0; k < Niter; k++)
    {
        d[k] = creal(diag[k]);    // Supposed to be real
        e[k] = creal(offdiag[k]); // Supposed to be real
        for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
    }

    /**************************************************************
     ***   CALL LAPACK FOR TRIDIAGONAL MATRIX DIAGONALIZATION   ***
     **************************************************************/

    k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
    if (k != 0) LAPACK_PROBLEM(k,"GROUND_STATE");

    GSenergy = d[0]; // lowest energy eigenvalue

    // UPDATE C WITH THE COEFFICIENTS OF GROUND STATE
    j = 0;
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



double BOSEBOSE_GS(int Niter, CompoundSpace MixSpace, Carray C, Carray HoA,
                   Carray HoB, double g [])
{

/** COMPUTE GROUND STATE OF TWO SPECIES BOSONIC MIXTURE USING LANCZOS 
    INPUT/OUTPUT : C input guess / output ground  state  coefficients
    RETURN : Lowest eigenvalue found                              **/

    int
        i,
        k,
        j,
        nc,
        Npar,
        predictedIter;
    double
        GSenergy;
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

    printf("\n * EVALUATING GROUND STATE OF TWO SPECIES BOSONIC MIXTURE\n");

    Npar = MixSpace->Na + MixSpace->Nb;

    nc = MixSpace->size;
    if (nc == 1)
    {
        // In this case there must be only one config. because there
        // is only one single particle state for each species
        diag = carrDef(1);
        C[0] = 1;
        mixture_actH(MixSpace,HoA,HoB,g,C,diag);
        GSenergy = creal(C[0]*diag[0]);
        free(diag);
        return GSenergy;
    }

    // TRY TO ALLOCATE MATRIX. IF THE MEMORY REQUIRED EXCEED
    // THE TOLERANCE DEFINED RETURN NULL POINTER
    H = boseboseH(MixSpace,HoA,HoB,g);
    if (H == NULL) printf("   Without set (sparse) Hamiltonian matrix\n");
    else           printf("   Using (sparse) Hamiltonian matrix\n");

    // USE RESTARTS IF THE MEMORY FOR LANCZOS VECTORS WILL
    // EXCEED THE TOLERANCE DEFINED BEFORE COMPLETING  THE
    // DEFAULT NUMBER OF LANCZOS ITERATIONS 'MAX_LNCZS_IT'
    // RESTART EVERY 20 ITERATIONS ACCUMULATED
    if (nc > MAX_LNCZS_IT && Niter < MAX_LNCZS_IT) Niter = 20;

    // VARIABLES TO CALL LAPACK DIAGONALIZATION OF TRIDIAGONAL MATRIX
    d = rarrDef(Niter);             // diagonal elements
    e = rarrDef(Niter);             // off-diagonal elements
    eigvec = rarrDef(Niter*Niter);  // eigenvectors in row-major-order

    // OUTPUT OF LANCZOS ITERATIVE METHOD - TRIDIAGONAL DECOMPOSITION
    diag = carrDef(Niter);      // diagonal elements
    offdiag = carrDef(Niter);   // off-diagonal elements
    offdiag[Niter-1] = 0;

    // LANCZOS VECTORS (ORGANIZED BY ROWS OF THE FOLLOWING MATRIX)
    lvec = cmatDef(Niter,nc);

    // INITIATE DATE TO CALL LANCZOS. THE FIRST VECTOR IS THE INPUT GUESS
    for (i = 0; i < nc; i++) lvec[0][i] = C[i];

    /************************
     ***   CALL LANCZOS   ***
     ************************/

    predictedIter = Niter;
    if (H == NULL)
    {
        Niter = LNCZS_BBMIX(Niter,MixSpace,HoA,HoB,g,diag,offdiag,lvec);
    }
    else
    {
        Niter = LNCZS_HMAT(Niter,Npar,nc,H,diag,offdiag,lvec);
    }

    // TRANSFER DATA TO USE LAPACK ROUTINE
    for (k = 0; k < Niter; k++)
    {
        d[k] = creal(diag[k]);    // Supposed to be real
        e[k] = creal(offdiag[k]); // Supposed to be real
        for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
    }

    /**************************************************************
     ***   CALL LAPACK FOR TRIDIAGONAL MATRIX DIAGONALIZATION   ***
     **************************************************************/

    k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
    if (k != 0)  LAPACK_PROBLEM(k,"BOSEBOSE_GS");

    GSenergy = d[0]; // lowest eigenvalue

    // UPDATE C WITH THE COEFFICIENTS OF GROUND STATE
    j = 0;
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
    if (H != NULL) freeHmat(H);

    return GSenergy;
}



double BOSEFERMI_GS(int Niter, BFCompoundSpace MixSpace, Carray C, Carray HoB,
                    Carray HoF, double g [])
{

/** COMPUTE GROUND STATE OF BOSE-FERMI MIXTURE USING LANCZOS
    INPUT/OUTPUT : C input guess / output ground state coefficients
    RETURN : Lowest eigenvalue found                            **/

    int
        i,
        k,
        j,
        nc,
        Npar,
        predictedIter;
    double
        GSenergy;
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

    printf("\n * EVALUATING GROUND STATE OF BOSE-FERMI MIXTURE\n");

    Npar = MixSpace->Nf + MixSpace->Nb;

    nc = MixSpace->size;
    if (MixSpace->size == 1)
    {
        // In this case there must be only one config. because there
        // is only one single particle state for each species
        diag = carrDef(1);
        C[0] = 1;
        bosefermi_actH(MixSpace,HoB,HoF,g,C,diag);
        GSenergy = creal(C[0]*diag[0]);
        free(diag);
        return GSenergy;
    }

    // TRY TO ALLOCATE MATRIX. IF THE MEMORY REQUIRED EXCEED
    // THE TOLERANCE DEFINED RETURN NULL POINTER
    H = bosefermiH(MixSpace,HoB,HoF,g);
    if (H == NULL) printf("   Without set (sparse) Hamiltonian matrix\n");
    else           printf("   Using (sparse) Hamiltonian matrix\n");

    // USE RESTARTS IF THE MEMORY FOR LANCZOS VECTORS WILL
    // EXCEED THE TOLERANCE DEFINED BEFORE COMPLETING  THE
    // DEFAULT NUMBER OF LANCZOS ITERATIONS 'MAX_LNCZS_IT'
    // RESTART EVERY 20 ITERATIONS ACCUMULATED
    if (nc > MAX_LNCZS_IT && Niter < MAX_LNCZS_IT) Niter = 20;

    // VARIABLES TO CALL LAPACK ROUTINE
    d = rarrDef(Niter);
    e = rarrDef(Niter);
    eigvec = rarrDef(Niter * Niter);

    // OUTPUT OF LANCZOS ITERATIVE METHOD - TRIDIAGONAL DECOMPOSITION
    diag = carrDef(Niter);
    offdiag = carrDef(Niter);
    offdiag[Niter-1] = 0;

    // LANCZOS VECTORS (organized by rows of the following matrix)
    lvec = cmatDef(Niter,nc);
    // The first vector is the input guess
    for (i = 0; i < nc; i++) lvec[0][i] = C[i];

    /************************
     ***   CALL LANCZOS   ***
     ************************/

    predictedIter = Niter;
    if (H == NULL)
    {
        Niter = LNCZS_BFMIX(Niter,MixSpace,HoB,HoF,g,diag,offdiag,lvec);
    }
    else
    {
        Niter = LNCZS_HMAT(Niter,Npar,nc,H,diag,offdiag,lvec);
    }

    // TRANSFER DATA TO USE LAPACK ROUTINE
    for (k = 0; k < Niter; k++)
    {
        d[k] = creal(diag[k]);    // Supposed to be real
        e[k] = creal(offdiag[k]); // Supposed to be real
        for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
    }

    /**************************************************************
     ***   CALL LAPACK FOR TRIDIAGONAL MATRIX DIAGONALIZATION   ***
     **************************************************************/

    k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
    if (k != 0) LAPACK_PROBLEM(k,"BOSEFERMI_GS");

    GSenergy = d[0]; // lowest eigenvalue

    // UPDATE C WITH THE COEFFICIENTS OF GROUND STATE
    j = 0;
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
    if (H != NULL) freeHmat(H);

    return GSenergy;
}



void SCANNING(int n_cases, char prefix [], unsigned int full_output)
{
    int
        i,
        j,
        nc,
        Npar,
        lmax,
        lan_it,
        total_mom,
        coefMemory;

    double
        l,
        g,
        E0;

    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100];

    Carray
        C,
        Ho;

    FILE
        * e_file,
        * out_file;

    // set up name of file containing input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set up name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_energy.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Energy per particle | Num. of Particles | L | g");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of configurational
        // space in line 'i' of input file
        parLine(in_fname,i,&Npar,&lmax,&total_mom,&g);
        // number of configuration(nc) - config. space dimension
        nc = BFixedMom_mcsize(Npar,lmax,total_mom);
        printf("\n\n\nCOMPUTING GROUND STATE %d/%d | ",i+1,n_cases);
        printf("N = %d lmax = %d L = %d g = %.10lf\n",Npar,lmax,total_mom,g);
        printf("space dimension %d",nc);
        sepline();
        // set up (default)  number of Lanczos iterations
        // respecting the memory allowed (MEMORY_TOL) and
        // maximum number of iterations (MAX_LNCZS_IT)
        coefMemory = nc * sizeof(double complex);
        if (nc > MAX_LNCZS_IT)
        {
            lan_it = 10;
            while (lan_it*coefMemory < MEMORY_TOL && lan_it < MAX_LNCZS_IT)
            {
                lan_it++;
            }
        }
        else lan_it = nc; // Small config. space
        // set up one-body matrix elements
        Ho = carrDef(2*lmax+1);
        for (j = 0; j < 2*lmax+1; j++)
        {
            l = 2 * PI * (j-lmax);
            Ho[j] = 0.5*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        initGuess(nc,C);
        // call (main)routine to compute the ground state
        E0 = GROUND_STATE(lan_it,Npar,lmax,total_mom,C,Ho,g);
        printf("\nAverage energy per particle = %.10lf",E0/Npar);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"\n%.10E %d %d %.10E",E0/Npar,Npar,total_mom,g);
        if (full_output)
        {
            // set output file name
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            out_file = openFileWrite(out_fname);
            fprintf(out_file,"(%.10E+0.0j)",E0/Npar); // energy per particle
            carrAppend(out_file,nc,C);  // ground state coefficients
            fclose(out_file);
        }
        free(C);
        free(Ho);
    }

    fclose(e_file);
}



void MIXTURE_SCANNING(int n_cases, char prefix [], unsigned int full_output)
{
    int
        i,
        j,
        nc,
        NparA,
        lmaxA,
        NparB,
        lmaxB,
        lan_it,
        total_mom,
        coefMemory;

    double
        l,
        E0;

    double
        g[3];

    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100];

    Carray
        C,
        HoA,
        HoB;

    FILE
        * e_file,
        * out_file;

    CompoundSpace
        MixSpace;

    // set file name with input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_energy.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Energy per particle | Npar A | Npar B | L | ");
    fprintf(e_file,"ga | gb | gab");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of config. space from line 'i' of input file
        mixParLine(in_fname,i,&NparA,&lmaxA,&NparB,&lmaxB,&total_mom,g);
        MixSpace = AllocCompBasis(NparA,NparB,lmaxA,lmaxB,total_mom);
        nc = MixSpace->size;
        printf("\n\n\nCOMPUTING GROUND STATE(MIXTURE) %d/%d\n",i+1,n_cases);
        printf("NA = %d lmaxA = %d | ",NparA,lmaxA);
        printf("NB = %d lmaxB = %d | ",NparB,lmaxB);
        printf("L = %d | space size = %d\n",total_mom,nc);
        printf("Interaction ga = %.10lf | gb = %.10lf | ",g[0],g[1]);
        printf("gab = %.10lf",g[2]);
        sepline();
        // set up (default)  number of Lanczos iterations
        // respecting the memory allowed (MEMORY_TOL) and
        // maximum number of iterations (MAX_LNCZS_IT)
        coefMemory = nc * sizeof(double complex);
        if (nc > MAX_LNCZS_IT)
        {
            lan_it = 10;
            while (lan_it*coefMemory < MEMORY_TOL && lan_it < MAX_LNCZS_IT)
            {
                lan_it++;
            }
        }
        else lan_it = nc; // Small config. space
        // set up one-body matrix elements
        HoA = carrDef(2*lmaxA+1);
        for (j = 0; j < 2*lmaxA+1; j++)
        {
            l = 2 * PI * (j-lmaxA);
            HoA[j] = 0.5*l*l;
        }
        HoB = carrDef(2*lmaxB+1);
        for (j = 0; j < 2*lmaxB+1; j++)
        {
            l = 2 * PI * (j-lmaxB);
            HoB[j] = 0.5*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        initGuess(nc,C);
        // call routine for the ground state
        E0 = BOSEBOSE_GS(lan_it,MixSpace,C,HoA,HoB,g);
        printf("\nAverage energy per particle = %.10lf",E0/(NparA+NparB));
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"\n%.10E ",E0/(NparA+NparB));
        fprintf(e_file,"%d %d %d ",NparA,NparB,total_mom);
        fprintf(e_file,"%.10E %.10E %.10E",g[0],g[1],g[2]);
        if (full_output)
        {
            // set output filename for coefficients
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            out_file = openFileWrite(out_fname);
            fprintf(out_file,"(%.10E+0.0j)",E0/(NparA+NparB));
            carrAppend(out_file,nc,C);
            fclose(out_file);
        }
        free(C);
        free(HoA);
        free(HoB);
        freeCompSpace(MixSpace);
    }

    fclose(e_file);
}



void BOSEFERMI_SCANNING(int n_cases, char prefix [], unsigned int full_output)
{
    int
        i,
        j,
        nc,
        NparB,
        lmaxB,
        NparF,
        lmaxF,
        lan_it,
        total_mom,
        coefMemory;

    double
        l,
        E0;

    double
        g[3];

    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100];

    Carray
        C,
        HoB,
        HoF;

    FILE
        * e_file,
        * out_file;

    BFCompoundSpace
        MixSpace;

    // set up file name with input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_energy.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Energy per particle | N bosons | N fermions | L | ");
    fprintf(e_file,"gb | gbf");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of config. space from line 'i' of input file
        mixParLine(in_fname,i,&NparB,&lmaxB,&NparF,&lmaxF,&total_mom,g);
        MixSpace = BoseFermiBasis(NparB,NparF,lmaxB,lmaxF,total_mom);
        nc = MixSpace->size;
        printf("\n\n\nCOMPUTING GROUND STATE(BOSE-FERMI) %d/%d\n",i+1,n_cases);
        printf("%d Bosons in lmax = %d | ",NparB,lmaxB);
        printf("%d Fermions in lmax = %d | ",NparF,lmaxF);
        printf("L = %d | space size = %d\n",total_mom,nc);
        printf("Interaction gb = %.10lf | gbf = %.10lf",g[0],g[2]);
        sepline();
        // set up (default)  number of Lanczos iterations
        // respecting the memory allowed (MEMORY_TOL) and
        // maximum number of iterations (MAX_LNCZS_IT)
        coefMemory = nc * sizeof(double complex);
        if (nc > MAX_LNCZS_IT)
        {
            lan_it = 2;
            while (lan_it*coefMemory < MEMORY_TOL && lan_it < MAX_LNCZS_IT)
            {
                lan_it++;
            }
        }
        else lan_it = nc; // Small config. space
        // set up one-body matrix elements
        HoB = carrDef(2*lmaxB+1);
        for (j = 0; j < 2*lmaxB+1; j++)
        {
            l = 2 * PI * (j-lmaxB);
            HoB[j] = 0.5*l*l;
        }
        HoF = carrDef(2*lmaxF+1);
        for (j = 0; j < 2*lmaxF+1; j++)
        {
            l = 2 * PI * (j-lmaxF);
            HoF[j] = 0.5*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        initGuess(nc,C);
        // call routine for the ground state
        E0 = BOSEFERMI_GS(lan_it,MixSpace,C,HoB,HoF,g);
        printf("\nAverage energy per particle = %.10lf",E0/(NparB+NparF));
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"\n%.10E ",E0/(NparB+NparF));
        fprintf(e_file,"%d %d %d ",NparB,NparF,total_mom);
        fprintf(e_file,"%.10E %.10E",g[0],g[2]);
        if (full_output)
        {
            // set file name for coefficients
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            out_file = openFileWrite(out_fname);
            fprintf(out_file,"(%.10E+0.0j)",E0/(NparB+NparF));
            carrAppend(out_file,nc,C);
            fclose(out_file);
        }
        free(C);
        free(HoB);
        free(HoF);
        freeBoseFermiSpace(MixSpace);
    }

    fclose(e_file);
}



#endif
