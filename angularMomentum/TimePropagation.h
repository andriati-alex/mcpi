#ifndef _TimePropagation_h
#define _TimePropagation_h

#include "LanczosRoutines.h"

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



#endif
