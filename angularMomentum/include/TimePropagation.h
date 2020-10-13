#ifndef _TimePropagation_h
#define _TimePropagation_h

#include "LanczosRoutines.h"
#include "MomentumObservables.h"

double TIME_EVOLUTION(int Nsteps, double Tstep, char Tinfo, double Trec,
       char prefix [], int jobId, int Npar, int lmax, int total_mom,
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
        m,
        nc,
        step,
        Nrec,
        Niter;
    double
        GSenergy,
        previousE;
    double complex
        dt;
    Iarray
        * ht;
    Rarray
        d,
        e,
        eigvec,
        avgocc;
    Carray
        aux,
        diag,
        offdiag,
        Clanczos;
    Cmatrix
        lvec;
    HConfMat
        H;
    char
        fname[100],
        jobIdstr[5];
    FILE
        * dataf,
        * timef;

    Niter = LNCZS_TIME_IT;

    // Configurational basis setup
    nc = BFixedMom_mcsize(Npar,lmax,total_mom);
    ht = BAssembleHT(Npar,lmax,total_mom,nc);
    H = boseH(Npar,lmax,nc,ht,Ho,g);

    avgocc = rarrDef(2*lmax+1);

    sprintf(jobIdstr,"%d",jobId);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_timecoef_job");
    strcat(fname,jobIdstr);
    strcat(fname,".dat");
    // open output file for coefficients of Fock space expansion
    dataf = openFileWrite(fname);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_timeobs_job");
    strcat(fname,jobIdstr);
    strcat(fname,".dat");
    // open output file for average occuaptions and time instants
    timef = openFileWrite(fname);

    // record initial data
    carr_inline(dataf,nc,C);

    // record time instants and average occupations
    fprintf(timef,"# time  |  Vector of average occupations\n");
    single_avgocc(Npar,lmax,nc,ht,C,avgocc);
    fprintf(timef,"%8.5lf",0.0);
    for (m = 0; m < 2*lmax+1; m++)
    {
        fprintf(timef," %13.10lf",avgocc[m]);
    }
    fprintf(timef,"\n");

    Nrec = 1;

    if (Tinfo == 'i' || Tinfo == 'I') dt = Tstep;
    else                              dt = I*Tstep;

    if (H == NULL) printf("\nWithout set Many-Body Hamiltonian matrix\n");
    else           printf("\nUsing (sparse) Hamiltonian matrix\n");

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
    aux = carrDef(Niter);
    Clanczos = carrDef(Niter);

    // Compute initial energy (GUESS)
    if (H == NULL) actH(lmax,nc,ht,Ho,g,C,lvec[0]);
    else           matmul(nc,H->rows,H->cols,H->vals,C,lvec[0]);
    previousE = creal(carrDot(nc,C,lvec[0]));

    // Initiate time evolution
    for (step = 0; step < Nsteps; step++)
    {
        /***   CALL LANCZOS FOR TRIDIAGONAL DECOMPOSITION   ***/
        for (i = 0; i < nc; i++) lvec[0][i] = C[i];
        if (H == NULL) LNCZS_TIME(nc,lmax,ht,Ho,g,diag,offdiag,lvec);
        else           LNCZS_HMAT_TIME(nc,H,diag,offdiag,lvec);

        // Transfer data to use lapack routine
        for (k = 0; k < Niter; k++)
        {
            d[k] = creal(diag[k]);    // Supposed to be real
            e[k] = creal(offdiag[k]); // Supposed to be real
            for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
        }

        /***   CALL LAPACK FOR TRIDIAGONAL LANCZOS OUTPUT   ***/
        k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
        if (k != 0)  LAPACK_PROBLEM(k,"stopLanczos");

        /*** USE TIME EVOLUTION IN LANCZOS SPACE AND THEN TRANSFORM BACK ***/
        Clanczos[0] = 1.0;
        for (k = 1; k < Niter; k++) Clanczos[k] = 0;

        for (k = 0; k < Niter; k++)
        {   // Solve in diagonal basis and for this apply eigvec trasformation
            aux[k] = eigvec[0*Niter + k] * Clanczos[0];
            // aux[k] = aux[k] * cexp(- I * d[k] * dt); // REAL TIME CASE
            aux[k] = aux[k] * cexp(- d[k] * dt);
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
        if (Tinfo == 'i' || Tinfo == 'I') carrNormalize(nc,C);

        // update energy
        if (H == NULL) actH(lmax,nc,ht,Ho,g,C,lvec[0]);
        else           matmul(nc,H->rows,H->cols,H->vals,C,lvec[0]);
        GSenergy = creal(carrDot(nc,C,lvec[0]));

        if ((step+1) % (Nsteps / 1000) == 0)
        {
            printf("\n %6d/%d   ",step+1,Nsteps);
            printf(" %10.6lf",(step+1)*Tstep);
            printf(" %16.9lf",GSenergy/Npar);
        }

        if (Tinfo == 'i' || Tinfo == 'I')
        {
            if ((step+1) % (Nsteps / 20) == 0)
            {
                if (1.0 - GSenergy/previousE < 1E-9)
                {
                    // ACHIEVED CONVERGENCE CRITERIA
                    fclose(dataf);
                    fclose(timef);
                    free(avgocc);
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
                    if (H != NULL) freeHmat(H);

                    printf("\n\nACHIEVED CONVERGENCE CRITERION");

                    return GSenergy;
                }
                else previousE = GSenergy;
            }
        }
        // record data if complete another record period
        if (fabs((step+1)*Tstep/Trec - Nrec) <= 1E-10)
        {
            single_avgocc(Npar,lmax,nc,ht,C,avgocc);
            fprintf(timef,"%8.5lf",(step+1)*Tstep);
            for (m = 0; m < 2*lmax+1; m++)
            {
                fprintf(timef," %13.10lf",avgocc[m]);
            }
            fprintf(timef,"\n");
            carr_inline(dataf,nc,C);
            Nrec++;
        }
    }



    fclose(dataf);
    fclose(timef);
    free(avgocc);
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
    if (H != NULL) freeHmat(H);

    return GSenergy;
}



double BB_TIME_EVOLUTION(int Nsteps, double Tstep, char Tinfo, double Trec,
       char prefix [], int jobId, CompoundSpace S, Carray C, Carray HoA,
       Carray HoB, double g [])
{

    int
        i,
        k,
        j,
        m,
        nc,
        Npar,
        step,
        Nrec,
        Niter;
    double
        avgLA,
        avgLB,
        variA,
        variB,
        covAB,
        GSenergy,
        previousE;
    double complex
        dt;
    Rarray
        d,
        e,
        eigvec,
        avgoccA,
        avgoccB;
    Carray
        aux,
        diag,
        offdiag,
        Clanczos;
    Cmatrix
        lvec;
    HConfMat
        H;
    char
        fname[100],
        jobIdstr[5];
    FILE
        * timef,
        * dataf;

    Niter = LNCZS_TIME_IT;
    nc = S->size;
    Npar = S->Na + S->Nb;

    avgoccA = rarrDef(2*S->lmaxA+1);
    avgoccB = rarrDef(2*S->lmaxB+1);

    sprintf(jobIdstr,"%d",jobId);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_timecoef_job");
    strcat(fname,jobIdstr);
    strcat(fname,".dat");
    // open output file for coefficients of Fock space expansion
    dataf = openFileWrite(fname);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_timeobs_job");
    strcat(fname,jobIdstr);
    strcat(fname,".dat");
    // open output file for some observables and time instants
    timef = openFileWrite(fname);

    // record initial data
    carr_inline(dataf,nc,C);

    // record some observables related to angular momentum
    mixture_avgocc(S,C,'A',avgoccA);
    mixture_avgocc(S,C,'B',avgoccB);
    avgLA = mixture_avgmom(S,C,'A');
    avgLB = mixture_avgmom(S,C,'B');
    variA = mixture_momvariance(S,C,'A');
    variB = mixture_momvariance(S,C,'B');
    covAB = mixture_momcov(S,C);
    fprintf(timef,"# time | avg La | avg Lb | var La | var Lb | cov LaLb | ");
    fprintf(timef,"vector of occ A | vector of occ B\n");
    fprintf(timef,"%8.5lf %13.10lf %13.10lf ",0.0,avgLA,avgLB);
    fprintf(timef,"%13.10lf %13.10lf %13.10lf",variA,variB,covAB);
    for (m = 0; m < 2*S->lmaxA+1; m++)
    {
        fprintf(timef," %13.10lf",avgoccA[m]);
    }
    for (m = 0; m < 2*S->lmaxB+1; m++)
    {
        fprintf(timef," %13.10lf",avgoccB[m]);
    }
    fprintf(timef,"\n");

    Nrec = 1;

    // select type of time propagation
    if (Tinfo == 'i' || Tinfo == 'I') dt = Tstep;
    else                              dt = I*Tstep;

    // TRY TO ALLOCATE MATRIX. IF THE MEMORY REQUIRED EXCEED
    // THE TOLERANCE DEFINED RETURN NULL POINTER
    H = boseboseH(S,HoA,HoB,g);
    if (H == NULL) printf("\nWithout set Hamiltonian matrix\n");
    else           printf("\nUsing (sparse) Hamiltonian matrix\n");

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
    aux = carrDef(Niter);
    Clanczos = carrDef(Niter);

    // Compute initial energy (GUESS)
    if (H == NULL) mixture_actH(S,HoA,HoB,g,C,lvec[0]);
    else           matmul(nc,H->rows,H->cols,H->vals,C,lvec[0]);
    previousE = creal(carrDot(nc,C,lvec[0]));

    // Initiate time evolution
    for (step = 0; step < Nsteps; step++)
    {
        /***   CALL LANCZOS FOR TRIDIAGONAL DECOMPOSITION   ***/
        for (i = 0; i < nc; i++)  lvec[0][i] = C[i];
        if (H == NULL) LNCZS_BBMIX_TIME(S,HoA,HoB,g,diag,offdiag,lvec);
        else           LNCZS_HMAT_TIME(nc,H,diag,offdiag,lvec);

        // Transfer data to use lapack routine
        for (k = 0; k < Niter; k++)
        {
            d[k] = creal(diag[k]);    // Supposed to be real
            e[k] = creal(offdiag[k]); // Supposed to be real
            for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
        }

        /***   CALL LAPACK FOR TRIDIAGONAL LANCZOS OUTPUT   ***/
        k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
        if (k != 0)  LAPACK_PROBLEM(k,"stopLanczos");

        /*** USE TIME EVOLUTION IN LANCZOS SPACE AND THEN TRANSFORM BACK ***/
        Clanczos[0] = 1.0;
        for (k = 1; k < Niter; k++) Clanczos[k] = 0;

        for (k = 0; k < Niter; k++)
        {   // Solve in diagonal basis and for this apply eigvec trasformation
            aux[k] = eigvec[0*Niter + k] * Clanczos[0];
            // aux[k] = aux[k] * cexp(- I * d[k] * dt); // REAL TIME CASE
            aux[k] = aux[k] * cexp(- d[k] * dt);
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
        if (Tinfo == 'i' || Tinfo == 'I') carrNormalize(nc,C);

        // update energy
        if (H == NULL) mixture_actH(S,HoA,HoB,g,C,lvec[0]);
        else           matmul(nc,H->rows,H->cols,H->vals,C,lvec[0]);
        GSenergy = creal(carrDot(nc,C,lvec[0]));

        if ((step+1) % (Nsteps / 1000) == 0)
        {
            printf("\n %6d/%d   ",step+1,Nsteps);
            printf(" %10.6lf",(step+1)*Tstep);
            printf(" %16.9lf",GSenergy/Npar);
        }

        if (Tinfo == 'i' || Tinfo == 'I')
        {
            if ((step+1) % (Nsteps / 20) == 0)
            {
                if (1.0 - GSenergy/previousE < 1E-9)
                {
                    // ACHIEVED CONVERGENCE CRITERIA
                    fclose(dataf);
                    fclose(timef);
                    free(avgoccA);
                    free(avgoccB);
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
                    if (H != NULL) freeHmat(H);
                    printf("\n\nACHIEVED CONVERGENCE CRITERION\n");
                    return GSenergy;
                }
                else previousE = GSenergy;
            }
        }

        // record data if complete another record period
        if (fabs((step+1)*Tstep/Trec - Nrec) <= 1E-10)
        {
            mixture_avgocc(S,C,'A',avgoccA);
            mixture_avgocc(S,C,'B',avgoccB);
            avgLA = mixture_avgmom(S,C,'A');
            avgLB = mixture_avgmom(S,C,'B');
            variA = mixture_momvariance(S,C,'A');
            variB = mixture_momvariance(S,C,'B');
            covAB = mixture_momcov(S,C);

            fprintf(timef,"%8.5lf ",(step+1)*Tstep);
            fprintf(timef,"%13.10lf %13.10lf ",avgLA,avgLB);
            fprintf(timef,"%13.10lf %13.10lf %13.10lf",variA,variB,covAB);
            for (m = 0; m < 2*S->lmaxA+1; m++)
            {
                fprintf(timef," %13.10lf",avgoccA[m]);
            }
            for (m = 0; m < 2*S->lmaxB+1; m++)
            {
                fprintf(timef," %13.10lf",avgoccB[m]);
            }
            fprintf(timef,"\n");
            carr_inline(dataf,nc,C);
            Nrec++;
        }
    }

    fclose(dataf);
    fclose(timef);
    free(avgoccA);
    free(avgoccB);
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
    if (H != NULL) freeHmat(H);

    return GSenergy;
}



double BF_TIME_EVOLUTION(int Nsteps, double Tstep, char Tinfo, double Trec,
       FILE * dataf, BFCompoundSpace S, Carray C, Carray HoB,
       Carray HoF, double g [])
{

    int
        i,
        k,
        j,
        nc,
        Npar,
        step,
        Nrec,
        Niter;

    double
        GSenergy,
        previousE;

    double complex
        dt;

    Rarray
        d,
        e,
        eigvec;

    Carray
        aux,
        diag,
        offdiag,
        Clanczos;

    Cmatrix
        lvec;

    HConfMat
        H;

    Niter = LNCZS_TIME_IT;
    nc = S->size;
    Npar = S->Nb + S->Nf;

    if (Tinfo == 'i' || Tinfo == 'I') dt = Tstep;
    else                              dt = I*Tstep;

    if (Tinfo == 'r' || Tinfo == 'R')
    {
        fprintf(dataf,"(0.0+0.0) ");
        carr_inline(dataf,nc,C);
    }
    Nrec = 1;

    // TRY TO ALLOCATE MATRIX. IF THE MEMORY REQUIRED EXCEED
    // THE TOLERANCE DEFINED RETURN NULL POINTER
    H = bosefermiH(S,HoB,HoF,g);
    if (H == NULL) printf("Without set Hamiltonian matrix\n");
    else           printf("Using (sparse) Hamiltonian matrix\n");

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
    aux = carrDef(Niter);
    Clanczos = carrDef(Niter);

    // Compute initial energy (GUESS)
    if (H == NULL) bosefermi_actH(S,HoB,HoF,g,C,lvec[0]);
    else           matmul(nc,H->rows,H->cols,H->vals,C,lvec[0]);
    previousE = creal(carrDot(nc,C,lvec[0]));

    // Initiate time evolution
    for (step = 0; step < Nsteps; step++)
    {
        /***   CALL LANCZOS FOR TRIDIAGONAL DECOMPOSITION   ***/
        for (i = 0; i < nc; i++)  lvec[0][i] = C[i];
        if (H == NULL) LNCZS_BFMIX_TIME(S,HoB,HoF,g,diag,offdiag,lvec);
        else           LNCZS_HMAT_TIME(nc,H,diag,offdiag,lvec);

        // Transfer data to use lapack routine
        for (k = 0; k < Niter; k++)
        {
            d[k] = creal(diag[k]);    // Supposed to be real
            e[k] = creal(offdiag[k]); // Supposed to be real
            for (j = 0; j < Niter; j++) eigvec[k * Niter + j] = 0;
        }

        /***   CALL LAPACK FOR TRIDIAGONAL LANCZOS OUTPUT   ***/
        k = LAPACKE_dstev(LAPACK_ROW_MAJOR,'V',Niter,d,e,eigvec,Niter);
        if (k != 0)  LAPACK_PROBLEM(k,"stopLanczos");

        /*** USE TIME EVOLUTION IN LANCZOS SPACE AND THEN TRANSFORM BACK ***/
        Clanczos[0] = 1.0;
        for (k = 1; k < Niter; k++) Clanczos[k] = 0;

        for (k = 0; k < Niter; k++)
        {   // Solve in diagonal basis and for this apply eigvec trasformation
            aux[k] = eigvec[0*Niter + k] * Clanczos[0];
            // aux[k] = aux[k] * cexp(- I * d[k] * dt); // REAL TIME CASE
            aux[k] = aux[k] * cexp(- d[k] * dt);
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
        if (Tinfo == 'i' || Tinfo == 'I') carrNormalize(nc,C);

        // update energy
        if (H == NULL) bosefermi_actH(S,HoB,HoF,g,C,lvec[0]);
        else           matmul(nc,H->rows,H->cols,H->vals,C,lvec[0]);
        GSenergy = creal(carrDot(nc,C,lvec[0]));

        if ((step+1) % (Nsteps / 1000) == 0)
        {
            printf("\n %6d/%d   ",step+1,Nsteps);
            printf(" %10.6lf",(step+1)*Tstep);
            printf(" %16.9lf",GSenergy/Npar);
        }

        if (Tinfo == 'i' || Tinfo == 'I')
        {
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
                    if (H != NULL) freeHmat(H);
                    printf("\n\nACHIEVED CONVERGENCE CRITERION\n");
                    return GSenergy;
                }
                else previousE = GSenergy;
            }
        }
        else
        {
            // record data if complete another record period
            if (fabs((step+1)*Tstep/Trec - Nrec) <= 1E-10)
            {
                fprintf(dataf,"(%.5lf+0.0) ",(step+1)*Tstep);
                carr_inline(dataf,nc,C);
                Nrec++;
            }
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
    if (H != NULL) freeHmat(H);

    return GSenergy;
}



void TIME_SCANNING(int n_cases, char prefix [], double Trec, char Tinfo)
{
    int
        i,
        j,
        nc,
        Npar,
        lmax,
        Nsteps,
        total_mom;
    double
        l,
        g,
        E0,
        Tstep,
        boost;
    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100],
        init_fname[100];
    Carray
        C,
        Ho;
    FILE
        * e_file;

    // set up name of file containing input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set up name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_timesetup.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Energy per particle | Num. of Particles | lmax ");
    fprintf(e_file,"| total mom. | boost | g");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of configurational
        // space in line 'i' of input file
        parLine_time(in_fname,i,&Npar,&lmax,&total_mom,&boost,&g,
                     &Nsteps,&Tstep);
        // number of configuration(nc) - config. space dimension
        nc = BFixedMom_mcsize(Npar,lmax,total_mom);
        // set up one-body matrix elements
        Ho = carrDef(2*lmax+1);
        for (j = 0; j < 2*lmax+1; j++)
        {
            l = (j-lmax) - boost;
            Ho[j] = 0.5*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        if (Tinfo == 'r' || Tinfo == 'R')
        {
            // set file name with input data
            strcpy(init_fname,"inp_data/");
            strcat(init_fname,prefix);
            strcat(init_fname,"_job");
            strcat(init_fname,strnum);
            strcat(init_fname,".dat");
            // read the input data
            carr_input_txt(init_fname,nc,C);
        }
        else
        {
            initGuess(nc,C);
        }
        printf("\n\n\n\n\nCOMPUTING TIME EVOLUTION %d/%d",i+1,n_cases);
        sepline();
        // call (main)routine to compute the ground state
        E0 = TIME_EVOLUTION(Nsteps,Tstep,Tinfo,Trec,prefix,i+1,
                            Npar,lmax,total_mom,C,Ho,g);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"\n%.10E %d %d %d ",E0/Npar,Npar,lmax,total_mom);
        fprintf(e_file,"%.10E %.10E",boost,g);
        if (Tinfo == 'i' || Tinfo == 'I')
        {
            // set output file name
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_itimefinal_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // record output ground-state from imaginary time
            carr_txt(out_fname,nc,C);
        }
        free(C);
        free(Ho);
    }

    fclose(e_file);
}



void adjustInputProd(int n_cases, char prefix [])
{

    int
        i,
        j,
        k,
        n,
        m,
        ia,
        La,
        Lb,
        NparA,
        NparB,
        lmaxA,
        lmaxB,
        sizeA,
        sizeB,
        total_size;
    double
        g,      // useless
        boost;  // useless - just to use parLine function
    char
        strnum[5],
        fname[100],
        fnameA[100],
        fnameB[100];
    Iarray
        * htA,
        * htB;
    Carray
        Ca,
        Cb,
        C;
    CompoundSpace
        S;

    // set file name with input data
    strcpy(fnameA,prefix);
    strcat(fnameA,"_specA.inp");
    strcpy(fnameB,prefix);
    strcat(fnameB,"_specB.inp");

    for (i = 0; i < n_cases; i++)
    {
        sprintf(strnum,"%d",i+1);
        parLine(fnameA,i,&NparA,&lmaxA,&La,&boost,&g);
        parLine(fnameB,i,&NparB,&lmaxB,&Lb,&boost,&g);
        // allocate multi-config. space information
        sizeA = BFixedMom_mcsize(NparA,lmaxA,La);
        sizeB = BFixedMom_mcsize(NparA,lmaxA,Lb);
        htA = BAssembleHT(NparA,lmaxA,La,sizeA);
        htB = BAssembleHT(NparB,lmaxB,Lb,sizeB);
        S = AllocCompBasis(NparA,NparB,lmaxA,lmaxB,La+Lb);
        total_size = S->size;
        // allocate vector of coefficients
        Ca = carrDef(sizeA);
        Cb = carrDef(sizeB);
        C = carrDef(total_size);
        // set file name with input data of species A
        strcpy(fname,"inp_data/");
        strcat(fname,prefix);
        strcat(fname,"_job");
        strcat(fname,strnum);
        strcat(fname,"_specA.dat");
        // read the input data
        carr_input_txt(fname,sizeA,Ca);
        // set file name with input data of species B
        strcpy(fname,"inp_data/");
        strcat(fname,prefix);
        strcat(fname,"_job");
        strcat(fname,strnum);
        strcat(fname,"_specB.dat");
        // read the input data
        carr_input_txt(fname,sizeB,Cb);
        // transfer data to the compound space coefficients
        for (k = 0; k < total_size; k++) C[k] = 0.0 + 0.0*I;
        for (k = 0; k < sizeA; k++)
        {
            n = BBsubIndex(S,htA[k]);
            ia = BgetIndex(S->lmaxA,S->sub[n].sizeA,S->sub[n].hta,htA[k]);
            for (j = 0; j < sizeB; j++)
            {
                m = BBgetIndex_B(S,n,ia,htB[j]);
                C[m] = Ca[k]*Cb[j];
            }
        }
        assert_norm(total_size,C);
        // record data in text file of coefficients in compound space
        strcpy(fname,"inp_data/");
        strcat(fname,prefix);
        strcat(fname,"_job");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,total_size,C);
        // free memory
        for (k = 0; k < sizeA; k++) free(htA[k]);
        free(htA);
        for (k = 0; k < sizeB; k++) free(htB[k]);
        free(htB);
        free(Ca);
        free(Cb);
        free(C);
        freeCompSpace(S);
    }
}



void TIME_MIXTURE_SCANNING(int n_cases, char prefix [], double Trec, char Tinfo,
     unsigned int inputType)
{
    int
        i,
        j,
        nc,
        NparA,
        lmaxA,
        NparB,
        lmaxB,
        total_mom,
        Nsteps;

    double
        l,
        E0,
        Tstep,
        boost,
        MassImbal;

    double
        g[3];

    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100],
        init_fname[100];

    Carray
        C,
        HoA,
        HoB;

    FILE
        * e_file;

    CompoundSpace
        MixSpace;

    if (inputType == 1) adjustInputProd(n_cases,prefix);

    // set file name with input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_timesetup.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Energy per particle | N bosons A | lmax A ");
    fprintf(e_file,"| N bosons B | lmax B | Total Mom. | boost B ");
    fprintf(e_file,"| Ma / Mb | ga | gb | gab");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of config. space from line 'i' of input file
        mixParLine_time(in_fname,i,&NparA,&lmaxA,&NparB,&lmaxB,&total_mom,
                &boost,&MassImbal,g,&Nsteps,&Tstep);
        MixSpace = AllocCompBasis(NparA,NparB,lmaxA,lmaxB,total_mom);
        nc = MixSpace->size;
        // set up one-body matrix elements
        HoA = carrDef(2*lmaxA+1);
        for (j = 0; j < 2*lmaxA+1; j++)
        {
            l = j-lmaxA;
            HoA[j] = 0.5*l*l;
        }
        HoB = carrDef(2*lmaxB+1);
        for (j = 0; j < 2*lmaxB+1; j++)
        {
            l = (j-lmaxB) - boost/MassImbal;
            HoB[j] = 0.5*MassImbal*l*l;
        }
        // Allocate many-body state coefficients
        C = carrDef(nc);
        if (Tinfo == 'r' || Tinfo == 'R')
        {
            // set file name with input data
            strcpy(init_fname,"inp_data/");
            strcat(init_fname,prefix);
            strcat(init_fname,"_job");
            strcat(init_fname,strnum);
            strcat(init_fname,".dat");
            // read the input data
            carr_input_txt(init_fname,nc,C);
        }
        else
        {
            initGuess(nc,C);
        }
        printf("\n\n\n\n\nTIME EVOLUTION(MIXTURE) %d/%d",i+1,n_cases);
        sepline();
        // call routine for the ground state
        E0 = BB_TIME_EVOLUTION(Nsteps,Tstep,Tinfo,Trec,prefix,i+1,
             MixSpace,C,HoA,HoB,g);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"\n%.10E ",E0/(NparA+NparB));
        fprintf(e_file,"%d %d %d %d %d ",NparA,lmaxA,NparB,lmaxB,total_mom);
        fprintf(e_file,"%.10E %.10E ",boost,MassImbal);
        fprintf(e_file,"%.10E %.10E %.10E",g[0],g[1],g[2]);
        // set output filename for coefficients
        if (Tinfo == 'i' || Tinfo == 'I')
        {
            // set output file name
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_itimefinal_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            carr_txt(out_fname,nc,C);
        }
        free(C);
        free(HoA);
        free(HoB);
        freeCompSpace(MixSpace);
    }

    fclose(e_file);
}



void TIME_BOSEFERMI_SCANNING(int n_cases, char prefix [],
     double Trec, char Tinfo)
{
    int
        i,
        j,
        nc,
        NparB,
        lmaxB,
        NparF,
        lmaxF,
        total_mom,
        Nsteps;

    double
        l,
        E0,
        Tstep,
        boost,
        MassImbal;

    double
        g[3];

    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100],
        init_fname[100];

    Carray
        C,
        HoF,
        HoB;

    FILE
        * e_file,
        * c_file;

    BFCompoundSpace
        MixSpace;

    // set file name with input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_timesetup.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Energy per particle | N bosons | lmax Bosons ");
    fprintf(e_file,"| N Fermions | lmax Fermions | Total Mom. ");
    fprintf(e_file,"| Boost Fermions | Mb / Mf | gb | gbf");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of config. space from line 'i' of input file
        mixParLine_time(in_fname,i,&NparB,&lmaxB,&NparF,&lmaxF,&total_mom,
                &boost,&MassImbal,g,&Nsteps,&Tstep);
        MixSpace = BoseFermiBasis(NparB,NparF,lmaxB,lmaxF,total_mom);
        nc = MixSpace->size;
        // set up one-body matrix elements
        HoB = carrDef(2*lmaxB+1);
        for (j = 0; j < 2*lmaxB+1; j++)
        {
            l = j-lmaxB;
            HoB[j] = 0.5*l*l;
        }
        HoF = carrDef(2*lmaxF+1);
        for (j = 0; j < 2*lmaxF+1; j++)
        {
            l = (j-lmaxF) - boost/MassImbal;
            HoF[j] = 0.5*MassImbal*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        if (Tinfo == 'r' || Tinfo == 'R')
        {
            // set file name with input data
            strcpy(init_fname,"inp_data/");
            strcat(init_fname,prefix);
            strcat(init_fname,"_job");
            strcat(init_fname,strnum);
            strcat(init_fname,".dat");
            // read the input data
            carr_input_txt(init_fname,nc,C);
            // set output file name
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_coef_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file for coefficients of Fock states
            c_file = openFileWrite(out_fname);
        }
        else
        {
            initGuess(nc,C);
            c_file = NULL;
        }
        printf("\n\n\n\n\nCOMPUTING TIME EVOLUTION OF ");
        printf("BOSE-FERMI MIXTURE %d/%d",i+1,n_cases);
        sepline();
        // call routine for the ground state
        E0 = BF_TIME_EVOLUTION(Nsteps,Tstep,Tinfo,Trec,c_file,
             MixSpace,C,HoB,HoF,g);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"\n%.10E ",E0/(NparB+NparF));
        fprintf(e_file,"%d %d %d %d %d ",NparB,lmaxB,NparF,lmaxF,total_mom);
        fprintf(e_file,"%.10E %.10E ",boost,MassImbal);
        fprintf(e_file,"%.10E %.10E",g[0],g[2]);
        // set output filename for coefficients
        if (Tinfo == 'i' || Tinfo == 'I')
        {
            // set output file name
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            carr_txt(out_fname,nc,C);
        }
        if (c_file != NULL) fclose(c_file);
        free(C);
        free(HoB);
        free(HoF);
        freeBoseFermiSpace(MixSpace);
    }

    fclose(e_file);
}



#endif
