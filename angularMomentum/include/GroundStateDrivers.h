#ifndef _GroundStateDrivers_h
#define _GroundStateDrivers_h

#include "LanczosRoutines.h"
#include "MomentumObservables.h"



/** FUNCTIONS CALLED IN THE EXECUTABLE PROGRAMS
    ===========================================
    These functions provide an interface to read parameters
    call Lanczos algorithm and record output data.      **/






double GROUND_STATE(int Niter, int Npar, int lmax, int total_mom, int nc,
       Iarray * ht, Carray C, Carray Ho, double g)
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

    // empty system
    if (Npar == 0) return 0;

    // Assert the config. space is not void
    if (nc == 0)
    {
        printf("\n\nDIAGONALIZATION PROCESS ABORTED : The parameters ");
        printf("requested, Npar = %d, lmax = %d, L = %d",Npar,lmax,total_mom);
        printf(" generated empty config. space.\n\n");
        exit(EXIT_FAILURE);
    }

    H = boseH(Npar,lmax,nc,ht,Ho,g);
    if (H == NULL) printf("\nWithout set Hamiltonian matrix\n");
    else           printf("\nUsing (sparse) Hamiltonian matrix\n");

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
    if (H == NULL) printf("\nWithout set Hamiltonian matrix\n");
    else           printf("\nUsing (sparse) Hamiltonian matrix\n");

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
    if (H == NULL) printf("\nWithout set Hamiltonian matrix\n");
    else           printf("\nUsing (sparse) Hamiltonian matrix\n");

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
        E0,
        boost;
    char
        strnum[10],
        e_fname[100],
        in_fname[100],
        out_fname[100];
    Iarray
        * ht;
    Rarray
        avgocc;
    Carray
        C,
        Ho;
    FILE *
        e_file;

    // set up name of file containing input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set up name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_setup.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# Num. of Particles | lmax ");
    fprintf(e_file,"| total mom. | boost | g | Energy per Particle\n");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of configurational
        // space in line 'i' of input file
        parLine(in_fname,i,&Npar,&lmax,&total_mom,&boost,&g);
        // number of configuration(nc) - config. space dimension
        nc = BFixedMom_mcsize(Npar,lmax,total_mom);
        ht = BAssembleHT(Npar,lmax,total_mom,nc);
        printf("\n\n\n\n\nCOMPUTING GROUND STATE %d/%d",i+1,n_cases);
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
            l = (j-lmax) - boost;
            Ho[j] = 0.5*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        initGuess(nc,C);
        // call (main)routine to compute the ground state
        E0 = GROUND_STATE(lan_it,Npar,lmax,total_mom,nc,ht,C,Ho,g);
        // compute average occupations in single particle states
        avgocc = rarrDef(2*lmax+1);
        single_avgocc(Npar,lmax,nc,ht,C,avgocc);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"%d %d %d ",Npar,lmax,total_mom);
        fprintf(e_file,"%.10E %.10E %.10E",boost,g,E0/Npar);
        for (j = 0; j < 2*lmax+1; j++) fprintf(e_file," %.10lf",avgocc[j]);
        fprintf(e_file,"\n");
        if (full_output)
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
        free(avgocc);
        free(C);
        free(Ho);
        for (j = 0; j < nc; j++) free(ht[j]);
        free(ht);
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
        E0,
        boost,
        MassImbal,
        avgLA,
        avgLB,
        variA,
        variB,
        covAB;
    double
        g[3];
    Rarray
        avgoccA,
        avgoccB;
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
        * e_file;
    CompoundSpace
        MixSpace;

    // set file name with input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_setup.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# N bosons A | lmax A ");
    fprintf(e_file,"| N bosons B | lmax B | Total Mom. | boost B ");
    fprintf(e_file,"| Ma / Mb | ga | gb | gab | Energy per particle\n");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of config. space from line 'i' of input file
        mixParLine(in_fname,i,&NparA,&lmaxA,&NparB,&lmaxB,&total_mom,
                   &boost,&MassImbal,g);
        MixSpace = AllocCompBasis(NparA,NparB,lmaxA,lmaxB,total_mom);
        nc = MixSpace->size;
        printf("\n\n\n\n\nCOMPUTING GROUND STATE ");
        printf("OF BOSONIC MIXTURE %d/%d",i+1,n_cases);
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
            l = j-lmaxA;
            HoA[j] = 0.5*l*l;
        }
        HoB = carrDef(2*lmaxB+1);
        for (j = 0; j < 2*lmaxB+1; j++)
        {
            l = (j-lmaxB) - boost/MassImbal;
            HoB[j] = 0.5*MassImbal*l*l;
        }
        // initialize many-body state coefficients
        C = carrDef(nc);
        initGuess(nc,C);
        // call routine for the ground state
        E0 = BOSEBOSE_GS(lan_it,MixSpace,C,HoA,HoB,g);
        // compute some observables
        avgoccA = rarrDef(2*lmaxA+1);
        avgoccB = rarrDef(2*lmaxB+1);
        mixture_avgocc(MixSpace,C,'A',avgoccA);
        mixture_avgocc(MixSpace,C,'B',avgoccB);
        avgLA = mixture_avgmom(MixSpace,C,'A');
        avgLB = mixture_avgmom(MixSpace,C,'B');
        variA = mixture_momvariance(MixSpace,C,'A');
        variB = mixture_momvariance(MixSpace,C,'B');
        covAB = mixture_momcov(MixSpace,C);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"%d %d %d %d %d ",NparA,lmaxA,NparB,lmaxB,total_mom);
        fprintf(e_file,"%.10E %.10E ",boost,MassImbal);
        fprintf(e_file,"%.10E %.10E %.10E ",g[0],g[1],g[2]);
        fprintf(e_file,"%.10E ",E0/(NparA+NparB));
        fprintf(e_file,"%13.10lf %13.10lf ",avgLA,avgLB);
        fprintf(e_file,"%13.10lf %13.10lf %13.10lf",variA,variB,covAB);
        for (j = 0; j < 2*lmaxA+1; j++)
        {
            fprintf(e_file," %.10lf",avgoccA[j]);
        }
        for (j = 0; j < 2*lmaxB+1; j++)
        {
            fprintf(e_file," %.10lf",avgoccB[j]);
        }
        fprintf(e_file,"\n");
        if (full_output)
        {
            // set output filename for coefficients
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            carr_txt(out_fname,nc,C);
        }
        free(avgoccA);
        free(avgoccB);
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
        E0,
        boost,
        MassImbal;

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

    FILE *
        e_file;

    BFCompoundSpace
        MixSpace;

    // set up file name with input data
    strcpy(in_fname,prefix);
    strcat(in_fname,".inp");

    // set name of file containing energy
    strcpy(e_fname,"output/");
    strcat(e_fname,prefix);
    strcat(e_fname,"_setup.dat");
    e_file = openFileWrite(e_fname);
    fprintf(e_file,"# N bosons | lmax Bosons ");
    fprintf(e_file,"| N Fermions | lmax Fermions | Total Mom. ");
    fprintf(e_file,"| Boost Fermions | Mb / Mf | gb | gbf | ");
    fprintf(e_file,"Energy per particle\n");

    for (i = 0; i < n_cases; i++)
    {
        // number of line as string to append in output file name
        sprintf(strnum,"%d",i+1);
        // set up parameters of config. space from line 'i' of input file
        mixParLine(in_fname,i,&NparB,&lmaxB,&NparF,&lmaxF,&total_mom,
                   &boost,&MassImbal,g);
        MixSpace = BoseFermiBasis(NparB,NparF,lmaxB,lmaxF,total_mom);
        nc = MixSpace->size;
        printf("\n\n\n\n\nCOMPUTING GROUND STATE ");
        printf("OF BOSE-FERMI MIXTURE %d/%d",i+1,n_cases);
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
        initGuess(nc,C);
        // call routine for the ground state
        E0 = BOSEFERMI_GS(lan_it,MixSpace,C,HoB,HoF,g);
        // WRITE OUTPUT DATA IN A FILE
        fprintf(e_file,"%d %d %d %d %d ",NparB,lmaxB,NparF,lmaxF,total_mom);
        fprintf(e_file,"%.10E %.10E ",boost,MassImbal);
        fprintf(e_file,"%.10E %.10E ",g[0],g[2]);
        fprintf(e_file,"%.10E\n",E0/(NparB+NparF));
        if (full_output)
        {
            // set output filename for coefficients
            strcpy(out_fname,"output/");
            strcat(out_fname,prefix);
            strcat(out_fname,"_job");
            strcat(out_fname,strnum);
            strcat(out_fname,".dat");
            // open output file
            carr_txt(out_fname,nc,C);
        }
        free(C);
        free(HoB);
        free(HoF);
        freeBoseFermiSpace(MixSpace);
    }

    fclose(e_file);
}



#endif
