
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : February/03/2020

 -------------------------------------------------------------------------

 ****  ROUTINE TO COMPUTE GROUND STATE ENERGY OF LIEB GAS

 * COMPILE :
 * ---------

 * icc -o exe groundStateLieb.c \
            -L${MKLROOT}/lib/intel64 \
            -lmkl_intel_lp64 \
            -lmkl_gnu_thread \
            -lmkl_core \
            -qopenmp \
            -lgomp \
            -O3 \

 * gcc -o exe groundStateLieb.c \
            -L${MKLROOT}/lib/intel64 \
            -lmkl_intel_lp64 \
            -lmkl_gnu_thread \
            -lmkl_core \
            -fopenmp \
            -lgomp \
            -lm \
            -O3 \



 * HOW TO EXECUTE :
 * ----------------

 * ./exe Nparticles Norbitals

 * where Nparticles and Morbitals are command line arguments for the
 * number of particles and individual particle states respectively.
 * Perform 1/8 of the size of the configurational space as lanczos
 * iterations.  To  compare  with the Fermi energy a odd number of
 * particles is required to avoid degeneracy



 * DEPENDENCIES :
 * --------------

 * It uses the routine  'dstev' from LAPACK package.  Here the compilation
 * options are directed to Intel Math Kernel Library(MKL) distribution for
 * Linux based systems. For other distributions of  LAPACK  or operational
 * system one must search the correct way to link the libraries.
 * More information on how to link libraries from Intel MKL are available:

 https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
 (last consulted on Feb/01/2020)

 * ----------------------------------------------------------------------- */

#define PI 3.141592653589793

#include "onebodyMatrix.h"
#include "twobodyMatrix.h"
#include "hamiltonianMatrix.h"

#include <mkl.h>



double complex carrDot(int n, Carray v1, Carray v2)
{

/** Convetional scalar product for complex vectors **/

    int i;

    double complex z = 0;

    for (i = 0; i < n; i++) z = z + conj(v1[i]) * v2[i];

    return z;
}



double carrMod(int n, Carray v)
{
    int i;

    double mod = 0;

    for (i = 0; i < n; i++)
    {
        mod = mod + creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }

    return sqrt(mod);
}



void setupOrbitals(int Morb, int Mpos, double L, Rarray x, Cmatrix orbs)
{

/** Setup orbitals as plane waves in a periodic domain of length L **/

    int
        n,
        j;

    double
        k;

    for (n = 0; n < Morb; n++)
    {
        k = 2 * PI * (n - Morb / 2) / L;
        for (j = 0; j < Mpos; j++)
        {
            orbs[n][j] = cexp(I * k * x[j]) / sqrt(L);
        }
    }
}



void setupHo(int Morb, double L, Cmatrix Ho)
{

/** Configure hamiltonian matrix using the individual particle states
    The plane waves are eigenstates of the one-body hamiltonian, thus
    the matrix is evaluated analytically **/

    int
        n,
        j;

    double
        k;

    for (n = 0; n < Morb; n++)
    {
        k = 2 * PI * (n - Morb / 2) / L;
        Ho[n][n] = 0.5 * k * k;

        for (j = n + 1; j < Morb; j++)
        {
            Ho[n][j] = 0.0;
            Ho[j][n] = 0.0;
        }
    }

}



void setupHint(int Morb, double g, Carray Hint)
{

/** Configure two-body hamiltonian matrix elements for contact interaction
    using the plane waves as orbitals. From momentum conservation evaluate
    the elements analytically **/

    int
        n,
        s,
        q,
        l,
        nk,
        sk,
        qk,
        lk,
        M,
        M2,
        M3;

    M  = Morb;
    M2 = M * M;
    M3 = M * M2;

    for (n = 0; n < Morb; n++)
    {
        nk = (n - Morb / 2);

        for (s = 0; s < Morb; s++)
        {
            sk = (s - Morb / 2);

            for (q = 0; q < Morb; q++)
            {
                qk = (q - Morb / 2);

                for (l = 0; l < Morb; l++)
                {
                    lk = (l - Morb / 2);

                    // exploit momentum conservation
                    if (lk + qk - sk - nk == 0)
                    {
                        Hint[n + s * M + q * M2 + l * M3] = g;
                    }
                    else
                    {
                        Hint[n + s * M + q * M2 + l * M3] = 0.0;
                    }
                }
            }
        }
    }

}



int lanczos(int lm, int Npar, int Morb, Iarray IF, Iarray strideOT,
            Iarray strideTT, Iarray Map, Iarray MapOT, Iarray MapTT,
            Cmatrix Ho, Carray Hint, Carray diag, Carray offdiag, Cmatrix lvec)
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
        number of itertion done (just lm if the method does not breakdown) **/



    int
        i,
        j,
        k,
        nc;

    double
        tol,
        maxCheck;

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



    nc = NC(Npar,Morb);
    HC = carrDef(nc),
    ortho = carrDef(lm);



    // Initiate the method
    applyHconf_omp(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,IF,
                   lvec[0],Ho,Hint,HC);
    diag[0] = carrDot(nc, lvec[0], HC);
    for (j = 0; j < nc; j++) HC[j] = HC[j] - diag[0] * lvec[0][j];

    // Core iteration procedure
    for (i = 0; i < lm - 1; i++)
    {
        offdiag[i] = carrMod(nc, HC);

        if (maxCheck < creal(offdiag[i])) maxCheck = creal(offdiag[i]);

        // If method break return number of iterations achieved
        if (creal(offdiag[i]) / maxCheck < tol) return i;

        for (j = 0; j < nc; j++) lvec[i+1][j] = HC[j] / offdiag[i];
        // carrScalarMultiply(nc, HC, 1.0 / offdiag[i], lvec[i + 1]);

        applyHconf_omp(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,IF,
                   lvec[i+1],Ho,Hint,HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - offdiag[i] * lvec[i][j];
        }

        diag[i + 1] = carrDot(nc, lvec[i + 1], HC);

        for (j = 0; j < nc; j++)
        {
            HC[j] = HC[j] - diag[i+1]*lvec[i+1][j];
        }

        // Additional re-orthogonalization procedure
        for (j = 0; j < i + 2; j++) ortho[j] = carrDot(nc, lvec[j], HC);
        for (j = 0; j < nc; j++)
        {
            for (k = 0; k < i + 2; k++) HC[j] -= lvec[k][j] * ortho[k];
        }

        if ( (i + 1) % 100 == 0 )
        {
            printf("\n  %5.1lf%%",(100.0*i)/(lm-1));
        }
    }

    free(ortho);
    free(HC);

    return lm;
}



double ground(int Niter, int Npar, int Morb, Iarray IF, Iarray strideOT,
            Iarray strideTT, Iarray Map, Iarray MapOT, Iarray MapTT,
            Cmatrix Ho, Carray Hint, Carray C)
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
        predictedIter,
        nc;

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



    nc = NC(Npar,Morb);



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

    Niter = lanczos(Niter,Npar,Morb,IF,strideOT,strideTT,Map,MapOT,MapTT,
            Ho,Hint,diag,offdiag,lvec);
    
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

    k = LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', Niter, d, e, eigvec, Niter);
    if (k != 0)
    {
        printf("\n\nERROR IN DIAGONALIZATION\n\n");
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



int main(int argc, char * argv[])
{

    /*** NUMBER OF THREADS USED - CHOOSE ACCORDINGLY TO CPU LIMITS ***/
    omp_set_num_threads(omp_get_max_threads() / 2);

    int
        i,
        nc,
        Npar,
        Morb,
        Mpos,
        nthreads;

    double
        L,
        g,
        Eo,
        Ef,
        dx,
        sum,
        end_omp,
        start_omp,
        time_used;

    Iarray
        Map,
        MapOT,
        MapTT,
        IFmat,
        NCmat,
        strideOT,
        strideTT;

    Rarray
        x;

    Carray
        C,
        HC,
        Hint;

    Cmatrix
        orbs,
        Ho;



    nthreads = omp_get_max_threads() / 2;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }



    // grid configuration of domain
    g = 100.0;
    L = 1.0;
    Mpos = 1001;
    dx = L / (Mpos - 1);

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    nc = NC(Npar,Morb);

    printf("\nNumber of particles : %3d", Npar);
    printf("\nNumber of orbitals  : %3d", Morb);
    printf("\nNumber of configurations : %d", nc);



    // Structures to handle the configurations(Fock states)
    IFmat = setupFocks(Npar,Morb);
    NCmat = setupNCmat(Npar,Morb);
    strideTT = iarrDef(nc);
    strideOT = iarrDef(nc);
    Map = OneOneMap(Npar,Morb,NCmat,IFmat);
    MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
    MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

    // Individual particle states and matrices
    x = rarrDef(Mpos);
    orbs = cmatDef(Morb,Mpos);
    Ho = cmatDef(Morb,Morb);
    Hint = carrDef(Morb*Morb*Morb*Morb);



    // Coefficients of expansion in fock basis
    C = carrDef(nc);
    HC = carrDef(nc);



    // Initialize coefficients
    sum = 0.0;
    for (i = 0; i < nc; i++)
    {
        C[i] = sin(20*((double)i)/nc) * (i%13) + ((i%8)-(i%3)) * I;
        sum = sum + creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
    }
    // normalize to 1
    for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);



    // Initialize orbitals in discretized points x
    for (i = 0; i < Mpos; i++) x[i] = - L / 2 + i * dx;
    setupOrbitals(Morb,Mpos,L,x,orbs);



    // Initialize one- and two-body matrices
    setupHo(Morb,L,Ho);
    setupHint(Morb,g,Hint);



    printf("\n\n======================================");

    time_used = 0;

    start_omp = omp_get_wtime();

    Eo = ground(nc/8,Npar,Morb,IFmat,strideOT,strideTT,Map,MapOT,MapTT,Ho,Hint,C);

    end_omp = omp_get_wtime();
    time_used = time_used + ((double) (end_omp - start_omp));

    printf("\n\nTime to find ground state(%d threads) : ",nthreads);
    printf("%.3lfs", time_used);

    // Fermi Energy
    Ef = 0;
    for (i = 1; i <= Npar/2; i++)
    {
        Ef = Ef + (2 * PI * i / L) * (2 * PI * i / L);
    }

    printf("\n\ng = %.5lf | Eo = %.5lf | Efermi = %.5lf",g,Eo/Npar,Ef/Npar);


    printf("\n\n======================================\n\n");

    free(x);

    free(C);
    free(HC);
    free(Map);

    for(i = 0; i < Morb; i++) free(orbs[i]);
    free(orbs);

    for(i = 0; i < Morb; i++) free(Ho[i]);
    free(Ho);

    free(Hint);

    free(IFmat);
    free(NCmat);

    free(strideOT);
    free(strideTT);
    free(MapOT);
    free(MapTT);

    printf("\n\nDone.\n\n");
    return 0;
}
