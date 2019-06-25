#include ".h"
#include <omp.h>





void cprint(double complex z)
{

/** printf complex number in coordinates with 2 decimal digits**/

    printf("(%9.2E,%9.2E )", creal(z), cimag(z));
}





void carr_txt(char fname [], int M, Carray v)
{

/** Record a array of complex elements in a text file in a
  * suitable format to import as numpy array with python.
**/

    int
        j;

    double
        real,
        imag;

    FILE
        * data_file = fopen(fname, "w");



    if (data_file == NULL)
    {
        printf("\n\n\n\tERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < M - 1; j ++)
    {

        real = creal(v[j]);
        imag = cimag(v[j]);

        if (imag >= 0) fprintf(data_file, "(%.15E+%.15Ej)", real, imag);
        else           fprintf(data_file, "(%.15E%.15Ej)", real, imag);

        fprintf(data_file, "\n");
    }

    real = creal(v[M-1]);
    imag = cimag(v[M-1]);

    if (imag >= 0) fprintf(data_file, "(%.15E+%.15Ej)", real, imag);
    else           fprintf(data_file, "(%.15E%.15Ej)", real, imag);

    fclose(data_file);
}





void cmat_txt (char fname [], int m, int n, Cmatrix A)
{

    int
        i,
        j;

    double
        real,
        imag;

    FILE
        * f = fopen(fname, "w");



    if (f == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\n\tERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }



    for (i = 0; i < m - 1; i++)
    {

        for (j = 0; j < n; j++)
        {

            real = creal(A[i][j]);
            imag = cimag(A[i][j]);

            if (imag >= 0) fprintf(f, "(%.15E+%.15Ej) ", real, imag);
            else           fprintf(f, "(%.15E%.15Ej) ", real, imag);
        }

        fprintf(f, "\n");
    }

    for (j = 0; j < n; j++)
    {

        real = creal(A[m-1][j]);
        imag = cimag(A[m-1][j]);

        if (imag >= 0) fprintf(f, "(%.15E+%.15Ej) ", real, imag);
        else           fprintf(f, "(%.15E%.15Ej) ", real, imag);
    }

    fclose(f);
}





int fac(int n)
{   // return n !
    int n_fac = 1;
    for (int i = 1; i < n; i++) n_fac = n_fac * (i + 1);
    return n_fac;
}





int NC(int N, int M)
{   // Return # of possible configurations of N particles in M orbitals
    int i, n = 1;
    if  (M > N)
    {
        for (i = N + M - 1; i > M - 1; i --) n = n * i;
        return n / fac(N);
    }
    else
    {
        for (i = N + M - 1; i > N; i --) n = n * i;
        return n / fac(M - 1);
    }
}





void IndexToFock(int k, int N, int M, int ** NCmat, int * v)
{
    int x;

    int  i,
         m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    /* -------------------------------------------------------------
     * Put the particles in orbitals while has combinations to spend
    ----------------------------------------------------------------- */
    while ( k > 0 )
    {
        while ( k - NCmat[N][m] < 0 ) m = m - 1;
        x = k - NCmat[N][m];
        while ( x >= 0 )
        {
            v[m] = v[m] + 1; // One more particle in orbital m
            N = N - 1;       // Less one particle to setup
            k = x;
            x = x - NCmat[N][m];
        }
    }

    /* -------------------------------------------------------------
     * Put the particles in orbitals while has combinations to spend
    ----------------------------------------------------------------- */
    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}





int FockToIndex(int N, int M, int ** NCmat, int * v)
{
    int i, n;
    int k = 0;

    /* ---------------------------------------------------
     * Empty one by one orbital starting from the last one
    ------------------------------------------------------ */
    for (i = M - 1; i > 0; i--)
    {
        n = v[i]; // Number of particles in the orbital
        while (n > 0)
        {
            k = k + NCmat[N][i]; // number of combinations needed
            N = N - 1;           // decrease the number of particles
            n = n - 1;
        }
    }

    return k;
}





int ** MountNCmat(int N, int M)
{   // Matrix of all possible configurations with
    // # particles < N and M < # of orbitals
    int i,
        j;

    int ** NCmat = (int ** ) malloc((N + 1) * sizeof(int * ));

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i] = (int * ) malloc((M + 1) * sizeof(int));
    }

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i][0] = 0;
        NCmat[i][1] = 1;
        for (j = 2; j < M + 1; j++) NCmat[i][j] = NC( i , j );
    }

    return NCmat;
}





int ** MountFocks(int N, int M, int ** NCmat)
{   // All possible occupation vectors organized in rows
    int k;

    int ** ItoFock = (int **) malloc(NC(N, M) * sizeof(int *));

    for (k = 0; k < NC(N, M); k++)
    {
        ItoFock[k] = (int * ) malloc(M * sizeof(int));
        IndexToFock(k, N, M, NCmat, ItoFock[k]);
    }

    return ItoFock;
}





void JumpMapping(int N, int M, int ** NCmat, int ** IF, int * Map)
{
    int i,
        q,
        k,
        l,
        nc = NCmat[N][M],
        * v = (int *) malloc(M * sizeof(int));
    
    for (i = 0; i < nc * M * M; i++) Map[i] = 0;

    for (i = 0; i < nc; i++)
    {
        // Copy the occupation vector from C[i] coeff.
        for (q = 0; q < M; q++) v[q] = IF[i][q];

        for (k = 0; k < M; k++)
        {
            // Take one particle from k state
            if (v[k] < 1) continue;
            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                // Put one particle in l state
                v[k] -= 1;
                v[l] += 1;
                Map[i + k * nc + l * M * nc] = FockToIndex(N, M, NCmat, v);
                v[k] += 1;
                v[l] -= 1;
            }
        }
    }

    free(v);
}





void TwiceJumpMapping(int N, int M, int ** NCmat, int ** IF, int * Map)
{
    int i,
        q,
        k,
        s,
        l,
        ind,
        M2 = M*M,
        M3 = M*M2,
        M4 = M*M3,
        nc = NCmat[N][M],
        * v = (int *) malloc(M * sizeof(int));
    
    for (i = 0; i < nc * M4; i++) Map[i] = -1;

    for (i = 0; i < nc; i++)
    {
        // Copy the occupation vector from C[i] coeff.
        for (q = 0; q < M; q++) v[q] = IF[i][q];

        for (k = 0; k < M; k++)
        {
            // Take one particle from k state
            if (v[k] < 1) continue;

            for (s = 0; s < M; s++)
            {
                // take one particle from s state
                if (v[s] < 1) continue;

                for (l = 0; l < M; l++)
                {

                    // put one particle in l state
                    for (q = 0; q < M; q++)
                    {

                        // Put one particle in q state
                        ind = i + k*nc + s*M*nc + l*M2*nc + q*M3*nc;
                        v[k] -= 1;
                        v[s] -= 1;
                        v[l] += 1;
                        v[q] += 1;
                        Map[ind] = FockToIndex(N, M, NCmat, v);
                        v[k] += 1;
                        v[s] += 1;
                        v[l] -= 1;
                        v[q] -= 1;
                    }
                }
            }
        }
    }

    free(v);
}










/* ========================================================================
 *
 *                           <   a*_k   a_l   >
 *                    -------------------------------
 *
 * Once defined a set of Single-Particle Wave Functions (SPWF) a many
 * body state  can be expanded  in a  Occupation Number Configuration
 * Basis (ONCB) whose vector are also named Fock states. The one body
 * density matrix is known as the expected value of 1 creation  and 1
 * annihilation operators for a given many-body state.  Use the basis
 * to express the state and then compute using its coefficients (Cj).
 *
 * ======================================================================== */





void OBrho(int N, int M, int ** NCmat, int ** IF, Carray C, Cmatrix rho)
{

    int i,
        j,
        k,
        l,
        q,
        nc,
        * v;

    double
        mod2;

    double complex
        RHO;

    nc = NCmat[N][M];

    for (k = 0; k < M; k++)
    {

        // Diagonal elements
        // ------------------------------------------------------------------

        RHO = 0;

        #pragma omp parallel shared(k, nc, C, IF) private(i, mod2) \
        reduction(+:RHO)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO = RHO + mod2 * IF[i][k];
            }
        }

        rho[k][k] = RHO;
        // ------------------------------------------------------------------



        // Off-diagonal elements (different terms from operators)
        // ------------------------------------------------------------------
        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            #pragma omp parallel shared(l,k,M,N,nc,C,NCmat,IF) \
            private(i, j, q, v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for schedule(static)
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 1) continue;
                    for (q = 0; q < M; q++) v[q] = IF[i][q];
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    v[k] += 1;
                    v[l] -= 1;
                    RHO += conj(C[i]) * C[j] * sqrt((double)(v[l]+1) * v[k]);
                }

                free(v); // Each thread release its vector
            }

            rho[k][l] = RHO;
            rho[l][k] = conj(RHO); // hermiticity
        }
        // ------------------------------------------------------------------
    }

}





void OBrho_X(int N, int M, int * Map1, int ** NCmat, int ** IF, Carray C, Cmatrix rho)
{

    int i,
        j,
        k,
        l,
        q,
        nc,
        vl,
        vk;

    double
        mod2;

    double complex
        RHO;

    nc = NCmat[N][M];

    for (k = 0; k < M; k++)
    {

        // Diagonal elements
        // ------------------------------------------------------------------

        RHO = 0;

        #pragma omp parallel shared(k, nc, C, IF) private(i, mod2) \
        reduction(+:RHO)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO = RHO + mod2 * IF[i][k];
            }
        }

        rho[k][k] = RHO;
        // ------------------------------------------------------------------



        // Off-diagonal elements (different terms from operators)
        // ------------------------------------------------------------------
        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            #pragma omp parallel shared(l,k,M,N,nc,C,NCmat,IF) \
            private(i, j, q, vk, vl) reduction(+:RHO)
            {

                #pragma omp for schedule(static)
                for (i = 0; i < nc; i++)
                {
                    j = Map1[i + k * nc + l * M * nc];
                    if (j < 0) continue;
                    vk = IF[i][k];
                    vl = IF[i][l];
                    /* if (IF[i][k] < 1) continue;
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    v[k] += 1;
                    v[l] -= 1; */
                    RHO += conj(C[i]) * C[j] * sqrt((double)(vl+1) * vk);
                }

            }

            rho[k][l] = RHO;
            rho[l][k] = conj(RHO); // hermiticity
        }
        // ------------------------------------------------------------------
    }

}










/* ========================================================================
 *
 *                    <   a*_k   a*_s   a_q   a_l   >
 *                    -------------------------------
 *
 * Once defined a set of Single-Particle Wave Functions (SPWF) a many
 * body state  can be expanded  in a  Occupation Number Configuration
 * Basis (ONCB) whose vector are also named Fock states. The two body
 * density matrix is known as the expected value of 2 creation  and 2
 * annihilation operators for a given many-body state.  Use the basis
 * to express the state and then compute using its coefficients (Cj).
 *
 * ======================================================================== */





void TBrho(int N, int M, int ** NCmat, int ** IF, Carray C, Carray rho)
{

    int i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        t,
        nc,
        M2,
        M3,
        * v;
    
    double
        mod2,   // |Cj| ^ 2
        sqrtOf; // Factors from the action of creation/annihilation

    double complex
        RHO;





    // Auxiliar to memory access
    M2 = M * M;
    M3 = M * M * M;

    nc = NCmat[N][M];



    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {

        RHO = 0;

        #pragma omp parallel for private(i,mod2) reduction(+:RHO)
        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO  = RHO + mod2 * IF[i][k] * (IF[i][k] - 1);
        }

        rho[k + M * k + M2 * k + M3 * k] = RHO;
    }



    /* ---------------------------------------------
     * Rule 2: Creation on k s / Annihilation on k s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {

            RHO = 0;

            #pragma omp parallel for private(i,mod2) reduction(+:RHO)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[i][k] * IF[i][s];
            }

            // commutation of bosonic operators is used
            // to fill elements by exchange  of indexes
            rho[k + s * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + k * M3] = RHO;
            rho[k + s * M + s * M2 + k * M3] = RHO;
        }
    }



    /* ---------------------------------------------
     * Rule 3: Creation on k k / Annihilation on q q
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = sqrt((double)(v[k]-1)*v[k]*(v[q]+1)*(v[q]+2));
                    v[k] -= 2;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            // Use 2-index-'hermiticity'
            rho[k + k * M + q * M2 + q * M3] = RHO;
            rho[q + q * M + k * M2 + k * M3] = conj(RHO);
        }
    }



    /* ---------------------------------------------
     * Rule 4: Creation on k k / Annihilation on k l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            rho[k + k * M + k * M2 + l * M3] = RHO;
            rho[k + k * M + l * M2 + k * M3] = RHO;
            rho[l + k * M + k * M2 + k * M3] = conj(RHO);
            rho[k + l * M + k * M2 + k * M3] = conj(RHO);
        }
    }



    /* ---------------------------------------------
     * Rule 5: Creation on k s / Annihilation on s s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = v[s] * sqrt((double)v[k]*(v[s]+1));
                    v[k] -= 1;
                    v[s] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            rho[k + s * M + s * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + s * M3] = RHO;
            rho[s + s * M + s * M2 + k * M3] = conj(RHO);
            rho[s + s * M + k * M2 + s * M3] = conj(RHO);
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (k = q + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (l = q + 1; l < M; l++)
        {
            for (k = l + 1; k < M; k++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.0: Creation on k s / Annihilation on s l (s < k < l)
    ------------------------------------------------------------------- */
    for (s = 0; s < M; s++)
    {
        for (k = s + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.1: Creation on k s / Annihilation on s l (k < s < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (l = s + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.2: Creation on k s / Annihilation on s l (k < l < s)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            for (s = l + 1; s < M; s++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k]*(v[l]+1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* ---------------------------------------------
     * Rule 8: Creation on k s / Annihilation on q l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = 0; s < M; s++)
        {
            if (s == k) continue;

            for (q = 0; q < M; q++)
            {
                if (q == s || q == k) continue;

                for (l = 0; l < M; l ++)
                {

                    RHO = 0;

                    if (l == k || l == s || l == q) continue;

                    #pragma omp parallel private(i,j,t,sqrtOf,v) \
                    reduction(+:RHO)
                    {
                        v = (int *) malloc(M * sizeof(int));
                        #pragma omp for
                        for (i = 0; i < nc; i++)
                        {
                            if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                            for (t = 0; t < M; t++) v[t] = IF[i][t];
                            sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                            v[k] -= 1;
                            v[s] -= 1;
                            v[q] += 1;
                            v[l] += 1;
                            j = FockToIndex(N, M, NCmat, v);
                            RHO += conj(C[i]) * C[j] * sqrtOf;
                        }
                        free(v);
                    }

                    rho[k + s * M + q * M2 + l * M3] = RHO;
                }   // Finish l loop
            }       // Finish q loop
        }           // Finish s loop
    }               // Finish k loop


    /*       ------------------- END OF ROUTINE -------------------       */
}





void TBrho_X(int N, int M, int * Map1, int * Map2, int ** NCmat, int ** IF, Carray C, Carray rho)
{

    int i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        t,
        nc,
        M2,
        M3,
        vk,
        vl,
        vq,
        vs,
        * v;
    
    double
        mod2,   // |Cj| ^ 2
        sqrtOf; // Factors from the action of creation/annihilation

    double complex
        RHO;





    // Auxiliar to memory access
    M2 = M * M;
    M3 = M * M * M;

    nc = NCmat[N][M];



    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {

        RHO = 0;

        #pragma omp parallel for private(i,mod2) reduction(+:RHO)
        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO  = RHO + mod2 * IF[i][k] * (IF[i][k] - 1);
        }

        rho[k + M * k + M2 * k + M3 * k] = RHO;
    }



    /* ---------------------------------------------
     * Rule 2: Creation on k s / Annihilation on k s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {

            RHO = 0;

            #pragma omp parallel for private(i,mod2) reduction(+:RHO)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[i][k] * IF[i][s];
            }

            // commutation of bosonic operators is used
            // to fill elements by exchange  of indexes
            rho[k + s * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + k * M3] = RHO;
            rho[k + s * M + s * M2 + k * M3] = RHO;
        }
    }



    /* ---------------------------------------------
     * Rule 3: Creation on k k / Annihilation on q q
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    j = i + k*nc + k*M*nc + q*M2*nc + q*M3*nc;
                    if (Map2[j] < 0) continue;
                    j = Map2[j];
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = sqrt((double)(v[k]-1)*v[k]*(v[q]+1)*(v[q]+2));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            // Use 2-index-'hermiticity'
            rho[k + k * M + q * M2 + q * M3] = RHO;
            rho[q + q * M + k * M2 + k * M3] = conj(RHO);
        }
    }



    /* ---------------------------------------------
     * Rule 4: Creation on k k / Annihilation on k l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    j = i + k*nc + l*M*nc;
                    j = Map1[j];
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            rho[k + k * M + k * M2 + l * M3] = RHO;
            rho[k + k * M + l * M2 + k * M3] = RHO;
            rho[l + k * M + k * M2 + k * M3] = conj(RHO);
            rho[k + l * M + k * M2 + k * M3] = conj(RHO);
        }
    }



    /* ---------------------------------------------
     * Rule 5: Creation on k s / Annihilation on s s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {

            RHO = 0;

            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));

                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    j = i + k*nc + s*M*nc;
                    if (Map1[j] < 0) continue;
                    j = Map1[j];
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = v[s] * sqrt((double)v[k]*(v[s]+1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }

                free(v);
            }

            rho[k + s * M + s * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + s * M3] = RHO;
            rho[s + s * M + s * M2 + k * M3] = conj(RHO);
            rho[s + s * M + k * M2 + s * M3] = conj(RHO);
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        j = i + k*nc + k*M*nc + l*M2*nc + q*M3*nc;
                        if (Map2[j] < 0) continue;
                        j = Map2[j];
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (k = q + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        j = i + k*nc + k*M*nc + l*M2*nc + q*M3*nc;
                        if (Map2[j] < 0) continue;
                        j = Map2[j];
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (l = q + 1; l < M; l++)
        {
            for (k = l + 1; k < M; k++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        j = i + k*nc + k*M*nc + l*M2*nc + q*M3*nc;
                        if (Map2[j] < 0) continue;
                        j = Map2[j];
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = conj(RHO);
                rho[q + l * M + k * M2 + k * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.0: Creation on k s / Annihilation on s l (s < k < l)
    ------------------------------------------------------------------- */
    for (s = 0; s < M; s++)
    {
        for (k = s + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        j = i + k*nc + l*M*nc;
                        if (Map1[j] < 0) continue;
                        j = Map1[j];
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.1: Creation on k s / Annihilation on s l (k < s < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (l = s + 1; l < M; l++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));

                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        j = i + k*nc + l*M*nc;
                        if (Map1[j] < 0) continue;
                        j = Map1[j];
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* -----------------------------------------------------------
     * Rule 7.2: Creation on k s / Annihilation on s l (k < l < s)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            for (s = l + 1; s < M; s++)
            {

                RHO = 0;

                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        j = i + k*nc + l*M*nc;
                        if (Map1[j] < 0) continue;
                        j = Map1[j];
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt((double)v[k]*(v[l]+1));
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + s * M2 + k * M3] = conj(RHO);
                rho[s + l * M + k * M2 + s * M3] = conj(RHO);
                rho[l + s * M + k * M2 + s * M3] = conj(RHO);
            }
        }
    }



    /* ---------------------------------------------
     * Rule 8: Creation on k s / Annihilation on q l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (q = 0; q < M; q++)
            {
                if (q == s || q == k) continue;

                for (l = q + 1; l < M; l ++)
                {

                    RHO = 0;

                    if (l == k || l == s) continue;

                    #pragma omp parallel private(i,j,t,sqrtOf,v) \
                    reduction(+:RHO)
                    {
                        v = (int *) malloc(M * sizeof(int));
                        #pragma omp for
                        for (i = 0; i < nc; i++)
                        {
                            j = i + k*nc + s*M*nc + l*M2*nc + q*M3*nc;
                            if (Map2[j] < 0) continue;
                            j = Map2[j];
                            for (t = 0; t < M; t++) v[t] = IF[i][t];
                            sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                            RHO += conj(C[i]) * C[j] * sqrtOf;
                        }
                        free(v);
                    }

                    rho[k + s * M + q * M2 + l * M3] = RHO;
                    rho[s + k * M + q * M2 + l * M3] = RHO;
                    rho[k + s * M + l * M2 + q * M3] = RHO;
                    rho[s + k * M + l * M2 + q * M3] = RHO;
                }   // Finish l loop
            }       // Finish q loop
        }           // Finish s loop
    }               // Finish k loop


    /*       ------------------- END OF ROUTINE -------------------       */
}










int main(int argc, char * argv[])
{

    omp_set_num_threads(2);


    int
        i,
        j,
        nc,
        Npar,
        Morb,
        * Map1,
        * Map2,
        ** IF,
        ** NCmat;

    double
        sum,
        start,
        time_used;

    Carray
        C,
        rho2,
        rho2_X;

    Cmatrix
        rho1_X,
        rho1;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    nc = NC(Npar,Morb);

    rho1 = (doublec **) malloc(Morb * sizeof(doublec *));
    for (i = 0; i < Morb; i++)
    {
        rho1[i] = (doublec *) malloc(Morb * sizeof(doublec));
    }

    rho1_X = (doublec **) malloc(Morb * sizeof(doublec *));
    for (i = 0; i < Morb; i++)
    {
        rho1_X[i] = (doublec *) malloc(Morb * sizeof(doublec));
    }

    rho2 = (doublec * ) malloc(Morb*Morb*Morb*Morb*sizeof(doublec));
    rho2_X = (doublec * ) malloc(Morb*Morb*Morb*Morb*sizeof(doublec));

    NCmat = MountNCmat(Npar,Morb);
    IF = MountFocks(Npar,Morb,NCmat);

    C = (doublec *) malloc(nc * sizeof(doublec));
    sum = 0;
    for (i = 0; i < nc; i++)
    {
        C[i] = sin( 30 * ((double) i) / nc);
        sum = sum + creal(C[i]) * creal(C[i]) +  cimag(C[i]) * cimag(C[i]);
    }
    for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);

    Map1 = (int *) malloc(nc * Morb * Morb * sizeof(int));
    JumpMapping(Npar,Morb,NCmat,IF,Map1);
    Map2 = (int *) malloc(nc*Morb*Morb*Morb*Morb*sizeof(int));
    TwiceJumpMapping(Npar,Morb,NCmat,IF,Map2);










    printf("\n\nNumber of particles : %3d", Npar);
    printf(  "\nNumber of orbitals  : %3d", Morb);
    printf(  "\nNumber of configurations : %d", nc);

    printf("\n\n=============================================\n\n");

    printf("All configurations :\n");
    for (i = 0; i < nc; i++)
    {
        printf("\n%8d  [", i);
        for (j = 0; j < Morb - 1; j++) printf(" %3d |", IF[i][j]);
        printf(" %3d ]", IF[i][Morb-1]);
        //printf("  %5d ", Map1[i + 1*nc + (Morb-1)*Morb*nc]);
        printf(" || C(%d) = ",FockToIndex(Npar,Morb,NCmat,IF[i]));
        printf("%.1lf + %.1lfi", creal(C[i]), cimag(C[i]));
    }





    start = omp_get_wtime();
    for (i = 0; i < 10; i++) OBrho(Npar,Morb,NCmat,IF,C,rho1);
    time_used = (double) (omp_get_wtime() - start) / 10;
    printf("\n\nTime to setup rho1 : %.1lfms", time_used * 1000);

    start = omp_get_wtime();
    for (i = 0; i < 10; i++) OBrho_X(Npar,Morb,Map1,NCmat,IF,C,rho1_X);
    time_used = (double) (omp_get_wtime() - start) / 10;
    printf("\n\nTime to setup rho1_X : %.1lfms", time_used * 1000);

    start = omp_get_wtime();
    for (i = 0; i < 5; i++) TBrho(Npar,Morb,NCmat,IF,C,rho2);
    time_used = (double) (omp_get_wtime() - start) / 5;
    printf("\n\nTime to setup rho2 : %.1lfms", time_used * 1000);

    start = omp_get_wtime();
    for (i = 0; i < 5; i++) TBrho_X(Npar,Morb,Map1,Map2,NCmat,IF,C,rho2_X);
    time_used = (double) (omp_get_wtime() - start) / 5;
    printf("\n\nTime to setup rho2_X : %.1lfms", time_used * 1000);

    carr_txt("rho2.dat",Morb*Morb*Morb*Morb,rho2);
    carr_txt("rho2_X.dat",Morb*Morb*Morb*Morb,rho2_X);
    cmat_txt("rho1.dat",Morb,Morb,rho1);
    cmat_txt("rho1_X.dat",Morb,Morb,rho1_X);

    free(rho2);
    free(C);

    for (i = 0; i < nc; i++) free(IF[i]);
    free(IF);

    for (i = 0; i < Npar + 1; i++) free(NCmat[i]);
    free(NCmat);

    printf("\n\n");
    return 0;
}
