#include <time.h>
#include "Hdiagonal.h"

#define PI 3.141592653589793

/* icc -o exe groundStateLieb.c \
            -L${MKLROOT}/lib/intel64 \
            -lmkl_intel_lp64 \
            -lmkl_gnu_thread \
            -lmkl_core \
            -qopenmp \
            -O3 \
*/



int positiveMomentum(int M, Iarray occ)
{

/** given the occupations ('occ') in the IPS with angular momentum  l >= 0,
    compute the total momentum, where 'M' is the number of IPS with l >= 0 **/

    int
        i,
        L;

    L = 0;
    for (i = 1; i < M; i++) L = L + i * occ[i];

    return L;
}



void count_rec(int L, int l, int * vec, int vec_size, int * countAdd, int Nmax)
{

/** RECURSION TO BUILD THE ADDITIONS
    --------------------------------
    The first argument 'L' is how much need to be added yet.  The
    second 'l' is the last number added. 'vec' stack up the terms
    already used. **/

    int
        x,
        Npar,
        init;

    if (L == 0)
    {
        Npar = 0;
        for (x = 0; x < vec_size; x++) Npar = Npar + vec[x];
        if (Npar <= Nmax) *countAdd = *countAdd + 1;
        return;
    }

    if (L > l) init = l;
    else       init = L;

    for (x = init; x > 0; x--)
    {
        vec[x-1] = vec[x-1] + 1;
        count_rec(L-x,x,vec,vec_size,countAdd,Nmax);
        vec[x-1] = vec[x-1] - 1;
    }
}



int countPos(int L, int lmax, int Nmax)
{

/** 'L' is the total angular momentum yet to be subtracted.  This function
    count how many configurations using the -lmax < l < 0 IPS there are to
    provide angular momenum '-L' constrained to use up to 'Nmax' particles.
    This can be mapped to the problem of how many sums of positive integer
    numbers we can write to result in 'L'. **/

    int
        i,
        l,
        Npos,
        start_l;

    int
        * vec;

    // if there is nothing to be subtracted
    if (L == 0) return 1;

    vec = iarrDef(lmax);
    Npos = 0;

    if (L > lmax) start_l = lmax;
    else          start_l = L;

    for (l = start_l; l > 0; l--)
    {
        for (i = 0; i < lmax; i++) vec[i] = 0;
        vec[l-1] = 1;
        count_rec(L-l,l,vec,lmax,&Npos,Nmax); // recursion
    }

    free(vec);
    return Npos;
}



int nFocks(int N, int lmax, int totalL)
{

/** given 'N' particles to be arranged in (2*lmax+1) IPS corresponging to
    the angular momemtum single-particle eigenstates in increasing  order
    -lmax, -lmax+1, ... , lmax-1, lmax, RETURN the size of the multiconf.
    space with 'totalL' angular momentum.

    ALGORITHM DESCRIPTION
    =====================
    First the IPS are divided in two classes,  the  first one with angular
    momentum l >= 0, and the second one with l < 0.  Then  we  sweep  over
    all possible configurations of 'N' particles in the first class of IPS,
    using the default prescription of general hashing without restrictions.
    Each of the these configurations produce an total angular momentum 'Lp'
    Then, we use the remaining particles in the IPS with l = 0 to set them
    on the second class of IPS, with l < 0, in order to produce 'Lm', such
    that,  Lp + Lm = totalL.  The  problem can be mapped in how many there
    are to obtain the result Lp - totalL as sum of positive integer numbers
    that are less than 'lmax'. **/

    int
        i,
        M,
        N0,
        Lp,
        Nconf,
        maxIndex;

    Iarray
        occ_Lp;

    M = lmax + 1;       // Number of IPS with lmax >= l >= 0
    maxIndex = NC(N,M); // number of configurations in the space with l >= 0
    occ_Lp = iarrDef(M);

    Nconf = 0;
    for (i = 0; i < maxIndex; i++)
    {
        IndexToFock(i,N,M,occ_Lp);
        Lp = positiveMomentum(M,occ_Lp);
        N0 = occ_Lp[0]; // particles available to set in l < 0 IPS
        // Check with the available particles if there is a  configuration
        // using l < 0 IPS in order to give total angular momentum 'totalL'
        if (Lp - N0*lmax <= totalL && Lp >= totalL)
        {
            Nconf = Nconf + countPos(Lp-totalL,lmax,N0);
        }
    }

    free(occ_Lp);
    return Nconf;
}



void htBranch_rec(int L, int l, int * vec, int vec_size, int * countAdd,
        int Nmax, int stride, Iarray * ht)
{

    int
        i,
        x,
        Npar,
        init,
        confIndex;

    if (L == 0)
    {
        Npar = 0;
        // Compute number of particles required
        for (x = 0; x < vec_size; x++) Npar = Npar + vec[x];
        // Check if the number of particles required is less than
        // the available in the IPS with l = 0
        if (Npar <= Nmax)
        {
            // setup a new configuration in the hashing table
            confIndex = *countAdd + stride;
            ht[confIndex][vec_size] = Nmax - Npar;
            for (i = 0; i < vec_size; i++)
            {
                ht[confIndex][i] = vec[vec_size-1-i];
            }
            // update next configuration index in the hashing table
            *countAdd = *countAdd + 1;
        }
        return;
    }

    if (L > l) init = l;
    else       init = L;

    for (x = init; x > 0; x--)
    {
        vec[x-1] = vec[x-1] + 1;
        htBranch_rec(L-x,x,vec,vec_size,countAdd,Nmax,stride,ht);
        vec[x-1] = vec[x-1] - 1;
    }
}



int htBranch(int L, int lmax, int Nmax, int last_i, Iarray * ht)
{

/** 'L' is the total angular momentum yet to be subtracted,  obtained
    from positive momentum configuration.  'lmax' is the maximum ang.
    momentum of IPS and 'Nmax' is the remaining particles to be used.
    'last_i'  is the last configuration index used of the total mult.
    conf. space, whose the hashing table 'ht' refers to.  'last_i' is
    a offset (stride) that must be added to each possibility found in
    order to access the correct index in the hashing table 'ht'.

    RETURN the number of new configurations assembled in the 'ht' **/

    int
        i,
        l,
        Npos,
        start_l;

    int
        * vec;

    if (L == 0)
    {
        ht[last_i][lmax] = Nmax;
        for (i = 0; i < lmax; i++) ht[last_i][i] = 0;
        return 1;
    }

    vec = iarrDef(lmax);
    Npos = 0;

    if (L > lmax) start_l = lmax;
    else          start_l = L;

    for (l = start_l; l > 0; l--)
    {
        for (i = 0; i < lmax; i++) vec[i] = 0;
        vec[l-1] = 1;
        htBranch_rec(L-l,l,vec,lmax,&Npos,Nmax,last_i,ht); // recursion
    }

    free(vec);
    return Npos;
}



Iarray * assembleHT(int N, int lmax, int totalL, int mcsize)
{

/** with the total number of configurations (occupation vectors) given in
    'mcsize',  allocate a matrix  (hashing table)  with 'mcsize' rows and
    2*lmax+1  columns,  whose its elements corresponds to the occupations
    in the IPS, in ascending order of angular momentum, -lmax, ..., lmax.
    Each configuration given in a row has total angular momentum 'totalL'
    and 'N' particles.
 
    RETURN hashing table as matrix with configurations along rows **/

    int
        i,
        j,
        k,
        L,
        M,
        conf_i,
        prev_i,
        posL_mcsize;

    Iarray
        occ;

    Iarray
        * htable;

    // Allocate the hashing table. The number of IPS is 2*lmax+1
    htable = (Iarray *) malloc(mcsize*sizeof(int *));
    for (i = 0; i < mcsize; i++) htable[i] = iarrDef(2*lmax+1);

    M = lmax + 1; // Number of IPS with l >= 0
    occ = iarrDef(M); // occupations in IPS with l >= 0
    posL_mcsize = NC(N,M); // Size of multiconf. space with positive 'l' IPS

    conf_i = 0;
    for (i = 0; i < posL_mcsize; i++)
    {
        IndexToFock(i,N,M,occ);
        L = positiveMomentum(M,occ);
        if (L - occ[0]*lmax <= totalL && L >= totalL)
        {
            // assemble positive momentum SPS occupations
            prev_i = conf_i;
            // assemble negative momentum SPS occupations
            // 'conf' track the configuration index in the hashing table
            conf_i = conf_i + htBranch(L-totalL,lmax,occ[0],conf_i,htable);
            // setup the positive angular momentum IPS occupations
            for (k = prev_i; k < conf_i; k++)
            {
                for (j = lmax+1; j < 2*lmax+1; j++)
                {
                    htable[k][j] = occ[j-lmax];
                }
            }
        }
    }

    free(occ);
    return htable;
}






int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        nnz,
        Npar, // number of particles
        lmax, // number of orbitals
        totalL,
        mcSize,
        lan_it,
        nthreads;

    double
        l,
        g,
        E0,
        prevE0,
        sum,
        real,
        imag,
        time_used;

    clock_t
        start,
        end;

    Iarray
        NNZrow,
        * ht;

    Carray
        C,
        Ho;

    HConfMat
        H;

    nthreads = omp_get_max_threads() / 2;
    omp_set_num_threads(nthreads);

    if (argc != 5)
    {
        printf("\n\nERROR: Need four command line arguments");
        printf("\n\t1. Number of particles");
        printf("\n\t2. max. IPS angular momentum");
        printf("\n\t3. Total angular momentum");
        printf("\n\t4. Contact interaction strength\n\n");
        exit(EXIT_FAILURE);
    }

    // READ COMMAND LINE ARGUMENTS TO SETUP THE CONFIGURATIONAL
    // SPACE WITH FIXED TOTAL ANGULAR MOMENTUM  'totalL'  GIVEN
    // NUMBER OF PARTICLES AND THE MAXIMUM  INDIVIDUAL MOMENTUM
    // 'lmax' A PARTICLE CAN HAVE. FINALLY, GET THE INTERACTION
    // STRENGTH PARAMETER 'g'
    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&lmax);
    sscanf(argv[3],"%d",&totalL);
    sscanf(argv[4],"%lf",&g);

    printf("\nConfiguring multiconf. space structures ...\n");

    start = clock(); // trigger to measure time
    mcSize = nFocks(Npar,lmax,totalL);
    ht = assembleHT(Npar,lmax,totalL,mcSize);
    end = clock();   // stop time measure
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("\nTime to setup Conf. space : %.1lf(ms)",time_used*1000);

    // Compute the total Number of NonZero elements in Hamiltonian matrix
    // 'nnz' and also record the Number of NonZero per row in NNZrow
    NNZrow = iarrDef(mcSize);
    nnz = NNZ_PerRow(Npar,lmax,mcSize,ht,NNZrow);

    k = 0; // maximum NNZ in a same row of the matrix
    j = 0; // Self-Consistency check with the returned value 'nnz'
    for (i = 0; i < mcSize; i++)
    {
        if (k < NNZrow[i]) k = NNZrow[i];
        j = j + NNZrow[i];
    }
    if (j != nnz)
    {
        printf("\n\nFATAL ERROR : FUNCTION TO COMPUTE ");
        printf("NONZERO ELEMENTS OF SPARSE MATRIX IS WRONG\n\n");
        exit(EXIT_FAILURE);
    }

    // Configure the single-particle hamiltonian matrix.
    // Since it is diagonal  it is store in a  1D array.
    Ho = carrDef(2*lmax+1);
    for (i = 0; i < 2*lmax+1; i++)
    {
        l = 2 * PI * (i-lmax);
        Ho[i] = 0.5*l*l;
    }

    start = clock(); // trigger to measure time
    H = assembleH(Npar,lmax,mcSize,ht,Ho,g);
    end = clock();   // stop time measure
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("\nTime to setup multiconf. Hamiltonian matrix : ");
    printf("%.1lf(ms)",time_used*1000);



    // DISPLAY INFO ABOUT THE CONFIGURATIONAL SPACE ASSEMBLED



    printf("\n\n----------------------------\n");
    printf("Multiconf. space information\n");
    printf("----------------------------\n");
    printf("Total number of states (with L = %2d subjected to ",totalL);
    printf("lmax = %d) is %d",lmax,mcSize);
    printf("\nNumber of nonzero entries in H : %d ",nnz);
    printf("Sparsity = %.5lf",1.0 - ((double) nnz)/mcSize/mcSize);
    printf("\nMaximum number of nonzero entries in a same row : %d",k);
    printf("\n=====================================");
    printf("=====================================\n");



    // COMPUTE THE GROUND STATE ENERGY AND STATE



    // define an initial guess for Lanczos tridiagonal decomposition
    C = carrDef(mcSize);
    sum = 0;
    for (i = 0; i < mcSize; i++)
    {
        real = (i % 3) - 1 + 4 * (i % 5);
        imag = 3 - (i % 7) + 2 * (i % 3);
        C[i] = real + I * imag;
        sum = sum + real*real + imag*imag;
    }
    // normalize to 1
    for (i = 0; i < mcSize; i++) C[i] = C[i] / sqrt(sum);

    if (mcSize < 101)
    {
        time_used = omp_get_wtime(); // trigger to measure time
        E0 = ground(mcSize-1,mcSize,H,C);
        time_used = omp_get_wtime() - time_used; // finish time measure
    }
    else
    {
        lan_it = 100;
        time_used = omp_get_wtime(); // trigger to measure time
        E0 = ground(100,mcSize,H,C);
        prevE0 = E0;
        for (lan_it = 100; lan_it + 100 < mcSize - 1; lan_it += 100)
        {
            E0 = ground(100,mcSize,H,C);
            if ((mcSize-lan_it)*(prevE0/E0-1)/100.0 < 1E-10)
            {
                lan_it = lan_it + 100;
                break;
            }
            prevE0 = E0;
        }
        time_used = omp_get_wtime() - time_used; // finish time measure
    }

    printf("\n\nGround state energy per particle : %.10lf ",E0/Npar);
    printf("(with %d Lanczos it.)",lan_it);
    printf("\nTime to compute ground state (with %d threads) : ",nthreads);
    printf("%.1lf(s)",time_used);

    printf("\n\nDone.\n\n");

    for (i = 0; i < mcSize; i++) free(ht[i]);
    free(C);
    free(ht);
    free(Ho);
    free(NNZrow);
    freeHmat(H);

    return 0;
}
