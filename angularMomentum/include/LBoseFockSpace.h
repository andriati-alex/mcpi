#ifndef _LBoseFockSpace_h
#define _LBoseFockSpace_h

#include "BoseFockSpace.h"



int positiveMomentum(int M, Iarray occ)
{

/** given the occupations 'occ' in the IPS with  angular  momentum  l >= 0,
    compute the total momentum, where 'M' is the number of IPS with l >= 0 **/

    int
        i,
        L;

    L = 0;
    for (i = 1; i < M; i++) L = L + i * occ[i];

    return L;
}



void BCount_rec(int L, int l, int * vec, int vec_size, int * countAdd, int Nav)
{

/** RECURSION TO BUILD THE ADDITIONS - ALL WAYS TO GET 'L' ADDING
    POSITIVE INTEGER NUMBERS WITH CONSTRAINTS
    =============================================================
    The first argument 'L' is how much has to be added yet. The second
    'l' is the last number added. 'vec' track the repetitions in  each
    number used. 'Nav' constrain the possible combination found to use
    only up to 'Nav' numbers in the sum,  corresponding  to the number
    of particles available. **/

    int
        x,
        Nreq,
        init;

    if (L == 0)
    {
        Nreq = 0;
        // Sum the repetitions(occupations) to obtain
        // the number of particles required
        for (x = 0; x < vec_size; x++) Nreq = Nreq + vec[x];
        // The result is only valid if 'Nreq' particles
        // are available from the IPS with l = 0
        if (Nreq <= Nav) *countAdd = *countAdd + 1;
        return;
    }

    if (L > l) init = l;
    else       init = L;

    for (x = init; x > 0; x--)
    {
        vec[x-1] = vec[x-1] + 1;
        BCount_rec(L-x,x,vec,vec_size,countAdd,Nav);
        vec[x-1] = vec[x-1] - 1;
    }
}



int BCountPos(int L, int lmax, int Nav)
{

/** 'L' is the total angular momentum yet to be subtracted.  This function
    counts how many configurations using the -lmax <= l < 0 IPS there  are
    to provide angular momenum '-L' constrained to use up to 'Nav' particles.
    This can be mapped on the problem of how many ways we can add positive
    integer numbers <= lmax which yield 'L', using up to 'Nav' numbers.
    This is solved by the recursive function. **/

    int
        i,
        l,
        Npos,
        start_l;

    int
        * vec;

    // if there is nothing to be subtracted
    if (L == 0) return 1;

    vec = iarrDef(lmax); // occupation numbers / repetition of number l
    Npos = 0;

    // Select the orbital to start
    if (L > lmax) start_l = lmax;
    else          start_l = L;

    for (l = start_l; l > 0; l--)
    {
        for (i = 0; i < lmax; i++) vec[i] = 0;
        vec[l-1] = 1;
        BCount_rec(L-l,l,vec,lmax,&Npos,Nav); // recursion
    }

    free(vec);
    return Npos;
}



int BFixedMom_mcsize(int N, int lmax, int totalL)
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
    are to obtain the result  (Lp - totalL)  as  sum  of  positive integer
    numbers that are <= 'lmax' restricted to used the particles  available
    in the IPS with l = 0. **/

    int
        i,
        M,
        N0,
        Lp,
        Nconf,
        maxIndex,
        MemReq;

    Iarray
        occ_Lp;

    M = lmax + 1;       // Number of IPS with lmax >= l >= 0
    maxIndex = NC(N,M); // number of configurations in the space with l >= 0
    occ_Lp = iarrDef(M);

    Nconf = 0;
    for (i = 0; i < maxIndex; i++)
    {
        indexToConfig(i,N,M,occ_Lp);
        Lp = positiveMomentum(M,occ_Lp);
        N0 = occ_Lp[0]; // particles available to set in l < 0 IPS
        // Check with the available particles if there is a  configuration
        // using l < 0 IPS in order to give total angular momentum 'totalL'
        if (Lp - N0*lmax <= totalL && Lp >= totalL)
        {
            Nconf = Nconf + BCountPos(Lp-totalL,lmax,N0);
        }

        // estimate the memory required for the Hashing table
        MemReq = Nconf * (2 * lmax + 1) * sizeof(int);
        if (MemReq > MEMORY_TOL)
        {
            printf("\n\nPROCESS ABORTED : When computing the size of the ");
            printf("multiconfg. space, the approximated memory demanded ");
            printf("to set the hashing table for %d particles ",N);
            printf("subject to lmax = %d is ",lmax);
            printf("%.1lf\n",((double) MEMORY_TOL)/1E9);
            printf("Stopped at function 'BFixedMom_mcsize'\n\n");
            exit(EXIT_FAILURE);
        }
    }

    free(occ_Lp);
    return Nconf;
}



void BNegMomOcc_rec(int L, int l, int * vec, int vec_size, int * countAdd,
                    int Nav, int stride, Iarray * ht)
{

/** This routine perform essentially the same task described in 'BCount_rec'
    BUT it additionally record a valid solution for the occupations in  the
    hashing table 'ht'. See 'BCount_rec' description above **/

    int
        i,
        x,
        Nreq,
        init,
        confIndex;

    // 'vec_size' is just 'lmax'
    if (L == 0)
    {
        Nreq = 0;
        // Compute number of particles required from 'vec' of occupations
        for (x = 0; x < vec_size; x++) Nreq = Nreq + vec[x];
        // Check if the number of particles required is less  than
        // the available in the IPS with l = 0   (in positive case
        // the solution is copied as negative momentum occupations
        // in the hashing table)
        if (Nreq <= Nav)
        {
            // setup a new configuration in the hashing table
            confIndex = *countAdd + stride;
            // Set the num. of particles in l = 0 after using 'Nreq'
            ht[confIndex][vec_size] = Nav - Nreq;
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
        BNegMomOcc_rec(L-x,x,vec,vec_size,countAdd,Nav,stride,ht);
        vec[x-1] = vec[x-1] - 1;
    }
}



int BNegMomOcc(int L, int lmax, int Nav, int last_i, Iarray * ht)
{

/** 'L' is the total angular momentum yet to be subtracted,  obtained
    from positive momentum configuration.  'lmax' is the maximum ang.
    momentum of IPS and 'Nav' are the available particles to be used.
    'last_i'  is the last configuration index used of the total mult.
    conf. space, whose the hashing table 'ht' refers to.  'last_i' is
    a offset (stride) that must be added to each possibility found in
    order to access the correct index in the hashing table 'ht'.

    The problem is solve by mapping in the algorithm to find  all the
    possibile ways to yield 'L' using addition of positive  integers,
    with the maximum number used limited by 'lmax' and constrained to
    use up to 'Nav' numbers.  Finally,  these  number only have to be
    interpreted inverting their signs to map them in -lmax , ... , -1

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
        // There is nothing to be subtracted in momenum then set
        // the available particles in l = 0 IPS and none in l < 0 IPS
        ht[last_i][lmax] = Nav;
        for (i = 0; i < lmax; i++) ht[last_i][i] = 0;
        return 1;
    }

    // the indexes of 'vec' depicted positive numbers being summed
    // and its entries the number of appearances in the addition
    vec = iarrDef(lmax);
    Npos = 0;

    if (L > lmax) start_l = lmax;
    else          start_l = L;

    for (l = start_l; l > 0; l--)
    {
        for (i = 0; i < lmax; i++) vec[i] = 0;
        vec[l-1] = 1;
        BNegMomOcc_rec(L-l,l,vec,lmax,&Npos,Nav,last_i,ht); // recursion
    }

    free(vec);
    return Npos;
}



Iarray * BAssembleHT(int N, int lmax, int totalL, int mcsize)
{

/** With the total number of configurations (occupation vectors) given in
    'mcsize',  allocate a matrix  (hashing table)  with 'mcsize' rows and
    2*lmax+1  columns,  whose its elements corresponds to the occupations
    in the IPS, in ascending order of angular momentum, -lmax, ..., lmax.
    Each configuration given in a row has total angular momentum 'totalL'
    and 'N' particles(bosons).
 
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

    if (mcsize == 0)
    {
        printf("\n\nWARNING : empty hashing table allocation requested ");
        printf("with %d particles, |l_max| = %d and L = %d\n\n",N,lmax,totalL);
    }

    // Allocate the hashing table. The number of IPS is 2*lmax+1
    htable = (Iarray *) malloc(mcsize*sizeof(int *));
    for (i = 0; i < mcsize; i++) htable[i] = iarrDef(2*lmax+1);

    M = lmax + 1;           // Number of IPS with l >= 0
    occ = iarrDef(M);       // occupations in IPS with l >= 0
    posL_mcsize = NC(N,M);  // Size of multiconf. space restricted to l >= 0

    conf_i = 0;
    for (i = 0; i < posL_mcsize; i++)
    {
        indexToConfig(i,N,M,occ);
        L = positiveMomentum(M,occ);
        // Check if with the available particles from l = 0 state (occ[0])
        // the maximum value of momentum that can be subtracted is enough
        // to produce the desired 'totalL' momentum
        if (L - occ[0]*lmax <= totalL && L >= totalL)
        {
            prev_i = conf_i;
            // assemble negative momentum IPS occupations 'conf_i'
            // track the configuration index in the hashing table
            conf_i = conf_i + BNegMomOcc(L-totalL,lmax,occ[0],conf_i,htable);
            // setup the positive angular momentum IPS occupations
            // commom to all possibilities found for  l  <  0  IPS
            // Note that l = 0 IPS was already set in 'BNegMomOcc'
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



int BgetIndex(int lmax, int ht_size, Iarray * ht, Iarray occ)
{

/** Search in the hashing table a configuration and RETURN its index
    The algorithm assumes an ordering implied by the construction of
    the hashing table(ht). In the way the ht was set, First one need
    to compare the occupations in positive momentum IPS,  where  the
    ordering is increasing according to number of particles  in  IPS
    with higher momentum. Second, one need to compare occupations in
    negative momentum IPS. This step involves the recursion part, as
    such, the states with larger occupations in most negative moment
    IPS appear first and later the ones with larger occupations near
    the l = 0 IPS **/

    int
        i,
        n,
        lower,
        upper;

    if (lmax == 0) return 0; // only one single particle state

    lower = 0;
    upper = ht_size;
    while (1)
    {
        n = (upper + lower) / 2;

        i = 2*lmax;
        while(occ[i] == ht[n][i])
        {
            i = i - 1;
            if (i == lmax) break;
        }
        if (i > lmax)
        {
            if (occ[i] > ht[n][i]) lower = n;
            else upper = n;
            continue;
        }

        i = 0;
        while(occ[i] == ht[n][i])
        {
            i = i + 1;
            if (i == lmax) break;
        }
        if (i < lmax)
        {
            if (occ[i] < ht[n][i]) lower = n;
            else upper = n;
            continue;
        }

        return n;
    }

}



#endif
