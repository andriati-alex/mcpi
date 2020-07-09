#ifndef _LFermiFockSpace_h
#define _LFermiFockSpace_h

#include "FermiFockSpace.h"



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



void FCount_rec(int L, int l, int * vec, int vec_size, int * countAdd, int Nav)
{

/** RECURSION TO BUILD THE ADDITIONS FOR FERMIONS
    ---------------------------------------------
    The first argument 'L' is how much need to be added yet. The
    second 'l' is the last number added. 'vec' stack up the  IPS
    already used. 'Nav' is the number of particles available  to
    use. Thus, as this is for fermions any outcome must strictly
    use 'Nav'.  This is the recursive implementation analogue to
    the problem of how many ways there are to add positive integer
    number which yeild the same result 'L' constrained to use 'Nav'
    terms in the sum and with the max. term in the sum given by 'l'  **/

    int
        x,
        Nreq,
        init;

    if (L == 0)
    {
        Nreq = 0;
        // Compute number of particles required from 'vec' of occ.
        for (x = 0; x < vec_size; x++) Nreq = Nreq + vec[x];
        // FOR FERMIONS the number of available  particles  MUST
        // be used, then for a valid config. it must be strictly
        // equal to the number required 'Nreq'
        if (Nreq == Nav) *countAdd = *countAdd + 1;
        return;
    }

    // DIFFERENCE FROM BOSONIC CASE. Since 'l' is the last IPS
    // occupied, the it cannot start again from 'l'
    if (L > l-1) init = l-1;
    else         init = L;

    // The maximum momentum generated is an Arithmetic Prog.
    for (x = init; (x*(x+1))/2 >= L; x--)
    {
        vec[x-1] = vec[x-1] + 1;
        FCount_rec(L-x,x,vec,vec_size,countAdd,Nav);
        vec[x-1] = vec[x-1] - 1;
    }
}



int FCountPos(int L, int lmax, int Npar)
{

/** 'L' is the total angular momentum yet to be subtracted.   This function
    count how many configurations using the -lmax < l < 0  IPS there are to
    provide angular momenum '-L' constrained to use up to 'Npar' particles.
    This can be mapped to the problem of how many sums of positive  integer
    numbers we can write to result in 'L', constrained to used 'Npar' terms
    and bounded by 'lmax'. **/

    int
        i,
        l,
        Npos,
        start_l;

    int
        * vec;

    // if there is nothing to be subtracted, i.e L = 0 already
    // ADDITIONAL CHECK FOR FERMIONS
    if (L == 0)
    {
        if (Npar > 0) return 0;
        else          return 1;
    }

    vec = iarrDef(lmax); // 'vec' of occupations in l < 0 IPS
    Npos = 0; // Number of possibilities found. Update during recursion

    if (L > lmax) start_l = lmax;
    else          start_l = L;

    // For fermions the maximum momentum is obtained by an Arithmetic Prog.
    for (l = start_l; (l*(1+l))/2 >= L; l--)
    {
        for (i = 0; i < lmax; i++) vec[i] = 0;
        vec[l-1] = 1;
        FCount_rec(L-l,l,vec,lmax,&Npos,Npar); // recursion
    }

    free(vec);
    return Npos;
}



int FFixedMom_mcsize(int N, int lmax, int totalL)
{

/** given 'N'  fermions to be arranged in (2*lmax+1) IPS corresponging to
    the angular momemtum single-particle eigenstates in increasing  order
    -lmax, -lmax+1, ... , lmax-1, lmax, RETURN the size of the multiconf.
    space with 'totalL' momentum.

    ALGORITHM DESCRIPTION
    =====================
    First the IPS are divided in two classes,  the  first one with angular
    momentum l >= 0, and the second one with l < 0.  Then  we  sweep  over
    all possible configurations of N-N0 particles in the first class of IPS,
    using the default prescription of general hashing without restrictions.
    Each of the these configurations produce an total angular momentum 'Lp'
    Then,  we use the remaining particles  'N0'  to set them on the  second
    class of IPS, with l < 0, in order to produce 'Lm', such that,  Lp + Lm
    = totalL. The problem can be mapped in how many there are to obtain the
    result Lp - totalL as sum of positive integer numbers that are less than
    'lmax' with the proper constraints for fermions. **/

    int
        i,
        M,
        N0,
        Lp,
        Nconf,
        posNC;

    Iarray
        occ_Lp;

    M = lmax + 1; // Number of IPS with l such that 0 <= l <= lmax
    occ_Lp = iarrDef(M);

    Nconf = 0;
    for (N0 = 0; N0 <= N; N0++)
    {
        // N0  is the number of particles reserved
        // to use in negative angular momentum IPS
        posNC = NC(N-N0,M);
        // posNC is the size of the spanned multiconfigurational
        // space with positive angular momentum IPS
        for (i = 0; i < posNC; i++)
        {
            IndexToFock(i,N-N0,M,occ_Lp);
            Lp = positiveMomentum(M,occ_Lp);
            // Use the remaining particles  'N0' to search for occupations
            // in negative momentum IPS that combined with the occupations
            // in positive momentum  'occ_Lp' give total momentum 'totalL'
            // Note that the maximum contribution of negative momentum is
            // given by an AP starting in -lmax
            if (Lp - (N0*(2*lmax+1-N0))/2 <= totalL && Lp >= totalL)
            {
                Nconf = Nconf + FCountPos(Lp-totalL,lmax,N0);
            }
        }
    }

    free(occ_Lp);
    return Nconf;
}



void FNegMomOcc_rec(int L, int l, int * vec, int vec_size, int * countAdd,
                    int Navailable, int stride, Iarray * ht)
{

    int
        i,
        x,
        Nreq,
        init,
        confIndex;

    if (L == 0)
    {
        Nreq = 0;
        // Compute number of particles required from 'vec' of occ.
        for (x = 0; x < vec_size; x++) Nreq = Nreq + vec[x];
        // Check if the number of particles required is the
        // same as the one available
        if (Nreq == Navailable)
        {
            // setup a new configuration in the hashing table
            confIndex = *countAdd + stride;
            for (i = 0; i < vec_size; i++)
            {
                ht[confIndex][i] = vec[vec_size-1-i];
            }
            // update next configuration index in the hashing table
            *countAdd = *countAdd + 1;
        }
        return;
    }

    // DIFFERENCE FROM BOSONIC CASE. Since 'l' is the last IPS
    // occupied, the it cannot start again from 'l'
    if (L > l-1) init = l-1;
    else         init = L;

    // The maximum momentum generated is an Arithmetic Prog. for fermions
    for (x = init; (x*(x+1))/2 >= L; x--)
    {
        vec[x-1] = vec[x-1] + 1;
        FNegMomOcc_rec(L-x,x,vec,vec_size,countAdd,Navailable,stride,ht);
        vec[x-1] = vec[x-1] - 1;
    }
}



int FNegMomOcc(int L, int lmax, int Npar, int last_i, Iarray * ht)
{

/** 'L' is the total angular momentum yet to be subtracted using 'Npar'
    particles in the negative momentum IPS, with their angular momentum
    from -lmax , ..., -1. The problem is addressed on how many ways are
    to an addition of  'Npar'  DIFFERENT positive integer numbers yeild
    the same result  'L' using  'lmax' as maximum term in the addition.
    For each possibility,  assemble a new config.  in the hashing table
    starting from 'last_i' which marks the last index setup previously.
    Return the number of config. setup, to update 'last_i'.

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
        if (Npar > 0)
        {
            return 0;
        }
        else
        {
            for (i = 0; i <= lmax; i++) ht[last_i][i] = 0;
            return 1;
        }
    }

    vec = iarrDef(lmax); // 'vec' of occupations in IPS with l < 0
    Npos = 0; // number of possible config. found

    if (L > lmax) start_l = lmax;
    else          start_l = L;

    // DIFFERENCE FROM THE BOSONIC CASE. The maximum angular momentum
    // is given by a Arithmetic Prog. because of the  Pauli Exclusion
    // principle
    for (l = start_l; (l*(1+l))/2 >= L; l--)
    {
        for (i = 0; i < lmax; i++) vec[i] = 0;
        vec[l-1] = 1;
        // trigger recursion. See the problem of all the additions
        // of positive integers that yield the same result
        FNegMomOcc_rec(L-l,l,vec,lmax,&Npos,Npar,last_i,ht);
    }

    free(vec);
    return Npos;
}



Iarray * FAssembleHT(int N, int lmax, int totalL, int mcsize)
{

/** with the total number of configurations (occupation vectors) given in
    'mcsize',  allocate a matrix  (hashing table)  with 'mcsize' rows and
    2*lmax+1  columns,  whose its elements corresponds to the occupations
    in the IPS, in ascending order of angular momentum, -lmax, ..., lmax.
    Each configuration given in a row has total angular momentum 'totalL'
    and 'N' fermions.
 
    RETURN hashing table as matrix with the occupation along rows **/

    int
        i,
        j,
        k,
        L,
        M,
        N0,
        conf_i,
        prev_i,
        posNC;

    Iarray
        occ;

    Iarray
        * htable;

    // Allocate the hashing table. The number of IPS is 2*lmax+1
    htable = (Iarray *) malloc(mcsize*sizeof(int *));
    for (i = 0; i < mcsize; i++) htable[i] = iarrDef(2*lmax+1);

    M = lmax + 1;       // Number of IPS with l >= 0
    occ = iarrDef(M);   // occupations in IPS with l >= 0

    conf_i = 0;
    for (N0 = 0; N0 <= N; N0++)
    {
        // N0  is the number of particles reserved
        // to use in negative angular momentum IPS
        posNC = NC(N-N0,M);
        for (i = 0; i < posNC; i++)
        {
            IndexToFock(i,N-N0,M,occ);
            L = positiveMomentum(M,occ);
            // Use the remaining particles  'N0' to search for occupations
            // in negative momentum IPS that combined with the occupations
            // in positive momentum  'occ_Lp' give total momentum 'totalL'
            // Note that the maximum contribution of negative momentum is
            // given by an AP starting in -lmax
            if (L - (N0*(2*lmax+1-N0))/2 <= totalL && L >= totalL)
            {
                prev_i = conf_i;
                // assemble negative momentum SPS occupations
                // 'conf' track the configuration index in the hashing table
                conf_i = conf_i + FNegMomOcc(L-totalL,lmax,N0,conf_i,htable);
                // setup the positive angular momentum IPS occupations
                // commom in a stride in the hashing table
                for (k = prev_i; k < conf_i; k++)
                {
                    for (j = lmax; j < 2*lmax+1; j++)
                    {
                        htable[k][j] = occ[j-lmax];
                    }
                }
            }
        }
    }

    free(occ);
    return htable;
}



int FgetIndex(int lmax, int ht_size, Iarray * ht, Iarray occ)
{

/** Search in the hashing table a configuration and  RETURN  its index
    According to assembling procedure the problem  is  split  in  two.
    First the occupations in positive momentum IPS are compared, using
    the two following criteria
        1. Comparing the number of particles in these state with l>= 0
        2. Comparing the occupations in larger momentum IPS
    In the assembling procedure, config. with more particles in positive
    momentum IPS comes before those ones with fewer particles in  l >= 0
    states. Moreover, to distinguish the config. which have the same number
    of particles in l >= 0 states, they are ranked according to the cost
    function, that implies that the occupations are compared starting from
    the largest momentum IPS up to l = 0.
    Second, the occupation in negative momentum IPS are compared starting
    from the most negative momentum IPS l = -lmax, up to l = -1      **/

    int
        i,
        j,
        n,
        lower,
        upper,
        Npos,
        NposConf;

    lower = 0;
    upper = ht_size;

    // First criterion - total number of particles in l >=0 IPS
    NposConf = 0;
    for (j = lmax; j < 2*lmax+1; j++) NposConf = NposConf + occ[j];

    while (1)
    {
        n = (upper + lower) / 2;

        i = 2*lmax;
        while(occ[i] == ht[n][i])
        {
            i = i - 1;
            if (i == lmax) break;
        }

        Npos = 0;
        for (j = lmax; j < 2*lmax+1; j++) Npos = Npos + ht[n][j];
        // Use first criterion
        if (Npos > NposConf)
        {
            lower = n;
            continue;
        }
        else
        {
            if (Npos < NposConf)
            {
                upper = n;
                continue;
            }
            // In case the first criterion does not distinguish
            // look at the occupations in the IPS which the occ
            // are different
            if (i > lmax)
            {
                if (occ[i] > ht[n][i]) lower = n;
                else upper = n;
                continue;
            }
        }

        // Last case, use occupations in the negative momentum IPS
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
