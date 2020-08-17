#ifndef _Hamiltonian_h
#define _Hamiltonian_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "LBoseBoseFockSpace.h"
#include "LBoseFermiFockSpace.h"



HConfMat allocEmptyMat(int mcsize, int nnze)
{

/** Allocate sparse matrix struct with the space dimension 'mcsize'
    and total Number of NonZero Elements 'nnze' in Hamiltonian Mat. **/

    HConfMat
        M;

    M = (struct _HConfMat *) malloc(sizeof(struct _HConfMat));
    M->nnze = nnze;

    M->cols = iarrDef(nnze);
    M->vals = carrDef(nnze);
    M->rows = iarrDef(mcsize+1);
    M->rows[0] = 0;

    return M;
}



void freeHmat(HConfMat M)
{
    free(M->rows);
    free(M->cols);
    free(M->vals);
    free(M);
}



void matmul(int n, Iarray rows, Iarray cols, Carray vals,
            Carray vin, Carray vout)
{

/** Parallelized routine to perform matrix-vector multiplication from
    a sparse matrix structure fields 'rows', 'cols' and 'vals'   **/

    int
        i,
        j,
        threadId,
        nthreads;

    double complex
        z;

#pragma omp parallel private(i,j,z,threadId,nthreads)
    {
        threadId = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        for (i = threadId; i < n; i += nthreads)
        {
            z = 0;
            for (j = rows[i]; j < rows[i+1]; j++)
            {
                z = z + vals[j] * vin[cols[j]];
            }
            vout[i] = z;
        }
    }
}



void InsertNZ_ind(int * next_i, int col, Iarray nz_ind)
{

/** Auxiliar function to record indexes where the Hamiltonian matrix
    has nonzero entries.  'next_i'  marks  the  next position in the
    array 'nz_ind' to set a new column index  IF IT IS NOT CONTAINED
    in the array 'nz_ind' yet. 'col' is the column index where a non
    zero entry was found. **/

    int
        i;

    for (i = 0; i < *next_i; i++)
    {
        // if the column index already appeared it does nothing
        if (nz_ind[i] == col) return;
    }
    // The column index 'col' was not introduced yet.
    nz_ind[*next_i] = col;
    *next_i = *next_i + 1;
}



int NNZ_PerRow(int N, int lmax, int mcsize, Iarray * ht, Iarray NNZrow)
{

/** COMPUTE NUMBER OF NONZERO ENTRIES IN THE HAMILTONIAN USING THE
    MULTICONFIGURATIONAL BASIS
    =================================================================
    To compute the Number of NonZero(NNZ) entries this routine  scans
    all the rules from creation and annihilation op. In this case the
    ones that conserve the total angular momentum.                **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        M,
        nnze,
        next_i,
        filledIPS;

    Iarray
        z,
        v;

    unsigned long
        MemReq;

    M = 2 * lmax + 1;   // total number of IPS
    v = iarrDef(M);     // vector of occupation numbers

    nnze = 0; // total Number of NonZero Elements

    // 'z' stack up all column indexes that contains a non-zero
    // entry in the Hamiltonian matrix.  Of  curse  the maximum
    // size of 'z' is the dimension of the space 'mcsize'.
    z = iarrDef(mcsize);

    for (i = 0; i < mcsize; i++)
    {
        // Identify how many IPS have some occupation
        filledIPS = 0;
        for (k = 0; k < M; k++)
        {
            v[k] = ht[i][k];
            if (v[k] > 0) filledIPS = filledIPS + 1;
        }

        // initialize with the diagonal index  (aways  present)
        // from the rules that remove and replace the particles
        // in the same states,  thus providing the same config.
        z[0] = i;
        next_i = 1;

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < M; q++)
            {
                l = 2 * k - q;
                if (q == k || l < 0 || l > q) continue;

                v[k] -= 2;
                v[l] += 1;
                v[q] += 1;
                j = BgetIndex(lmax,mcsize,ht,v);
                InsertNZ_ind(&next_i,j,z);
                v[k] += 2;
                v[l] -= 1;
                v[q] -= 1;
            }
        }

        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < M; q++)
        {
            for (k = 0; k < M; k++)
            {
                s = 2 * q - k;
                if (q == k || s < 0 || s > k) continue;
                if (v[k] < 1 || v[s] < 1) continue;

                v[k] -= 1;
                v[s] -= 1;
                v[q] += 2;
                j = BgetIndex(lmax,mcsize,ht,v);
                InsertNZ_ind(&next_i,j,z);
                v[k] += 1;
                v[s] += 1;
                v[q] -= 2;
            }
        }

        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = 0; q < M; q++)
                {
                    l = k + s - q;
                    if (q == s || q == k) continue;
                    if (l < 0  || l > q ) continue;

                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 1;
                    v[l] += 1;
                    j = BgetIndex(lmax,mcsize,ht,v);
                    InsertNZ_ind(&next_i,j,z);
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 1;
                    v[l] -= 1;
                }       // Finish q
            }           // Finish s
        }               // Finish k

        if (nnze >= INT_MAX - next_i)
        {
            // With 32-bit integers is impossible to enumerate
            // vector in sparse Hamiltonian matrix
            free(v);
            free(z);
            return -2;
        }
        // Add the number of non-zero entries in this row
        nnze = nnze + next_i;
        NNZrow[i] = next_i;

        MemReq = nnze*(sizeof(int)+sizeof(double complex))+mcsize*sizeof(int);
        if (MemReq > MEMORY_TOL)
        {
            // According to the pre-defined memory tolerance
            // the program will not store the Hamiltonian matrix
            free(v);
            free(z);
            return -1;
        }
    }

    free(v);
    free(z);
    return nnze;
}



int NNZ_PerRow_BBMIX(CompoundSpace S, Iarray NNZrow)
{

/** COMPUTE NUMBER OF NONZERO ENTRIES IN THE HAMILTONIAN USING THE
    MULTICONFIGURATIONAL BASIS FOR 2 SPECIES OF BOSONS
    =================================================================
    To compute the Number of NonZero(NNZ) entries this routine  scans
    all the rules from creation and annihilation op. In this case the
    ones that conserve the total angular momentum.                **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        n,
        ia,
        ib,
        Nla,
        Nlb,
        low,
        up,
        nnze,
        next_i,
        mcsize;

    Iarray
        z,
        vA,
        vB;

    unsigned long
        MemReq;

    Nla = 2 * S->lmaxA + 1; // total number of Individual Particle States A
    Nlb = 2 * S->lmaxB + 1; // total number of Individual Particle States B
    mcsize = S->size;       // dimension of the multiconfig space

    vA = iarrDef(Nla);      // vector of occupation numbers species A
    vB = iarrDef(Nlb);      // vector of occupation numbers species B

    nnze = 0;               // total Number of NonZero Elements

    // 'z' stack up all column indexes that contains a non-zero
    // entry in the Hamiltonian matrix.  Of  curse  the maximum
    // size of 'z' is the dimension of the  space  'mcsize'  in
    // worst case.
    z = iarrDef(mcsize);

    for (i = 0; i < mcsize; i++)
    {
        // copy current configuration occupations
        BBgetConfigs(i,S,vA,vB);

        // Intraspecies interactions neither change the subspace
        // nor influence the individual hashing index.
        n = BBsubIndex(S,vA);
        ia = BgetIndex(S->lmaxA,S->sub[n].sizeA,S->sub[n].hta,vA);
        ib = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,vB);

        // initialize with the diagonal index  (aways  present)
        // from the rules that remove and replace the particles
        // in the same states,  thus providing the same config.
        z[0] = i;
        next_i = 1;

        /*************************************************************
         **                                                         **
         **                                                         **
         **               INTERACTION AMONG A-SPECIES               **
         **                (only off-diagonal rules)                **
         **                                                         **
         **                                                         **
         *************************************************************/

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nla; k++)
        {
            if (vA[k] < 2) continue;
            for (q = 0; q < Nla; q++)
            {
                l = 2 * k - q;
                if (q == k || l < 0 || l > q) continue;

                vA[k] -= 2;
                vA[l] += 1;
                vA[q] += 1;
                // j = BgetIndex(lmax,mcsize,ht,v);
                j = BBgetIndex_A(S,n,ib,vA);
                InsertNZ_ind(&next_i,j,z);
                vA[k] += 2;
                vA[l] -= 1;
                vA[q] -= 1;
            }
        }

        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nla; q++)
        {
            for (k = 0; k < Nla; k++)
            {
                s = 2 * q - k;
                if (q == k || s < 0 || s > k) continue;
                if (vA[k] < 1 || vA[s] < 1) continue;

                vA[k] -= 1;
                vA[s] -= 1;
                vA[q] += 2;
                // j = BgetIndex(lmax,mcsize,ht,v);
                j = BBgetIndex_A(S,n,ib,vA);
                InsertNZ_ind(&next_i,j,z);
                vA[k] += 1;
                vA[s] += 1;
                vA[q] -= 2;
            }
        }

        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nla; k++)
        {
            if (vA[k] < 1) continue;
            for (s = k + 1; s < Nla; s++)
            {
                if (vA[s] < 1) continue;
                for (q = 0; q < Nla; q++)
                {
                    l = k + s - q;
                    if (q == s || q == k) continue;
                    if (l < 0  || l > q ) continue;

                    vA[k] -= 1;
                    vA[s] -= 1;
                    vA[q] += 1;
                    vA[l] += 1;
                    // j = BgetIndex(lmax,mcsize,ht,v);
                    j = BBgetIndex_A(S,n,ib,vA);
                    InsertNZ_ind(&next_i,j,z);
                    vA[k] += 1;
                    vA[s] += 1;
                    vA[q] -= 1;
                    vA[l] -= 1;
                }       // Finish q
            }           // Finish s
        }               // Finish k

        /*************************************************************
         **                                                         **
         **                                                         **
         **               INTERACTION AMONG B-SPECIES               **
         **                (only off-diagonal rules)                **
         **                                                         **
         **                                                         **
         *************************************************************/

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 2) continue;
            for (q = 0; q < Nlb; q++)
            {
                l = 2 * k - q;
                if (q == k || l < 0 || l > q) continue;

                vB[k] -= 2;
                vB[l] += 1;
                vB[q] += 1;
                // j = BgetIndex(lmax,mcsize,ht,v);
                j = BBgetIndex_B(S,n,ia,vB);
                InsertNZ_ind(&next_i,j,z);
                vB[k] += 2;
                vB[l] -= 1;
                vB[q] -= 1;
            }
        }

        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nlb; q++)
        {
            for (k = 0; k < Nlb; k++)
            {
                s = 2 * q - k;
                if (q == k || s < 0 || s > k) continue;
                if (vB[k] < 1 || vB[s] < 1) continue;

                vB[k] -= 1;
                vB[s] -= 1;
                vB[q] += 2;
                // j = BgetIndex(lmax,mcsize,ht,v);
                j = BBgetIndex_B(S,n,ia,vB);
                InsertNZ_ind(&next_i,j,z);
                vB[k] += 1;
                vB[s] += 1;
                vB[q] -= 2;
            }
        }

        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 1) continue;
            for (s = k + 1; s < Nlb; s++)
            {
                if (vB[s] < 1) continue;
                for (q = 0; q < Nlb; q++)
                {
                    l = k + s - q;
                    if (q == s || q == k) continue;
                    if (l < 0  || l > q ) continue;

                    vB[k] -= 1;
                    vB[s] -= 1;
                    vB[q] += 1;
                    vB[l] += 1;
                    // j = BgetIndex(lmax,mcsize,ht,v);
                    j = BBgetIndex_B(S,n,ia,vB);
                    InsertNZ_ind(&next_i,j,z);
                    vB[k] += 1;
                    vB[s] += 1;
                    vB[q] -= 1;
                    vB[l] -= 1;
                }       // Finish q
            }           // Finish s
        }               // Finish k

        /************************************************************
         **                                                        **
         **                                                        **
         **                INTERSPECIES INTERACTION                **
         **                                                        **
         **                                                        **
         ************************************************************/
        // PHYSICAL INTERACTION
        for (k = 0; k < Nla; k++)
        {
            for (s = 0; s < Nlb; s++)
            {
                // NO MOMENTUM EXCHANGED CONTRIBUTES ONLY IN DIAGONAL

                // MOMENTUM EXCHANGED BY q UNITS
                if (s + k < Nlb) low = -k;
                else             low = -(Nlb - 1 - s);
                if (k + s < Nla) up  = s;
                else             up = Nla-1-k;
                for(q = low; q <= up; q++)
                {
                    if (q == 0) continue; // no momentum exchanged
                    if (vA[k+q]*vB[s-q] == 0) continue;
                    vA[k+q] -= 1;
                    vA[k]   += 1;
                    vB[s-q] -= 1;
                    vB[s]   += 1;
                    j = BBgetIndex(S,vA,vB);
                    InsertNZ_ind(&next_i,j,z);
                    vA[k+q] += 1;
                    vA[k]   -= 1;
                    vB[s-q] += 1;
                    vB[s]   -= 1;
                }
            }
        }

        if (nnze >= INT_MAX - next_i)
        {
            // With 32-bit integers is impossible to enumerate
            // vector in sparse Hamiltonian matrix
            free(vA);
            free(vB);
            free(z);
            return -2;
        }
        // Add the number of non-zero entries in this row
        nnze = nnze + next_i;
        NNZrow[i] = next_i;

        MemReq = nnze*(sizeof(int)+sizeof(double complex))+mcsize*sizeof(int);
        if (MemReq > MEMORY_TOL)
        {
            // According to the pre-defined memory tolerance
            // the program will not store the Hamiltonian matrix
            free(vA);
            free(vB);
            free(z);
            return -1;
        }
    }

    free(vA);
    free(vB);
    free(z);
    return nnze;
}



int NNZ_PerRow_BFMIX(BFCompoundSpace S, Iarray NNZrow)
{

/** COMPUTE NUMBER OF NONZERO ENTRIES IN THE HAMILTONIAN USING  THE
    MULTICONFIGURATIONAL BASIS FOR 2 SPECIES OF FERMIONS AND BOSONS
    =================================================================
    To compute the Number of NonZero(NNZ) entries this routine  scans
    all the rules from creation and annihilation op. In this case the
    ones that conserve the total angular momentum.                **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        n,
        iF,
        Nlf,
        Nlb,
        low,
        up,
        nnze,
        next_i,
        mcsize;
    Iarray
        z,
        vB;
    Farray
        vF;
    unsigned long
        MemReq;

    Nlb = 2 * S->lmaxB + 1; // total number of Individual Particle States A
    Nlf = 2 * S->lmaxF + 1; // total number of Individual Particle States B
    mcsize = S->size;       // dimension of the multiconfig space

    vF = farrDef(Nlf);      // vector of occupation numbers Fermions
    vB = iarrDef(Nlb);      // vector of occupation numbers Bosons

    nnze = 0;               // total Number of NonZero Elements

    // 'z' stack up all column indexes that contains a non-zero
    // entry in the Hamiltonian matrix.  Of  curse  the maximum
    // size of 'z' is the dimension of the  space  'mcsize'  in
    // worst case.
    z = iarrDef(mcsize);

    for (i = 0; i < mcsize; i++)
    {
        // copy current configuration occupations
        BFgetConfigs(i,S,vB,vF);

        // Intraspecies interactions neither change the subspace
        // nor influence the individual hashing index.
        n = BFsubIndex(S,vB);
        iF = FgetIndex(S->lmaxF,S->sub[n].sizeF,S->sub[n].htf,vF);

        // initialize with the diagonal index  (aways  present)
        // from the rules that remove and replace the particles
        // in the same states,  thus providing the same config.
        z[0] = i;
        next_i = 1;

        /*************************************************************
         **                                                         **
         **                                                         **
         **               INTERACTION  AMONG BOSONS                 **
         **               (only off-diagonal rules)                 **
         **                                                         **
         **                                                         **
         *************************************************************/

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 2) continue;
            for (q = 0; q < Nlb; q++)
            {
                l = 2 * k - q;
                if (q == k || l < 0 || l > q) continue;

                vB[k] -= 2;
                vB[l] += 1;
                vB[q] += 1;
                // j = BgetIndex(lmax,mcsize,ht,v);
                j = BFgetIndex_B(S,n,iF,vB);
                InsertNZ_ind(&next_i,j,z);
                vB[k] += 2;
                vB[l] -= 1;
                vB[q] -= 1;
            }
        }

        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nlb; q++)
        {
            for (k = 0; k < Nlb; k++)
            {
                s = 2 * q - k;
                if (q == k || s < 0 || s > k) continue;
                if (vB[k] < 1 || vB[s] < 1) continue;

                vB[k] -= 1;
                vB[s] -= 1;
                vB[q] += 2;
                // j = BgetIndex(lmax,mcsize,ht,v);
                j = BFgetIndex_B(S,n,iF,vB);
                InsertNZ_ind(&next_i,j,z);
                vB[k] += 1;
                vB[s] += 1;
                vB[q] -= 2;
            }
        }

        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 1) continue;
            for (s = k + 1; s < Nlb; s++)
            {
                if (vB[s] < 1) continue;
                for (q = 0; q < Nlb; q++)
                {
                    l = k + s - q;
                    if (q == s || q == k) continue;
                    if (l < 0  || l > q ) continue;

                    vB[k] -= 1;
                    vB[s] -= 1;
                    vB[q] += 1;
                    vB[l] += 1;
                    // j = BgetIndex(lmax,mcsize,ht,v);
                    j = BFgetIndex_B(S,n,iF,vB);
                    InsertNZ_ind(&next_i,j,z);
                    vB[k] += 1;
                    vB[s] += 1;
                    vB[q] -= 1;
                    vB[l] -= 1;
                }       // Finish q
            }           // Finish s
        }               // Finish k

        /************************************************************
         **                                                        **
         **                                                        **
         **                INTERSPECIES INTERACTION                **
         **                                                        **
         **                                                        **
         ************************************************************/
        for (s = 0; s < Nlf; s++)
        {
            for (k = 0; k < Nlb; k++)
            {
                // NO MOMENTUM EXCHANGED CONTRIBUTES ONLY IN DIAGONAL

                // MOMENTUM EXCHANGED BY q UNITS
                if (s + k < Nlf) low = -k;
                else             low = -(Nlf - 1 - s);
                if (k + s < Nlb) up  = s;
                else             up = Nlb-1-k;
                for(q = low; q <= up; q++)
                {
                    if (q == 0) continue; // no momentum exchanged
                    if (vB[k+q]*vF[s-q] == 0 || vF[s] > 0) continue;
                    vB[k+q] -= 1;
                    vB[k]   += 1;
                    vF[s-q]  = 0;
                    vF[s]    = 1;
                    j = BFgetIndex(S,vB,vF);
                    InsertNZ_ind(&next_i,j,z);
                    vB[k+q] += 1;
                    vB[k]   -= 1;
                    vF[s-q]  = 1;
                    vF[s]    = 0;
                }
            }
        }

        if (nnze >= INT_MAX - next_i)
        {
            // With 32-bit integers is impossible to enumerate
            // vector in sparse Hamiltonian matrix
            free(vF);
            free(vB);
            free(z);
            return -2;
        }
        // Add the number of non-zero entries in this row
        nnze = nnze + next_i;
        NNZrow[i] = next_i;

        MemReq = nnze*(sizeof(int)+sizeof(double complex))+mcsize*sizeof(int);
        if (MemReq > MEMORY_TOL)
        {
            // According to the pre-defined memory tolerance
            // the program will not store the Hamiltonian matrix
            free(vF);
            free(vB);
            free(z);
            return -1;
        }
    }

    free(vF);
    free(vB);
    free(z);
    return nnze;
}



HConfMat boseH(int N, int lmax, int mcsize, Iarray * ht, Carray Ho, double g)
{

/** SET UP THE HAMILTONIAN MATRIX (SPARSE) FOR SINGLE BOSONIC SPECIES
    IN THE MULTICONFIGURATIONAL BASIS WITH FIXED TOTAL MOMENTUM   **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        p,
        t,
        Nips,
        nnze;

    double
        bosef;

    double complex
        z,
        w;

    Iarray
        v,
        NNZrow;

    HConfMat
        M;

    Nips = 2 * lmax + 1;    // total number of IPS
    v = iarrDef(Nips);      // vector of occupation numbers

    // Compute number of  nonzero  entries to allocate sparse
    // matrix structure including the nonzero entries per row
    NNZrow = iarrDef(mcsize);
    nnze = NNZ_PerRow(N,lmax,mcsize,ht,NNZrow);

    if (nnze < 0)
    {
        // Two problems may have occurred in computing number of nonzero
        // elements of Hamiltonian matrix : 1. 'nnze' exceeded INT_MAX
        // 2. Memory required exceeded the tolerance
        free(v);
        free(NNZrow);
        return NULL;
    }

    M = allocEmptyMat(mcsize,nnze);

    // initialize with arbitrary values
    for (i = 0; i < nnze; i++)
    {
        M->vals[i] = 0;
        M->cols[i] = -1;
    }

    // Configure the strides to find in the sparse storage
    // vectors of values and columns where each row begins
    M->rows[0] = 0;
    for (i = 0; i < mcsize; i++) M->rows[i+1] = M->rows[i] + NNZrow[i];

    // 't' variable tracks what was the last 'index' used
    // in vector of values and columns in sparse matrix
    t = 0;

    for (i = 0; i < mcsize; i++)
    {

        w = 0;
        z = 0;

        // copy current configuration occupations from hashing table
        for (k = 0; k < Nips; k++) v[k] = ht[i][k];

        // THE RULES 0, 1 AND 2 CORRESPOND TO DIAGONAL ELEMENTS OF THE
        // HAMILTONIAN MATRIX, AS SUCH, THE COLUMN IS 'i'.  FOR  THESE
        // CASES NO SEARCH FOR NEW CONF. INDEX IS REQUIRED  SINCE  ALL
        // PARTICLES REMOVED ARE REPLACED IN THE SAME STATES

        // Rule : non-interacting part - creation and annihilation at k
        for (k = 0; k < Nips; k++)
        {
            if (v[k] < 1) continue;
            w = w + Ho[k] * v[k];
        }

        // Rule : Creation on k k / Annihilation on k k
        for (k = 0; k < Nips; k++)
        {
            bosef = v[k] * (v[k] - 1);
            z = z + g * bosef;
        }

        // Rule : Creation on k s / Annihilation on k s
        for (k = 0; k < Nips; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < Nips; s++)
            {
                bosef = v[k] * v[s];
                z = z + 4 * g * bosef;
            }
        }

        // FINISH THE DIAGONAL AND SETUP NON-DIAGONAL ENTRIES
        M->cols[t] = i;
        M->vals[t] = w + z/2;
        t = t + 1;

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nips; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < Nips; q++)
            {
                // Avoid forbidden rules
                if (q == k || 2 * k - q < 0 || 2 * k - q >= Nips) continue;
                l = 2 * k - q;
                if (l > q) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                // assemble the configuration according to action of op.
                v[k] -= 2;
                v[l] += 1;
                v[q] += 1;
                // Compute col index in 'j'
                j = BgetIndex(lmax,mcsize,ht,v);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry q-l
                z = 2 * g * bosef;
                // Check in the current row if the col 'j' has
                // already appeared. Then add its contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + z/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = z/2;
                    M->cols[p] = j;
                    t = t + 1;
                }
                // correct the configuration
                v[k] += 2;
                v[l] -= 1;
                v[q] -= 1;
            }
        }

        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nips; q++)
        {
            for (k = 0; k < Nips; k++)
            {
                // Avoid forbidden rules
                if (q == k || 2 * q - k < 0 || 2 * q - k >= Nips) continue;
                s = 2 * q - k;
                if (v[k] < 1 || v[s] < 1 || s > k) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));
                // assemble the configuration according to action of op.
                v[k] -= 1;
                v[s] -= 1;
                v[q] += 2;
                // Compute col index 'j'
                j = BgetIndex(lmax,mcsize,ht,v);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry k-s
                z = 2 * g * bosef;
                // Check in the current row if the col 'j' has
                // already been initialized.  Add contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + z/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // we get p = t introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = z/2;
                    M->cols[p] = j;
                    t = t + 1;
                    // update the index of the last column initialized
                }
                v[k] += 1;
                v[s] += 1;
                v[q] -= 2;
            }
        }

        // Rule : Creation on k s / Annihilation on s l
        // ONLY IN CASE k = l, BUT this case is included
        // in rule 2

        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nips; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < Nips; s++)
            {
                if (v[s] < 1) continue;
                for (q = 0; q < Nips; q++)
                {
                    if (q == s || q == k) continue;
                    if (k + s - q < 0 || k + s - q >= Nips ) continue;
                    l = k + s - q;
                    if (l <= q) continue;

                    bosef = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 1;
                    v[l] += 1;
                    j = BgetIndex(lmax,mcsize,ht,v);
                    z = 4 * g * bosef;
                    // Check in the current row if the col 'j' has
                    // already been initialized.  Add contribution
                    for (p = M->rows[i]; p < t; p++)
                    {
                        if (j == M->cols[p])
                        {
                            M->vals[p] = M->vals[p] + z/2;
                            break;
                        }
                    }
                    // it the col 'j' has not been initialized yet
                    // introduce the contribution 'z'
                    if (p == t)
                    {
                        M->vals[p] = z/2;
                        M->cols[p] = j;
                        t = t + 1;
                        // update the index of the last column initialized
                    }
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 1;
                    v[l] -= 1;

                }       // Finish q
            }           // Finish s
        }               // Finish k

        if (t != M->rows[i+1])
        {
            printf("\n\nERROR : The number of elements set in the ");
            printf("row %d does not correspond to the stride predicted ",i);
            printf("when computing the number of NonZero Elements(NZE)\n");
            printf("Expected %d but %d NZE were found\n\n",M->rows[i+1],t);
            exit(EXIT_FAILURE);
        }
    }

    free(v);
    free(NNZrow);
    return M;
}



HConfMat boseboseH(CompoundSpace S, Carray HoA, Carray HoB, double g [])
{

/** ASSEMBLE (SPARSE) HAMILTONIAN MATRIX OF 2-SPECIES BOSONIC MIXTURE
    =================================================================
    First compute the number of non-zero elements. If the memory
    required is less than the memory allowed by HMAT_MEM_TOL  it
    return a structure with the sparse matrix                **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        n,
        p,
        t,
        ia,
        ib,
        up,
        low,
        Nla,
        Nlb,
        nnze;
    double
        ga,
        gb,
        gab,
        bosef;
    double complex
        w,
        za,
        zb,
        zab;
    Iarray
        vB,
        vA,
        NNZrow;
    HConfMat
        M;

    Nla = 2 * S->lmaxA + 1; // total number of Individual Particle States A
    Nlb = 2 * S->lmaxB + 1; // total number of Individual Particle States B

    ga = g[0];
    gb = g[1];
    gab = g[2];

    vA = iarrDef(Nla);  // vector of occupation numbers A
    vB = iarrDef(Nlb);  // vector of occupation numbers B

    // Compute number of  nonzero  entries to allocate sparse
    // matrix structure including the nonzero entries per row
    NNZrow = iarrDef(S->size);
    nnze = NNZ_PerRow_BBMIX(S,NNZrow);

    if (nnze < 0)
    {
        // Two problems may have occurred in computing number of nonzero
        // elements of Hamiltonian matrix : 1. 'nnze' exceeded INT_MAX
        // 2. Memory required exceeded the tolerance
        free(vA);
        free(vB);
        free(NNZrow);
        return NULL;
    }

    M = allocEmptyMat(S->size,nnze);

    // initialize with arbitrary values
    for (i = 0; i < nnze; i++)
    {
        M->vals[i] = 0;
        M->cols[i] = -1;
    }

    // Configure the strides to find in the sparse storage
    // vectors of values and columns where each row begins
    M->rows[0] = 0;
    for (i = 0; i < S->size; i++) M->rows[i+1] = M->rows[i] + NNZrow[i];

    // 't' variable tracks what was the last 'index' used
    // in vector of values and columns in sparse matrix
    t = 0;

    for (i = 0; i < S->size; i++)
    {
        w = 0;
        za = 0;
        zb = 0;
        zab = 0;

        // copy current configuration occupations
        BBgetConfigs(i,S,vA,vB);

        // Intraspecies interactions neither change the subspace
        // nor influence the individual hashing index.
        n = BBsubIndex(S,vA);
        ia = BgetIndex(S->lmaxA,S->sub[n].sizeA,S->sub[n].hta,vA);
        ib = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,vB);

        // KINECT ENERGY - creation and annihilation at k
        for (k = 0; k < Nla; k++) w = w + HoA[k]*vA[k];
        for (k = 0; k < Nlb; k++) w = w + HoB[k]*vB[k];

        // Rule : Creation on k k / Annihilation on k k
        for (k = 0; k < Nla; k++)
        {
            bosef = vA[k] * (vA[k] - 1);
            za = za + ga * bosef;
        }
        // Rule : Creation on k s / Annihilation on k s
        for (k = 0; k < Nla; k++)
        {
            if (vA[k] < 1) continue;
            for (s = k + 1; s < Nla; s++)
            {
                bosef = vA[k] * vA[s];
                za = za + 4 * ga * bosef;
            }
        }
        // Rule : Creation on k k / Annihilation on k k
        for (k = 0; k < Nlb; k++)
        {
            bosef = vB[k] * (vB[k] - 1);
            zb = zb + gb * bosef;
        }
        // Rule : Creation on k s / Annihilation on k s
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 1) continue;
            for (s = k + 1; s < Nlb; s++)
            {
                bosef = vB[k] * vB[s];
                zb = zb + 4 * gb * bosef;
            }
        }
        for (k = 0; k < Nla; k++)
        {
            for (s = 0; s < Nlb; s++)
            {
                // NO MOMENTUM EXCHANGED
                bosef = vA[k]*vB[s];
                zab = zab + gab * bosef;
            }
        }

        M->cols[t] = i;
        M->vals[t] = w + za/2 + zb/2 + zab;
        t = t + 1;



        /*************************************************************
         **                                                         **
         **                                                         **
         **               INTERACTION AMONG A-SPECIES               **
         **                (off-diagonal terms only)                **
         **                                                         **
         **                                                         **
         *************************************************************/

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nla; k++)
        {
            if (vA[k] < 2) continue;
            for (q = 0; q < Nla; q++)
            {
                l = 2 * k - q;
                // Avoid repeating/forbidden rules
                if (q == k || l < 0 || l >= q) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)vA[k]*(vA[k]-1)*(vA[q]+1)*(vA[l]+1));
                // assemble the configuration according to action of op.
                vA[k] -= 2;
                vA[l] += 1;
                vA[q] += 1;
                // Compute col index in 'j'
                j = BBgetIndex_A(S,n,ib,vA);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry q-l
                // justifying the choice for only l > q
                za = 2 * ga * bosef;
                // Check in the current row if the col 'j' has
                // already appeared. Then add its contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + za/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = za/2;
                    M->cols[p] = j;
                    t = t + 1;
                }
                // correct the configuration
                vA[k] += 2;
                vA[l] -= 1;
                vA[q] -= 1;
            }
        }
        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nla; q++)
        {
            for (k = 0; k < Nla; k++)
            {
                // Avoid repeating/forbidden rules
                s = 2 * q - k;
                if (q == k || s < 0 || s >= k) continue;
                if (vA[k] < 1 || vA[s] < 1) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)vA[k]*vA[s]*(vA[q]+1)*(vA[q]+2));
                // assemble the configuration according to action of op.
                vA[k] -= 1;
                vA[s] -= 1;
                vA[q] += 2;
                // Compute col index 'j'
                j = BBgetIndex_A(S,n,ib,vA);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry k-s
                za = 2 * ga * bosef;
                // Check in the current row if the col 'j' has
                // already been initialized.  Add contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + za/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // we get p = t introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = za/2;
                    M->cols[p] = j;
                    t = t + 1;
                    // update the index of the last column initialized
                }
                vA[k] += 1;
                vA[s] += 1;
                vA[q] -= 2;
            }
        }
        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nla; k++)
        {
            if (vA[k] < 1) continue;
            for (s = k + 1; s < Nla; s++)
            {
                if (vA[s] < 1) continue;
                for (q = 0; q < Nla; q++)
                {
                    if (q == s || q == k) continue;
                    l = k + s - q;
                    if (l < 0 || l >= q ) continue;

                    bosef = sqrt((double)vA[k]*vA[s]*(vA[q]+1)*(vA[l]+1));
                    vA[k] -= 1;
                    vA[s] -= 1;
                    vA[q] += 1;
                    vA[l] += 1;
                    j = BBgetIndex_A(S,n,ib,vA);
                    za = 4 * ga * bosef;
                    // Check in the current row if the col 'j' has
                    // already been initialized.  Add contribution
                    for (p = M->rows[i]; p < t; p++)
                    {
                        if (j == M->cols[p])
                        {
                            M->vals[p] = M->vals[p] + za/2;
                            break;
                        }
                    }
                    // it the col 'j' has not been initialized yet
                    // introduce the contribution 'z'
                    if (p == t)
                    {
                        M->vals[p] = za/2;
                        M->cols[p] = j;
                        t = t + 1;
                        // update the index of the last column initialized
                    }
                    vA[k] += 1;
                    vA[s] += 1;
                    vA[q] -= 1;
                    vA[l] -= 1;

                }       // Finish q
            }           // Finish s
        }               // Finish k



        /*************************************************************
         **                                                         **
         **                                                         **
         **               INTERACTION AMONG B-SPECIES               **
         **                                                         **
         **                                                         **
         *************************************************************/

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 2) continue;
            for (q = 0; q < Nlb; q++)
            {
                l = 2 * k - q;
                // Avoid repeating/forbidden rules
                if (q == k || l < 0 || l >= q) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)vB[k]*(vB[k]-1)*(vB[q]+1)*(vB[l]+1));
                // assemble the configuration according to action of op.
                vB[k] -= 2;
                vB[l] += 1;
                vB[q] += 1;
                // Compute col index in 'j'
                j = BBgetIndex_B(S,n,ia,vB);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry q-l
                // justifying the choice for only l > q
                zb = 2 * gb * bosef;
                // Check in the current row if the col 'j' has
                // already appeared. Then add its contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + zb/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = zb/2;
                    M->cols[p] = j;
                    t = t + 1;
                }
                // correct the configuration
                vB[k] += 2;
                vB[l] -= 1;
                vB[q] -= 1;
            }
        }
        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nlb; q++)
        {
            for (k = 0; k < Nlb; k++)
            {
                // Avoid repeating/forbidden rules
                s = 2 * q - k;
                if (q == k || s < 0 || s >= k) continue;
                if (vB[k] < 1 || vB[s] < 1) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[q]+2));
                // assemble the configuration according to action of op.
                vB[k] -= 1;
                vB[s] -= 1;
                vB[q] += 2;
                // Compute col index 'j'
                j = BBgetIndex_B(S,n,ia,vB);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry k-s
                zb = 2 * gb * bosef;
                // Check in the current row if the col 'j' has
                // already been initialized.  Add contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + zb/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // we get p = t introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = zb/2;
                    M->cols[p] = j;
                    t = t + 1;
                    // update the index of the last column initialized
                }
                vB[k] += 1;
                vB[s] += 1;
                vB[q] -= 2;
            }
        }
        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 1) continue;
            for (s = k + 1; s < Nlb; s++)
            {
                if (vB[s] < 1) continue;
                for (q = 0; q < Nlb; q++)
                {
                    if (q == s || q == k) continue;
                    l = k + s - q;
                    if (l < 0 || l >= q ) continue;

                    bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[l]+1));
                    vB[k] -= 1;
                    vB[s] -= 1;
                    vB[q] += 1;
                    vB[l] += 1;
                    j = BBgetIndex_B(S,n,ia,vB);
                    zb = 4 * gb * bosef;
                    // Check in the current row if the col 'j' has
                    // already been initialized.  Add contribution
                    for (p = M->rows[i]; p < t; p++)
                    {
                        if (j == M->cols[p])
                        {
                            M->vals[p] = M->vals[p] + zb/2;
                            break;
                        }
                    }
                    // it the col 'j' has not been initialized yet
                    // introduce the contribution 'z'
                    if (p == t)
                    {
                        M->vals[p] = zb/2;
                        M->cols[p] = j;
                        t = t + 1;
                        // update the index of the last column initialized
                    }
                    vB[k] += 1;
                    vB[s] += 1;
                    vB[q] -= 1;
                    vB[l] -= 1;
                }       // Finish q
            }           // Finish s
        }               // Finish k



        /************************************************************
         **                                                        **
         **                                                        **
         **                INTERSPECIES INTERACTION                **
         **                                                        **
         **                                                        **
         ************************************************************/
        for (k = 0; k < Nla; k++)
        {
            for (s = 0; s < Nlb; s++)
            {
                // MOMENTUM EXCHANGED BY q UNITS
                if (s + k < Nlb) low = -k;
                else             low = -(Nlb - 1 - s);
                if (k + s < Nla) up  = s;
                else             up = Nla-1-k;
                for(q = low; q <= up; q++)
                {
                    if (q == 0) continue; // no momentum exchanged
                    if (vA[k+q]*vB[s-q] == 0) continue;
                    bosef = sqrt((double)vA[k+q]*vB[s-q]*(vA[k]+1)*(vB[s]+1));
                    vA[k+q] -= 1;
                    vA[k]   += 1;
                    vB[s-q] -= 1;
                    vB[s]   += 1;
                    j = BBgetIndex(S,vA,vB);
                    zab = gab * bosef;
                    // Check in the current row if the col 'j' has
                    // already been initialized.  Add contribution
                    for (p = M->rows[i]; p < t; p++)
                    {
                        if (j == M->cols[p])
                        {
                            M->vals[p] = M->vals[p] + zab;
                            break;
                        }
                    }
                    // it the col 'j' has not been initialized yet
                    // introduce the contribution 'z'
                    if (p == t)
                    {
                        M->vals[p] = zab;
                        M->cols[p] = j;
                        t = t + 1;
                        // update the index of the last column initialized
                    }
                    vA[k+q] += 1;
                    vA[k]   -= 1;
                    vB[s-q] += 1;
                    vB[s]   -= 1;
                }
            }
        }

        // SANITY CHECK
        if (t != M->rows[i+1])
        {
            printf("\n\nERROR : The number of elements set in the ");
            printf("row %d does not correspond to the stride predicted ",i);
            printf("when computing the number of NonZero Elements(NZE)\n");
            printf("Expected %d but %d NZE were found\n\n",M->rows[i+1],t);
            exit(EXIT_FAILURE);
        }
    }

    free(vA); free(vB); free(NNZrow);
    return M;
}



HConfMat bosefermiH(BFCompoundSpace S, Carray HoB, Carray HoF, double g [])
{

/** ASSEMBLE (SPARSE) HAMILTONIAN MATRIX OF 2-SPECIES BOSE-FERMI MIXTURE
    ====================================================================
    First compute the number of non-zero elements. If the memory
    required is less than the memory allowed by HMAT_MEM_TOL  it
    return a structure with the sparse matrix                **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        n,
        p,
        t,
        iF,
        up,
        low,
        Nlf,
        Nlb,
        nnze,
        fermif;
    double
        gb,
        gbf,
        bosef;
    double complex
        w,
        zb,
        zbf;
    Iarray
        vB,
        NNZrow;
    Farray
        vF;
    HConfMat
        M;

    Nlf = 2 * S->lmaxF + 1; // total number of Individual Particle States A
    Nlb = 2 * S->lmaxB + 1; // total number of Individual Particle States B

    gb = g[0];
    gbf = g[2];

    vF = farrDef(Nlf);  // vector of occupation numbers A
    vB = iarrDef(Nlb);  // vector of occupation numbers B

    // Compute number of  nonzero  entries to allocate sparse
    // matrix structure including the nonzero entries per row
    NNZrow = iarrDef(S->size);
    nnze = NNZ_PerRow_BFMIX(S,NNZrow);

    if (nnze < 0)
    {
        // Two problems may have occurred in computing number of nonzero
        // elements of Hamiltonian matrix : 1. 'nnze' exceeded INT_MAX
        // 2. Memory required exceeded the tolerance
        free(vF);
        free(vB);
        free(NNZrow);
        return NULL;
    }

    M = allocEmptyMat(S->size,nnze);

    // initialize with arbitrary values
    for (i = 0; i < nnze; i++)
    {
        M->vals[i] = 0;
        M->cols[i] = -1;
    }

    // Configure the strides to find in the sparse storage
    // vectors of values and columns where each row begins
    M->rows[0] = 0;
    for (i = 0; i < S->size; i++) M->rows[i+1] = M->rows[i] + NNZrow[i];

    // 't' variable tracks what was the last 'index' used
    // in vector of values and columns in sparse matrix
    t = 0;

    for (i = 0; i < S->size; i++)
    {
        w = 0;
        zb = 0;
        zbf = 0;

        // copy current configuration occupations
        BFgetConfigs(i,S,vB,vF);

        // Intraspecies interactions neither change the subspace
        // nor influence the individual hashing index.
        n = BFsubIndex(S,vB);
        iF = FgetIndex(S->lmaxF,S->sub[n].sizeF,S->sub[n].htf,vF);

        // KINECT ENERGY - creation and annihilation at k
        for (k = 0; k < Nlf; k++) w = w + HoF[k]*vF[k];
        for (k = 0; k < Nlb; k++) w = w + HoB[k]*vB[k];

        // Rule : Creation on k k / Annihilation on k k
        for (k = 0; k < Nlb; k++)
        {
            bosef = vB[k] * (vB[k] - 1);
            zb = zb + gb * bosef;
        }
        // Rule : Creation on k s / Annihilation on k s
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 1) continue;
            for (s = k + 1; s < Nlb; s++)
            {
                bosef = vB[k] * vB[s];
                zb = zb + 4 * gb * bosef;
            }
        }
        // INTERSPECIES INTERACTION WITHOUT MOMENTUM EXCHANGE
        for (k = 0; k < Nlb; k++)
        {
            for (s = 0; s < Nlf; s++)
            {
                bosef = vB[k]*vF[s];
                zbf = zbf + gbf * bosef;
            }
        }

        M->cols[t] = i;
        M->vals[t] = w + zb/2 + zbf;
        t = t + 1;

        /************************************************************
         **                                                        **
         **                                                        **
         **                INTERACTION AMONG BOSONS                **
         **                off-diagonal  terms only                **
         **                                                        **
         **                                                        **
         ************************************************************/

        // Rule : Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 2) continue;
            for (q = 0; q < Nlb; q++)
            {
                l = 2 * k - q;
                // Avoid repeating/forbidden rules
                if (q == k || l < 0 || l >= q) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)vB[k]*(vB[k]-1)*(vB[q]+1)*(vB[l]+1));
                // assemble the configuration according to action of op.
                vB[k] -= 2;
                vB[l] += 1;
                vB[q] += 1;
                // Compute col index in 'j'
                j = BFgetIndex_B(S,n,iF,vB);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry q-l
                // justifying the choice for only l > q
                zb = 2 * gb * bosef;
                // Check in the current row if the col 'j' has
                // already appeared. Then add its contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + zb/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = zb/2;
                    M->cols[p] = j;
                    t = t + 1;
                }
                // correct the configuration
                vB[k] += 2;
                vB[l] -= 1;
                vB[q] -= 1;
            }
        }
        // Rule : Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < Nlb; q++)
        {
            for (k = 0; k < Nlb; k++)
            {
                // Avoid repeating/forbidden rules
                s = 2 * q - k;
                if (q == k || s < 0 || s >= k) continue;
                if (vB[k] < 1 || vB[s] < 1) continue;

                // compute bosonic factor from op. action
                bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[q]+2));
                // assemble the configuration according to action of op.
                vB[k] -= 1;
                vB[s] -= 1;
                vB[q] += 2;
                // Compute col index 'j'
                j = BFgetIndex_B(S,n,iF,vB);
                // value to be set in the matrix at col 'j'
                // factor 2 counts for the symmetry k-s
                zb = 2 * gb * bosef;
                // Check in the current row if the col 'j' has
                // already been initialized.  Add contribution
                for (p = M->rows[i]; p < t; p++)
                {
                    if (j == M->cols[p])
                    {
                        M->vals[p] = M->vals[p] + zb/2;
                        break;
                    }
                }
                // it the col 'j' has not been initialized yet
                // we get p = t introduce the contribution 'z'
                if (p == t)
                {
                    M->vals[p] = zb/2;
                    M->cols[p] = j;
                    t = t + 1;
                    // update the index of the last column initialized
                }
                vB[k] += 1;
                vB[s] += 1;
                vB[q] -= 2;
            }
        }
        // Rule : Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < Nlb; k++)
        {
            if (vB[k] < 1) continue;
            for (s = k + 1; s < Nlb; s++)
            {
                if (vB[s] < 1) continue;
                for (q = 0; q < Nlb; q++)
                {
                    if (q == s || q == k) continue;
                    l = k + s - q;
                    if (l < 0 || l >= q ) continue;

                    bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[l]+1));
                    vB[k] -= 1;
                    vB[s] -= 1;
                    vB[q] += 1;
                    vB[l] += 1;
                    j = BFgetIndex_B(S,n,iF,vB);
                    zb = 4 * gb * bosef;
                    // Check in the current row if the col 'j' has
                    // already been initialized.  Add contribution
                    for (p = M->rows[i]; p < t; p++)
                    {
                        if (j == M->cols[p])
                        {
                            M->vals[p] = M->vals[p] + zb/2;
                            break;
                        }
                    }
                    // it the col 'j' has not been initialized yet
                    // introduce the contribution 'z'
                    if (p == t)
                    {
                        M->vals[p] = zb/2;
                        M->cols[p] = j;
                        t = t + 1;
                        // update the index of the last column initialized
                    }
                    vB[k] += 1;
                    vB[s] += 1;
                    vB[q] -= 1;
                    vB[l] -= 1;
                }       // Finish q
            }           // Finish s
        }               // Finish k

        /************************************************************
         **                                                        **
         **                                                        **
         **                INTERSPECIES INTERACTION                **
         **                                                        **
         **                                                        **
         ************************************************************/
        for (s = 0; s < Nlf; s++)
        {
            for (k = 0; k < Nlb; k++)
            {
                // MOMENTUM EXCHANGED BY q UNITS
                if (s + k < Nlf) low = -k;
                else             low = -(Nlf - 1 - s);
                if (k + s < Nlb) up = s;
                else             up = Nlb-1-k;
                for(q = low; q <= up; q++)
                {
                    // Avoid case of without momentum exchanged
                    if (q == 0) continue;
                    // avoid unpopulated states and creation on
                    // already occupied fermionic state
                    if (vB[k+q]*vF[s-q] == 0 || vF[s] > 0) continue;
                    // result of bosonic operators action
                    bosef  = sqrt((double)vB[k+q]*(vB[k]+1));
                    vB[k+q] -= 1;
                    vB[k]   += 1;
                    // fermi factor is taken into account moving the
                    // operators and inverting sign at each crossing
                    // among them until find the right place
                    fermif = 1;
                    for (l = 0; l < s-q; l++)
                    {
                        if (vF[l] == 1) fermif = (-1)*fermif;
                    }
                    vF[s-q] = 0;
                    for (l = 0; l < s; l++)
                    {
                        if (vF[l] == 1) fermif = (-1)*fermif;
                    }
                    vF[s] = 1;
                    j = BFgetIndex(S,vB,vF);
                    zbf = gbf * bosef * fermif;
                    // Check in the current row if the col 'j' has
                    // already been initialized.  Add contribution
                    for (p = M->rows[i]; p < t; p++)
                    {
                        if (j == M->cols[p])
                        {
                            M->vals[p] = M->vals[p] + zbf;
                            break;
                        }
                    }
                    // it the col 'j' has not been initialized yet
                    // introduce the contribution 'z'
                    if (p == t)
                    {
                        M->vals[p] = zbf;
                        M->cols[p] = j;
                        t = t + 1;
                        // update the index of the last column initialized
                    }
                    vB[k+q] += 1;
                    vB[k]   -= 1;
                    vF[s-q]  = 1;
                    vF[s]    = 0;
                }
            }
        }

        // SANITY CHECK
        if (t != M->rows[i+1])
        {
            printf("\n\nERROR : The number of elements set in the ");
            printf("row %d does not correspond to the stride predicted ",i);
            printf("when computing the number of NonZero Elements(NZE)\n");
            printf("Expected %d but %d NZE were found\n\n",M->rows[i+1],t);
            exit(EXIT_FAILURE);
        }
    }

    free(vF); free(vB); free(NNZrow);
    return M;
}



void actH(int lmax, int mcsize, Iarray * ht, Carray Ho, double g,
          Carray Cin, Carray Cout)
{

/** ROUTINE TO APPLY HAMILTONIAN IN A STATE EXPRESSED IN CONFIG. BASIS
    WITH MOMENTUM CONSERVATION AND CONTACT INTERACTION IN 1D
    ------------------------------------------------------------------
    Given the one-body matrices elements(diagonal only) in  individual
    particle momentum basis 'Ho' and the contact interacting  strength
    parameter 'g' act with the Hamiltonian in a state expressed by its
    coefficients 'Cin' and record the result in 'Cout'.  The  size  of
    the coefficients is the dimension of the multiconfig. space mcsize
    and the enumeration for the config. given in the hashing table 'ht' **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        Nips;

    unsigned int
        threadId,
        nthreads;

    double
        bosef;

    double complex
        w,
        z;

    Iarray
        v;

    Nips = 2 * lmax + 1;    // total number of Individual Particle States

    #pragma omp parallel private(i,j,k,s,q,l,bosef,z,w,v,threadId,nthreads)
    {
        threadId = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        v = iarrDef(Nips);      // vector of occupation numbers

        for (i = threadId; i < mcsize; i += nthreads)
        {
            w = 0;
            z = 0;

            // copy current configuration occupations from hashing table
            for (k = 0; k < Nips; k++) v[k] = ht[i][k];

            // THE THREE FIRST RULES CORRESPOND TO DIAGONAL ELEMENTS OF THE
            // HAMILTONIAN MATRIX, AS SUCH,  THE COLUMN IS 'i'.  FOR  THESE
            // CASES NO SEARCH FOR NEW CONFIG. INDEX IS REQUIRED, SINCE ALL
            // PARTICLES REMOVED ARE REPLACED IN THE SAME STATES

            // Rule : non-interacting part - creation and annihilation at k
            for (k = 0; k < Nips; k++)
            {
                if (v[k] < 1) continue;
                w = w + Ho[k] * v[k] * Cin[i];
            }

            // Rule : Creation on k k / Annihilation on k k
            for (k = 0; k < Nips; k++)
            {
                bosef = v[k] * (v[k] - 1);
                z = z + g * bosef * Cin[i];
            }

            // Rule : Creation on k s / Annihilation on k s
            for (k = 0; k < Nips; k++)
            {
                if (v[k] < 1) continue;
                for (s = k + 1; s < Nips; s++)
                {
                    bosef = v[k] * v[s];
                    z = z + 4 * g * bosef * Cin[i];
                }
            }

            // FINISH THE DIAGONAL AND SETUP NON-DIAGONAL ENTRIES

            // Rule : Creation on k k / Annihilation on q l
            // ONLY IN CASE q + l = 2 * k
            for (k = 0; k < Nips; k++)
            {
                if (v[k] < 2) continue;
                for (q = 0; q < Nips; q++)
                {
                    l = 2 * k - q;
                    // Avoid repeating/forbidden rules
                    if (q == k || l < 0 || l >= q) continue;

                    // compute bosonic factor from op. action
                    bosef = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                    // assemble the configuration according to action of op.
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    // Compute col index in 'j'
                    j = BgetIndex(lmax,mcsize,ht,v);
                    // factor 2 counts for the symmetry q-l
                    // justifying the choice for only l > q
                    z = z + 2 * g * bosef * Cin[j];
                    // correct the configuration
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }

            // Rule : Creation on k s / Annihilation on q q
            // ONLY IN CASE k + s = 2 * q
            for (q = 0; q < Nips; q++)
            {
                for (k = 0; k < Nips; k++)
                {
                    // Avoid repeating/forbidden rules
                    s = 2 * q - k;
                    if (q == k || s < 0 || s >= k) continue;
                    if (v[k] < 1 || v[s] < 1) continue;

                    // compute bosonic factor from op. action
                    bosef = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));
                    // assemble the configuration according to action of op.
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    // Compute col index 'j'
                    j = BgetIndex(lmax,mcsize,ht,v);
                    // factor 2 counts for the symmetry k-s
                    z = z + 2 * g * bosef * Cin[j];
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }

            // Rule : Creation on k s / Annihilation on s l
            // ONLY IN CASE k = l, BUT this case is included
            // in rule 2

            // Rule : Creation on k s / Annihilation on q l
            // ONLY IN CASE k + s = q + l
            for (k = 0; k < Nips; k++)
            {
                if (v[k] < 1) continue;
                for (s = k + 1; s < Nips; s++)
                {
                    if (v[s] < 1) continue;
                    for (q = 0; q < Nips; q++)
                    {
                        if (q == s || q == k) continue;
                        l = k + s - q;
                        if (l < 0 || l >= q ) continue;

                        bosef = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = BgetIndex(lmax,mcsize,ht,v);
                        z = z + 4 * g * bosef * Cin[j];
                        v[k] += 1;
                        v[s] += 1;
                        v[q] -= 1;
                        v[l] -= 1;

                    }       // Finish q
                }           // Finish s
            }               // Finish k

            Cout[i] = w + z/2;
        }

        free(v);
    }
}



void mixture_actH(CompoundSpace S, Carray HoA, Carray HoB, double g [],
                  Carray Cin, Carray Cout)
{

/** ROUTINE TO APPLY HAMILTONIAN IN A STATE EXPRESSED IN CONFIG. BASIS
    WITH MOMENTUM CONSERVATION AND CONTACT INTERACTION  IN  1D  FOR  A
    MIXTURE OF TWO BOSONIC SPECIES
    ------------------------------------------------------------------
    Given the one-body matrices elements (diagonal only)  in  single
    particle momentum basis 'HoA' and 'HoB', the contact interacting
    strength parameters 'g' sorted as gaa, gbb and gab, act with the
    Hamiltonian in a state expressed by its coefficients  'Cin'  and
    record the result in 'Cout'.  'S' is the structure corresponding
    to the basis of the two component system.                    **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        n,
        ia,
        ib,
        up,
        low,
        Nla,
        Nlb;

    unsigned int
        threadId,
        nthreads;

    double
        ga,
        gb,
        gab,
        bosef;

    double complex
        w,
        za,
        zb,
        zab;

    Iarray
        vB,
        vA;

    Nla = 2 * S->lmaxA + 1; // total number of Individual Particle States A
    Nlb = 2 * S->lmaxB + 1; // total number of Individual Particle States B

    ga = g[0];
    gb = g[1];
    gab = g[2];

    #pragma omp parallel \
    private(i,j,k,s,q,l,n,ia,ib,up,low,bosef,w,za,zb,zab,vA,vB,threadId,nthreads)
    {
        threadId = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        vA = iarrDef(Nla);  // vector of occupation numbers
        vB = iarrDef(Nlb);

        for (i = threadId; i < S->size; i += nthreads)
        {
            w = 0;
            za = 0;
            zb = 0;
            zab = 0;

            // copy current configuration occupations
            BBgetConfigs(i,S,vA,vB);

            // Intraspecies interactions neither change the subspace
            // nor influence the individual hashing index.
            n = BBsubIndex(S,vA);
            ia = BgetIndex(S->lmaxA,S->sub[n].sizeA,S->sub[n].hta,vA);
            ib = BgetIndex(S->lmaxB,S->sub[n].sizeB,S->sub[n].htb,vB);

            // KINECT ENERGY - creation and annihilation at k
            for (k = 0; k < Nla; k++) w = w + HoA[k]*vA[k]*Cin[i];
            for (k = 0; k < Nlb; k++) w = w + HoB[k]*vB[k]*Cin[i];



            /*************************************************************
             **                                                         **
             **                                                         **
             **               INTERACTION AMONG A-SPECIES               **
             **                                                         **
             **                                                         **
             *************************************************************/

            // Rule : Creation on k k / Annihilation on k k
            for (k = 0; k < Nla; k++)
            {
                bosef = vA[k] * (vA[k] - 1);
                za = za + ga * bosef * Cin[i];
            }
            // Rule : Creation on k s / Annihilation on k s
            for (k = 0; k < Nla; k++)
            {
                if (vA[k] < 1) continue;
                for (s = k + 1; s < Nla; s++)
                {
                    bosef = vA[k] * vA[s];
                    za = za + 4 * ga * bosef * Cin[i];
                }
            }
            // Rule : Creation on k k / Annihilation on q l
            // ONLY IN CASE q + l = 2 * k
            for (k = 0; k < Nla; k++)
            {
                if (vA[k] < 2) continue;
                for (q = 0; q < Nla; q++)
                {
                    l = 2 * k - q;
                    // Avoid repeating/forbidden rules
                    if (q == k || l < 0 || l >= q) continue;

                    // compute bosonic factor from op. action
                    bosef = sqrt((double)vA[k]*(vA[k]-1)*(vA[q]+1)*(vA[l]+1));
                    // assemble the configuration according to action of op.
                    vA[k] -= 2;
                    vA[l] += 1;
                    vA[q] += 1;
                    // Compute col index in 'j'
                    // j = BBgetIndex(S,vA,vB);
                    j = BBgetIndex_A(S,n,ib,vA);
                    // factor 2 counts for the symmetry q-l
                    // justifying the choice for only l > q
                    za = za + 2 * ga * bosef * Cin[j];
                    // correct the configuration
                    vA[k] += 2;
                    vA[l] -= 1;
                    vA[q] -= 1;
                }
            }
            // Rule : Creation on k s / Annihilation on q q
            // ONLY IN CASE k + s = 2 * q
            for (q = 0; q < Nla; q++)
            {
                for (k = 0; k < Nla; k++)
                {
                    // Avoid repeating/forbidden rules
                    s = 2 * q - k;
                    if (q == k || s < 0 || s >= k) continue;
                    if (vA[k] < 1 || vA[s] < 1) continue;

                    // compute bosonic factor from op. action
                    bosef = sqrt((double)vA[k]*vA[s]*(vA[q]+1)*(vA[q]+2));
                    // assemble the configuration according to action of op.
                    vA[k] -= 1;
                    vA[s] -= 1;
                    vA[q] += 2;
                    // Compute col index 'j'
                    // j = BBgetIndex(S,vA,vB);
                    j = BBgetIndex_A(S,n,ib,vA);
                    // factor 2 counts for the symmetry k-s
                    za = za + 2 * ga * bosef * Cin[j];
                    vA[k] += 1;
                    vA[s] += 1;
                    vA[q] -= 2;
                }
            }
            // Rule : Creation on k s / Annihilation on q l
            // ONLY IN CASE k + s = q + l
            for (k = 0; k < Nla; k++)
            {
                if (vA[k] < 1) continue;
                for (s = k + 1; s < Nla; s++)
                {
                    if (vA[s] < 1) continue;
                    for (q = 0; q < Nla; q++)
                    {
                        if (q == s || q == k) continue;
                        l = k + s - q;
                        if (l < 0 || l >= q ) continue;

                        bosef = sqrt((double)vA[k]*vA[s]*(vA[q]+1)*(vA[l]+1));
                        vA[k] -= 1;
                        vA[s] -= 1;
                        vA[q] += 1;
                        vA[l] += 1;
                        // j = BBgetIndex(S,vA,vB);
                        j = BBgetIndex_A(S,n,ib,vA);
                        za = za + 4 * ga * bosef * Cin[j];
                        vA[k] += 1;
                        vA[s] += 1;
                        vA[q] -= 1;
                        vA[l] -= 1;

                    }       // Finish q
                }           // Finish s
            }               // Finish k



            /*************************************************************
             **                                                         **
             **                                                         **
             **               INTERACTION AMONG B-SPECIES               **
             **                                                         **
             **                                                         **
             *************************************************************/

            // Rule : Creation on k k / Annihilation on k k
            for (k = 0; k < Nlb; k++)
            {
                bosef = vB[k] * (vB[k] - 1);
                zb = zb + gb * bosef * Cin[i];
            }
            // Rule : Creation on k s / Annihilation on k s
            for (k = 0; k < Nlb; k++)
            {
                if (vB[k] < 1) continue;
                for (s = k + 1; s < Nlb; s++)
                {
                    bosef = vB[k] * vB[s];
                    zb = zb + 4 * gb * bosef * Cin[i];
                }
            }
            // Rule : Creation on k k / Annihilation on q l
            // ONLY IN CASE q + l = 2 * k
            for (k = 0; k < Nlb; k++)
            {
                if (vB[k] < 2) continue;
                for (q = 0; q < Nlb; q++)
                {
                    l = 2 * k - q;
                    // Avoid repeating/forbidden rules
                    if (q == k || l < 0 || l >= q) continue;

                    // compute bosonic factor from op. action
                    bosef = sqrt((double)vB[k]*(vB[k]-1)*(vB[q]+1)*(vB[l]+1));
                    // assemble the configuration according to action of op.
                    vB[k] -= 2;
                    vB[l] += 1;
                    vB[q] += 1;
                    // Compute col index in 'j'
                    // j = BBgetIndex(S,vA,vB);
                    j = BBgetIndex_B(S,n,ia,vB);
                    // factor 2 counts for the symmetry q-l
                    // justifying the choice for only l > q
                    zb = zb + 2 * gb * bosef * Cin[j];
                    // correct the configuration
                    vB[k] += 2;
                    vB[l] -= 1;
                    vB[q] -= 1;
                }
            }
            // Rule : Creation on k s / Annihilation on q q
            // ONLY IN CASE k + s = 2 * q
            for (q = 0; q < Nlb; q++)
            {
                for (k = 0; k < Nlb; k++)
                {
                    // Avoid repeating/forbidden rules
                    s = 2 * q - k;
                    if (q == k || s < 0 || s >= k) continue;
                    if (vB[k] < 1 || vB[s] < 1) continue;

                    // compute bosonic factor from op. action
                    bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[q]+2));
                    // assemble the configuration according to action of op.
                    vB[k] -= 1;
                    vB[s] -= 1;
                    vB[q] += 2;
                    // Compute col index 'j'
                    // j = BBgetIndex(S,vA,vB);
                    j = BBgetIndex_B(S,n,ia,vB);
                    // factor 2 counts for the symmetry k-s
                    zb = zb + 2 * gb * bosef * Cin[j];
                    vB[k] += 1;
                    vB[s] += 1;
                    vB[q] -= 2;
                }
            }
            // Rule : Creation on k s / Annihilation on q l
            // ONLY IN CASE k + s = q + l
            for (k = 0; k < Nlb; k++)
            {
                if (vB[k] < 1) continue;
                for (s = k + 1; s < Nlb; s++)
                {
                    if (vB[s] < 1) continue;
                    for (q = 0; q < Nlb; q++)
                    {
                        if (q == s || q == k) continue;
                        l = k + s - q;
                        if (l < 0 || l >= q ) continue;

                        bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[l]+1));
                        vB[k] -= 1;
                        vB[s] -= 1;
                        vB[q] += 1;
                        vB[l] += 1;
                        // j = BBgetIndex(S,vA,vB);
                        j = BBgetIndex_B(S,n,ia,vB);
                        zb = zb + 4 * gb * bosef * Cin[j];
                        vB[k] += 1;
                        vB[s] += 1;
                        vB[q] -= 1;
                        vB[l] -= 1;

                    }       // Finish q
                }           // Finish s
            }               // Finish k



            /************************************************************
             **                                                        **
             **                                                        **
             **                INTERSPECIES INTERACTION                **
             **                                                        **
             **                                                        **
             ************************************************************/
            // PHYSICAL INTERACTION
            for (k = 0; k < Nla; k++)
            {
                for (s = 0; s < Nlb; s++)
                {
                    // NO MOMENTUM EXCHANGED
                    bosef = vA[k]*vB[s];
                    zab = zab + gab * bosef * Cin[i];

                    // MOMENTUM EXCHANGED BY q UNITS
                    if (s + k < Nlb) low = -k;
                    else             low = -(Nlb - 1 - s);
                    if (k + s < Nla) up  = s;
                    else             up = Nla-1-k;
                    for(q = low; q <= up; q++)
                    {
                        if (q == 0) continue; // no momentum exchanged
                        if (vA[k+q]*vB[s-q] == 0) continue;
                        bosef = sqrt((double)vA[k+q]*vB[s-q]*\
                                     (vA[k]+1)*(vB[s]+1));
                        vA[k+q] -= 1;
                        vA[k]   += 1;
                        vB[s-q] -= 1;
                        vB[s]   += 1;
                        j = BBgetIndex(S,vA,vB);
                        zab = zab + gab * bosef * Cin[j];
                        vA[k+q] += 1;
                        vA[k]   -= 1;
                        vB[s-q] += 1;
                        vB[s]   -= 1;
                    }
                }
            }

            Cout[i] = w + za/2 + zb/2 + zab;

        }   // FINISH LOOP OVER SPACE DIMENSION

        free(vA); free(vB);

    }   // FINISH PARALLEL REGION
}



void bosefermi_actH(BFCompoundSpace S, Carray HoB, Carray HoF, double g [],
                    Carray Cin, Carray Cout)
{

/** ROUTINE TO APPLY HAMILTONIAN IN A STATE EXPRESSED IN CONFIG. BASIS
    WITH MOMENTUM CONSERVATION AND CONTACT INTERACTION  IN  1D  FOR  A
    MIXTURE OF TWO BOSONS AND FERMIONS
    ------------------------------------------------------------------
    Given the one-body matrices elements (diagonal only)  in  single
    particle momentum basis 'HoB' and 'HoF', the contact interacting
    strength parameters 'g' sorted as gbb, gff and gbf, act with the
    Hamiltonian in a state expressed by its coefficients  'Cin'  and
    record the result in 'Cout'.  'S' is the structure corresponding
    to the basis of the two component system.                    **/

    int
        i,
        j,
        k,
        l,
        s,
        q,
        n,
        iB,
        iF,
        up,
        low,
        Nlb,
        Nlf,
        fermif;
    unsigned int
        threadId,
        nthreads;
    double
        gb,
        gbf,
        bosef;
    double complex
        w,
        zb,
        zbf;
    Iarray
        vB;
    Farray
        vF;

    Nlb = 2 * S->lmaxB + 1; // total number of single particle states Bosons
    Nlf = 2 * S->lmaxF + 1; // total number of single Particle States Fermions

    // There is no Fermi-Fermi contact interaction
    gb = g[0];
    gbf = g[2];

    #pragma omp parallel private(i,j,k,s,q,l,n,iB,iF,up,low,bosef,fermif,\
                                 w,zb,zbf,vB,vF,threadId,nthreads)
    {
        threadId = omp_get_thread_num();
        nthreads = omp_get_num_threads();

        vB = iarrDef(Nlb);  // vector of occupation numbers
        vF = farrDef(Nlf);

        for (i = threadId; i < S->size; i += nthreads)
        {
            w = 0;
            zb = 0;
            zbf = 0;

            // copy current configuration occupations
            BFgetConfigs(i,S,vB,vF);

            // Intraspecies interactions neither change the subspace
            // nor influence the individual hashing index.
            n = BFsubIndex(S,vB);
            iF = FgetIndex(S->lmaxF,S->sub[n].sizeF,S->sub[n].htf,vF);

            // KINECT ENERGY - creation and annihilation at k
            for (k = 0; k < Nlb; k++) w = w + HoB[k]*vB[k]*Cin[i];
            for (k = 0; k < Nlf; k++) w = w + HoF[k]*vF[k]*Cin[i];

            /************************************************************
             **                                                        **
             **                                                        **
             **                INTERACTION AMONG BOSONS                **
             **                                                        **
             **                                                        **
             ************************************************************/

            // Rule : Creation on k k / Annihilation on k k
            for (k = 0; k < Nlb; k++)
            {
                bosef = vB[k] * (vB[k] - 1);
                zb = zb + gb * bosef * Cin[i];
            }
            // Rule : Creation on k s / Annihilation on k s
            for (k = 0; k < Nlb; k++)
            {
                if (vB[k] < 1) continue;
                for (s = k + 1; s < Nlb; s++)
                {
                    bosef = vB[k] * vB[s];
                    zb = zb + 4 * gb * bosef * Cin[i];
                }
            }
            // Rule : Creation on k k / Annihilation on q l
            // ONLY IN CASE q + l = 2 * k
            for (k = 0; k < Nlb; k++)
            {
                if (vB[k] < 2) continue;
                for (q = 0; q < Nlb; q++)
                {
                    l = 2 * k - q;
                    // Avoid repeating/forbidden rules
                    if (q == k || l < 0 || l >= q) continue;
                    // compute bosonic factor from op. action
                    bosef = sqrt((double)vB[k]*(vB[k]-1)*(vB[q]+1)*(vB[l]+1));
                    // assemble the configuration according to action of op.
                    vB[k] -= 2;
                    vB[l] += 1;
                    vB[q] += 1;
                    // Compute col index in 'j'
                    // j = BBgetIndex(S,vA,vB);
                    j = BFgetIndex_B(S,n,iF,vB);
                    // factor 2 counts for the symmetry q-l
                    // justifying the choice for only l > q
                    zb = zb + 2 * gb * bosef * Cin[j];
                    // correct the configuration
                    vB[k] += 2;
                    vB[l] -= 1;
                    vB[q] -= 1;
                }
            }
            // Rule : Creation on k s / Annihilation on q q
            // ONLY IN CASE k + s = 2 * q
            for (q = 0; q < Nlb; q++)
            {
                for (k = 0; k < Nlb; k++)
                {
                    // Avoid repeating/forbidden rules
                    s = 2 * q - k;
                    if (q == k || s < 0 || s >= k) continue;
                    if (vB[k] < 1 || vB[s] < 1) continue;
                    // compute bosonic factor from op. action
                    bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[q]+2));
                    // assemble the configuration according to action of op.
                    vB[k] -= 1;
                    vB[s] -= 1;
                    vB[q] += 2;
                    // Compute col index 'j'
                    // j = BBgetIndex(S,vA,vB);
                    j = BFgetIndex_B(S,n,iF,vB);
                    // factor 2 counts for the symmetry k-s
                    zb = zb + 2 * gb * bosef * Cin[j];
                    vB[k] += 1;
                    vB[s] += 1;
                    vB[q] -= 2;
                }
            }
            // Rule : Creation on k s / Annihilation on q l
            // ONLY IN CASE k + s = q + l
            for (k = 0; k < Nlb; k++)
            {
                if (vB[k] < 1) continue;
                for (s = k + 1; s < Nlb; s++)
                {
                    if (vB[s] < 1) continue;
                    for (q = 0; q < Nlb; q++)
                    {
                        if (q == s || q == k) continue;
                        l = k + s - q;
                        if (l < 0 || l >= q ) continue;

                        bosef = sqrt((double)vB[k]*vB[s]*(vB[q]+1)*(vB[l]+1));
                        vB[k] -= 1;
                        vB[s] -= 1;
                        vB[q] += 1;
                        vB[l] += 1;
                        // j = BBgetIndex(S,vA,vB);
                        j = BFgetIndex_B(S,n,iF,vB);
                        zb = zb + 4 * gb * bosef * Cin[j];
                        vB[k] += 1;
                        vB[s] += 1;
                        vB[q] -= 1;
                        vB[l] -= 1;

                    }       // Finish q
                }           // Finish s
            }               // Finish k

            /************************************************************
             **                                                        **
             **                                                        **
             **                INTERSPECIES INTERACTION                **
             **                                                        **
             **                                                        **
             ************************************************************/
            for (s = 0; s < Nlf; s++)
            {
                for (k = 0; k < Nlb; k++)
                {
                    // NO MOMENTUM EXCHANGED CONTRIBUTION
                    bosef = vB[k]*vF[s];
                    zbf = zbf + gbf * bosef * Cin[i];

                    // MOMENTUM EXCHANGED BY q UNITS STILL INSIDE
                    // THE MULTICONFIGURATIONAL SPACE DEFINED
                    if (s + k < Nlf) low = -k;
                    else             low = -(Nlf - 1 - s);
                    if (k + s < Nlb) up = s;
                    else             up = Nlb-1-k;
                    for(q = low; q <= up; q++)
                    {
                        // Avoid case of without momentum exchanged
                        if (q == 0) continue;
                        // avoid unpopulated states and creation on
                        // already occupied fermionic state
                        if (vB[k+q]*vF[s-q] == 0 || vF[s] > 0) continue;
                        // result of bosonic operators action
                        bosef  = sqrt((double)vB[k+q]*(vB[k]+1));
                        vB[k+q] -= 1;
                        vB[k]   += 1;
                        // fermi factor is taken into account moving the
                        // operators and inverting sign at each crossing
                        // among them until find the right place
                        fermif = 1;
                        for (l = 0; l < s-q; l++)
                        {
                            if (vF[l] == 1) fermif = (-1)*fermif;
                        }
                        vF[s-q] = 0;
                        for (l = 0; l < s; l++)
                        {
                            if (vF[l] == 1) fermif = (-1)*fermif;
                        }
                        vF[s] = 1;
                        j = BFgetIndex(S,vB,vF);
                        zbf = zbf + gbf * bosef * fermif * Cin[j];
                        vB[k+q] += 1;
                        vB[k]   -= 1;
                        vF[s-q]  = 1;
                        vF[s]    = 0;
                    }
                }
            }

            Cout[i] = w + zb/2 + zbf;

        }   // FINISH LOOP OVER SPACE DIMENSION

        free(vF); free(vB);

    }   // FINISH PARALLEL REGION
}



#endif
