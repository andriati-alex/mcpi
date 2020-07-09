#ifndef _HMatrix_h
#define _HMatrix_h

#include <math.h>
#include "LBoseFockSpace.h"

struct _HConfMat
{

/* Sparse matrix structure for Hamiltonian matrix in Config. space */

    int
        nnze;

    Iarray
        rows,
        cols;

    Carray
        vals;
};

typedef struct _HConfMat * HConfMat;



HConfMat allocEmptyMat(int mcsize, int nnze)
{
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
        filledIPS,
        MemReq;

    Iarray
        z,
        v;

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

        // Add the number of non-zero entries in this row
        nnze = nnze + next_i;
        NNZrow[i] = next_i;

        MemReq = nnze*(sizeof(int)+sizeof(double complex))+mcsize*sizeof(int);
        if (MemReq > MEMORY_TOL)
        {
            printf("\n\nPROCESS ABORTED : The estimated memory required ");
            printf("to set the Hamiltonian matrix will exceed the ");
            printf("tolerance of %.1lf(GB).\n\n",((double) MEMORY_TOL)/1E9);
            exit(EXIT_FAILURE);
        }
    }

    free(v);
    free(z);
    return nnze;
}



HConfMat assembleH(int N, int lmax, int mcsize, Iarray * ht, Carray Ho, double g)
{

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



#endif
