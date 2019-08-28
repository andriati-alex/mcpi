#ifndef _hamiltonianMatrix_h
#define _hamiltonianMatrix_h

#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "configurationsMap.h"



void applyHconf (int N, int M, Iarray NCmat, Iarray IF, Carray C, Cmatrix Ho,
     Carray Hint, Carray out)
{
    // Apply the many-body hamiltonian in a state expressed in
    // number-occupation basis with coefficients defined by C.



    int // Index of coeficients
        i,
        j,
        nc;

    int // enumerate orbitals
        k,
        l,
        s,
        q;

    int // auxiliar variables
        M2 = M * M,
        M3 = M * M * M;

    Iarray
        v;

    double
        sqrtOf;

    double complex
        z,
        w;

    nc = NC(N,M);
    v = iarrDef(M);



    for (i = 0; i < nc; i++)
    {
        w = 0;
        z = 0;

        for (k = 0; k < M; k++) v[k] = IF[k + i*M];
    
        /* ================================================================ *
         *                                                                  *
         *                       One-body contribution                      *
         *                                                                  *
         * ================================================================ */

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;

            w = w + Ho[k][k] * v[k] * C[i];

            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = sqrt((double)v[k] * (v[l] + 1));
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N,M,NCmat,v);
                w = w + Ho[k][l] * sqrtOf * C[j];
                v[k] += 1;
                v[l] -= 1;
            }
        }


        /* ================================================================ *
         *                                                                  *
         *                       Two-body contribution                      *
         *                                                                  *
         * ================================================================ */


        /* ---------------------------------------------
         * Rule 1: Creation on k k / Annihilation on k k
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            sqrtOf = v[k] * (v[k] - 1);
            z += Hint[k + M * k + M2 * k + M3 * k] * C[i] * sqrtOf;
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 2: Creation on k s / Annihilation on k s
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = v[k] * v[s];
                z += 4 * Hint[k + s*M + k*M2 + s*M3] * sqrtOf * C[i];
                /*
                z += Hint[s + k*M + k*M2 + s*M3] * sqrtOf * C[i];
                z += Hint[s + k*M + s*M2 + k*M3] * sqrtOf * C[i];
                z += Hint[k + s*M + s*M2 + k*M3] * sqrtOf * C[i];
                */
            }
        }
        /* ---------------------------------------------------------------- */

        // ---------------------------------------------
        // Rule 3: Creation on k k / Annihilation on q q
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < M; q++)
            {
                if (q == k) continue;
                sqrtOf = sqrt((double)(v[k]-1) * v[k] * (v[q]+1) * (v[q]+2));
                v[k] -= 2;
                v[q] += 2;
                j = FockToIndex(N, M, NCmat, v);
                z += Hint[k + k * M + q * M2 + q * M3] * C[j] * sqrtOf;
                v[k] += 2;
                v[q] -= 2;
            }
        }

        // ---------------------------------------------
        // Rule 4: Creation on k k / Annihilation on k l
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                z += 2 * Hint[k + k * M + k * M2 + l * M3] * C[j] * sqrtOf;
                // z += Hint[k + k * M + l * M2 + k * M3] * C[j] * sqrtOf;
                v[k] += 1;
                v[l] -= 1;
            }
        }


        // ---------------------------------------------
        // Rule 5: Creation on k s / Annihilation on s s
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = 0; s < M; s++)
            {
                if (s == k || v[s] < 1) continue;
                sqrtOf = v[s] * sqrt((double)v[k] * (v[s] + 1));
                v[k] -= 1;
                v[s] += 1;
                j = FockToIndex(N, M, NCmat, v);
                z += 2 * Hint[k + s * M + s * M2 + s * M3] * C[j] * sqrtOf;
                // z += Hint[s + k * M + s * M2 + s * M3] * C[j] * sqrtOf;
                v[k] += 1;
                v[s] -= 1;
            }
        }


        // -----------------------------------------------------------
        // Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = k + 1; q < M; q++)
            {
                for (l = q + 1; l < M; l++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
        // -------------------------------------------------------------------
        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 2) continue;
                for (l = k + 1; l < M; l++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
        // -------------------------------------------------------------------
        for (q = 0; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                for (k = l + 1; k < M; k++)
                {
                    if (v[k] < 2) continue;
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }



        // -----------------------------------------------------------
        // Rule 7.0: Creation on k s / Annihilation on q q (q > k > s)
        // -------------------------------------------------------------------
        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 1) continue;
                for (s = k + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 7.1: Creation on k s / Annihilation on q q (k > q > s)
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (q = k + 1; q < M; q++)
            {
                for (s = q + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 7.2: Creation on k s / Annihilation on q q (k > s > q)
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = s + 1; q < M; q++)
                {
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }



        // ---------------------------------------------
        // Rule 8: Creation on k s / Annihilation on s l
        // -------------------------------------------------------------------
        for (s = 0; s < M; s++)
        {
            if (v[s] < 1) continue; // may improve performance
            for (k = 0; k < M; k++)
            {
                if (v[k] < 1 || k == s) continue;
                for (l = 0; l < M; l++)
                {
                    if (l == k || l == s) continue;
                    sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    z += 4 * Hint[k + s*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    
                    //z += Hint[s + k*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    //z += Hint[s + k*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    //z += Hint[k + s*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    
                    v[k] += 1;
                    v[l] -= 1;
                }
            }
        }


        // ---------------------------------------------
        // Rule 9: Creation on k s / Annihilation on q l
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = 0; s < M; s++)
            {
                if (v[s] < 1 || s == k) continue;
                for (q = 0; q < M; q++)
                {
                    if (q == s || q == k) continue;
                    for (l = 0; l < M; l ++)
                    {
                        if (l == k || l == s || l == q) continue;
                        sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        z += Hint[k + s*M + q*M2 + l*M3] * C[j] * sqrtOf;
                        v[k] += 1;
                        v[s] += 1;
                        v[q] -= 1;
                        v[l] -= 1;
                    }   // Finish l
                }       // Finish q
            }           // Finish s
        }               // Finish k



        out[i] = w + 0.5 * z;
    }

    free(v);

}





void applyHconf_X (int N, int M, Iarray Map, Iarray MapOT, Iarray strideOT,
     Iarray NCmat, Iarray IF, Carray C, Cmatrix Ho, Carray Hint, Carray out)
{
    // Apply the many-body hamiltonian in a state expressed in
    // number-occupation basis with coefficients defined by C.



    int // Index of coeficients
        i,
        j,
        nc,
        chunks;

    int // enumerate orbitals
        k,
        l,
        s,
        q;

    int // auxiliar variables
        strideOrb,
        M2 = M * M,
        M3 = M * M * M;

    Iarray
        v;

    double
        sqrtOf;

    double complex
        z,
        w;

    nc = NC(N,M);
    v = iarrDef(M);



    for (i = 0; i < nc; i++)
    {
        w = 0;
        z = 0;

        for (k = 0; k < M; k++) v[k] = IF[k + i*M];
    
        /* ================================================================ *
         *                                                                  *
         *                       One-body contribution                      *
         *                                                                  *
         * ================================================================ */

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;

            w = w + Ho[k][k] * v[k] * C[i];

            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = sqrt((double)v[k] * (v[l] + 1));
                j = Map[i + k * nc + l * M * nc];
                //v[k] -= 1;
                //v[l] += 1;
                //j = FockToIndex(N, M, NCmat, v);
                w = w + Ho[k][l] * sqrtOf * C[j];
                //v[k] += 1;
                //v[l] -= 1;
            }
        }


        /* ================================================================ *
         *                                                                  *
         *                       Two-body contribution                      *
         *                                                                  *
         * ================================================================ */


        /* ---------------------------------------------
         * Rule 1: Creation on k k / Annihilation on k k
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            sqrtOf = v[k] * (v[k] - 1);
            z += Hint[k + M * k + M2 * k + M3 * k] * C[i] * sqrtOf;
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 2: Creation on k s / Annihilation on k s
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = v[k] * v[s];
                z += 4 * Hint[k + s*M + k*M2 + s*M3] * sqrtOf * C[i];
                /*
                z += Hint[s + k*M + k*M2 + s*M3] * sqrtOf * C[i];
                z += Hint[s + k*M + s*M2 + k*M3] * sqrtOf * C[i];
                z += Hint[k + s*M + s*M2 + k*M3] * sqrtOf * C[i];
                */
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 3: Creation on k k / Annihilation on q q
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < M; q++)
            {
                if (q == k) continue;
                sqrtOf = sqrt((double)(v[k]-1) * v[k] * (v[q]+1) * (v[q]+2));

                chunks = 0;
                for (j = 0; j < k; j++)
                {
                    if (v[j] > 1) chunks++;
                }
                strideOrb = chunks * M * M;
                j = MapOT[strideOT[i] + strideOrb + q + q * M];

                //v[k] -= 2;
                //v[q] += 2;
                // j = FockToIndex(N, M, NCmat, v);
                z += Hint[k + k * M + q * M2 + q * M3] * C[j] * sqrtOf;
                //v[k] += 2;
                //v[q] -= 2;
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 4: Creation on k k / Annihilation on k l
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));
                j = Map[i + k * nc + l * M * nc];
                //v[k] -= 1;
                //v[l] += 1;
                //j = FockToIndex(N, M, NCmat, v);
                z += 2 * Hint[k + k * M + k * M2 + l * M3] * C[j] * sqrtOf;
                // z += Hint[k + k * M + l * M2 + k * M3] * C[j] * sqrtOf;
                //v[k] += 1;
                //v[l] -= 1;
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 5: Creation on k s / Annihilation on s s
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = 0; s < M; s++)
            {
                if (s == k || v[s] < 1) continue;
                sqrtOf = v[s] * sqrt((double)v[k] * (v[s] + 1));
                j = Map[i + k * nc + s * M * nc];
                //v[k] -= 1;
                //v[s] += 1;
                //j = FockToIndex(N, M, NCmat, v);
                z += 2 * Hint[k + s * M + s * M2 + s * M3] * C[j] * sqrtOf;
                // z += Hint[s + k * M + s * M2 + s * M3] * C[j] * sqrtOf;
                //v[k] += 1;
                //v[s] -= 1;
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = k + 1; q < M; q++)
            {
                for (l = q + 1; l < M; l++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (v[j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;
                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    //v[k] -= 2;
                    //v[l] += 1;
                    //v[q] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    //v[k] += 2;
                    //v[l] -= 1;
                    //v[q] -= 1;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
        ------------------------------------------------------------------- */
        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 2) continue;
                for (l = k + 1; l < M; l++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (v[j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;
                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    //v[k] -= 2;
                    //v[l] += 1;
                    //v[q] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    //v[k] += 2;
                    //v[l] -= 1;
                    //v[q] -= 1;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
        ------------------------------------------------------------------- */
        for (q = 0; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                for (k = l + 1; k < M; k++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (v[j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;
                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    //v[k] -= 2;
                    //v[l] += 1;
                    //v[q] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    //v[k] += 2;
                    //v[l] -= 1;
                    //v[q] -= 1;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 7.0: Creation on k s / Annihilation on q q (q > k > s)
        ------------------------------------------------------------------- */
        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 1) continue;
                for (s = k + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));

                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);

                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;

                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 7.1: Creation on k s / Annihilation on q q (k > q > s)
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (q = k + 1; q < M; q++)
            {
                for (s = q + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));

                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);

                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;

                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 7.2: Creation on k s / Annihilation on q q (k > s > q)
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = s + 1; q < M; q++)
                {
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));

                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);

                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;

                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 8: Creation on k s / Annihilation on s l
        ------------------------------------------------------------------- */
        for (s = 0; s < M; s++)
        {
            if (v[s] < 1) continue; // may improve performance
            for (k = 0; k < M; k++)
            {
                if (v[k] < 1 || k == s) continue;
                for (l = 0; l < M; l++)
                {
                    if (l == k || l == s) continue;
                    sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                    j = Map[i + k * nc + l * M * nc];
                    //v[k] -= 1;
                    //v[l] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 4 * Hint[k + s*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    /*
                    z += Hint[s + k*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    z += Hint[s + k*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    z += Hint[k + s*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    */
                    //v[k] += 1;
                    //v[l] -= 1;
                }
            }
        }


        /* ---------------------------------------------
         * Rule 9: Creation on k s / Annihilation on q l
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = 0; q < M; q++)
                {
                    if (q == s || q == k) continue;
                    for (l = q + 1; l < M; l ++)
                    {
                        if (l == k || l == s) continue;
                        sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));

                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);

                        z += 4 * Hint[k + s*M + q*M2 + l*M3] * C[j] * sqrtOf;

                        v[k] += 1;
                        v[s] += 1;
                        v[q] -= 1;
                        v[l] -= 1;

                    }   // Finish l
                }       // Finish q
            }           // Finish s
        }               // Finish k

        out[i] = w + 0.5 * z;
    }

    free(v);

}




void applyHconf_XX (int N, int M, Iarray Map, Iarray MapOT, Iarray MapTT,
     Iarray strideOT, Iarray strideTT, Iarray IF, Carray C, Cmatrix Ho,
     Carray Hint, Carray out)
{
    // Apply the many-body hamiltonian in a state expressed in
    // number-occupation basis with coefficients defined by C.



    int // Index of coeficients
        i,
        j,
        nc,
        chunks;

    int // enumerate orbitals
        h,
        g,
        k,
        l,
        s,
        q;

    int // auxiliar variables
        strideOrb,
        M2 = M * M,
        M3 = M * M * M;

    Iarray
        v;

    double
        sqrtOf;

    double complex
        z,
        w;

    nc = NC(N,M);


#pragma omp parallel private(i,j,k,s,q,l,h,g,chunks,strideOrb,z,w,sqrtOf,v)
    {

    v = iarrDef(M);

#pragma omp for schedule(static)
    for (i = 0; i < nc; i++)
    {
        w = 0;
        z = 0;

        for (k = 0; k < M; k++) v[k] = IF[k + i*M];
    
        /* ================================================================ *
         *                                                                  *
         *                       One-body contribution                      *
         *                                                                  *
         * ================================================================ */

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;

            w = w + Ho[k][k] * v[k] * C[i];

            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = sqrt((double)v[k] * (v[l] + 1));
                j = Map[i + k * nc + l * M * nc];
                //v[k] -= 1;
                //v[l] += 1;
                //j = FockToIndex(N, M, NCmat, v);
                w = w + Ho[k][l] * sqrtOf * C[j];
                //v[k] += 1;
                //v[l] -= 1;
            }
        }


        /* ================================================================ *
         *                                                                  *
         *                       Two-body contribution                      *
         *                                                                  *
         * ================================================================ */


        /* ---------------------------------------------
         * Rule 1: Creation on k k / Annihilation on k k
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            sqrtOf = v[k] * (v[k] - 1);
            z += Hint[k + M * k + M2 * k + M3 * k] * C[i] * sqrtOf;
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 2: Creation on k s / Annihilation on k s
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = v[k] * v[s];
                z += 4 * Hint[k + s*M + k*M2 + s*M3] * sqrtOf * C[i];
                /*
                z += Hint[s + k*M + k*M2 + s*M3] * sqrtOf * C[i];
                z += Hint[s + k*M + s*M2 + k*M3] * sqrtOf * C[i];
                z += Hint[k + s*M + s*M2 + k*M3] * sqrtOf * C[i];
                */
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 3: Creation on k k / Annihilation on q q
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < M; q++)
            {
                if (q == k) continue;
                sqrtOf = sqrt((double)(v[k]-1) * v[k] * (v[q]+1) * (v[q]+2));

                chunks = 0;
                for (j = 0; j < k; j++)
                {
                    if (v[j] > 1) chunks++;
                }
                strideOrb = chunks * M * M;
                j = MapOT[strideOT[i] + strideOrb + q + q * M];

                //v[k] -= 2;
                //v[q] += 2;
                // j = FockToIndex(N, M, NCmat, v);
                z += Hint[k + k * M + q * M2 + q * M3] * C[j] * sqrtOf;
                //v[k] += 2;
                //v[q] -= 2;
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 4: Creation on k k / Annihilation on k l
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));
                j = Map[i + k * nc + l * M * nc];
                //v[k] -= 1;
                //v[l] += 1;
                //j = FockToIndex(N, M, NCmat, v);
                z += 2 * Hint[k + k * M + k * M2 + l * M3] * C[j] * sqrtOf;
                // z += Hint[k + k * M + l * M2 + k * M3] * C[j] * sqrtOf;
                //v[k] += 1;
                //v[l] -= 1;
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 5: Creation on k s / Annihilation on s s
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = 0; s < M; s++)
            {
                if (s == k || v[s] < 1) continue;
                sqrtOf = v[s] * sqrt((double)v[k] * (v[s] + 1));
                j = Map[i + k * nc + s * M * nc];
                //v[k] -= 1;
                //v[s] += 1;
                //j = FockToIndex(N, M, NCmat, v);
                z += 2 * Hint[k + s * M + s * M2 + s * M3] * C[j] * sqrtOf;
                // z += Hint[s + k * M + s * M2 + s * M3] * C[j] * sqrtOf;
                //v[k] += 1;
                //v[s] -= 1;
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = k + 1; q < M; q++)
            {
                for (l = q + 1; l < M; l++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (v[j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;
                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    //v[k] -= 2;
                    //v[l] += 1;
                    //v[q] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    //v[k] += 2;
                    //v[l] -= 1;
                    //v[q] -= 1;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
        ------------------------------------------------------------------- */
        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 2) continue;
                for (l = k + 1; l < M; l++)
                {
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (v[j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;
                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    //v[k] -= 2;
                    //v[l] += 1;
                    //v[q] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    //v[k] += 2;
                    //v[l] -= 1;
                    //v[q] -= 1;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
        ------------------------------------------------------------------- */
        for (q = 0; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                for (k = l + 1; k < M; k++)
                {

                    if (v[k] < 2) continue;

                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (v[j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;
                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    //v[k] -= 2;
                    //v[l] += 1;
                    //v[q] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 2 * Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    // z += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    //v[k] += 2;
                    //v[l] -= 1;
                    //v[q] -= 1;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 7.0: Creation on k s / Annihilation on q q (q > k > s)
        ------------------------------------------------------------------- */
        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 1) continue;
                for (s = k + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));

                    chunks = 0;
                    for (h = 0; h < k; h++)
                    {
                        for (g = h + 1; g < M; g++)
                        {
                            if (v[h] > 0 && v[g] > 0) chunks++;
                        }
                    }

                    for (g = k + 1; g < s; g++)
                    {
                        if (v[g] > 0) chunks++;
                    }

                    strideOrb = chunks * M * M;

                    j = MapTT[strideTT[i] + strideOrb + q + q*M];

                    //v[k] -= 1;
                    //v[s] -= 1;
                    //v[q] += 2;
                    //j = FockToIndex(N, M, NCmat, v);

                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;

                    //v[k] += 1;
                    //v[s] += 1;
                    //v[q] -= 2;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 7.1: Creation on k s / Annihilation on q q (k > q > s)
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (q = k + 1; q < M; q++)
            {
                for (s = q + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));

                    chunks = 0;
                    for (h = 0; h < k; h++)
                    {
                        for (g = h + 1; g < M; g++)
                        {
                            if (v[h] > 0 && v[g] > 0) chunks++;
                        }
                    }

                    for (g = k + 1; g < s; g++)
                    {
                        if (v[g] > 0) chunks++;
                    }

                    strideOrb = chunks * M * M;

                    j = MapTT[strideTT[i] + strideOrb + q + q*M];

                    //v[k] -= 1;
                    //v[s] -= 1;
                    //v[q] += 2;
                    //j = FockToIndex(N, M, NCmat, v);

                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;

                    //v[k] += 1;
                    //v[s] += 1;
                    //v[q] -= 2;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* -----------------------------------------------------------
         * Rule 7.2: Creation on k s / Annihilation on q q (k > s > q)
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = s + 1; q < M; q++)
                {
                    sqrtOf = sqrt((double)v[k] * v[s] * (v[q]+1) * (v[q]+2));

                    chunks = 0;
                    for (h = 0; h < k; h++)
                    {
                        for (g = h + 1; g < M; g++)
                        {
                            if (v[h] > 0 && v[g] > 0) chunks++;
                        }
                    }

                    for (g = k + 1; g < s; g++)
                    {
                        if (v[g] > 0) chunks++;
                    }

                    strideOrb = chunks * M * M;

                    j = MapTT[strideTT[i] + strideOrb + q + q*M];

                    //v[k] -= 1;
                    //v[s] -= 1;
                    //v[q] += 2;
                    //j = FockToIndex(N, M, NCmat, v);

                    z += 2 * Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    // z += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;

                    //v[k] += 1;
                    //v[s] += 1;
                    //v[q] -= 2;
                }
            }
        }
        /* ---------------------------------------------------------------- */


        /* ---------------------------------------------
         * Rule 8: Creation on k s / Annihilation on s l
        ------------------------------------------------------------------- */
        for (s = 0; s < M; s++)
        {
            if (v[s] < 1) continue; // may improve performance
            for (k = 0; k < M; k++)
            {
                if (v[k] < 1 || k == s) continue;
                for (l = 0; l < M; l++)
                {
                    if (l == k || l == s) continue;
                    sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));
                    j = Map[i + k * nc + l * M * nc];
                    //v[k] -= 1;
                    //v[l] += 1;
                    //j = FockToIndex(N, M, NCmat, v);
                    z += 4 * Hint[k + s*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    /*
                    z += Hint[s + k*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    z += Hint[s + k*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    z += Hint[k + s*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    */
                    //v[k] += 1;
                    //v[l] -= 1;
                }
            }
        }


        /* ---------------------------------------------
         * Rule 9: Creation on k s / Annihilation on q l
        ------------------------------------------------------------------- */
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = 0; q < M; q++)
                {
                    if (q == s || q == k) continue;
                    for (l = q + 1; l < M; l ++)
                    {
                        if (l == k || l == s) continue;
                        sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));

                        chunks = 0;
                        for (h = 0; h < k; h++)
                        {
                            for (g = h + 1; g < M; g++)
                            {
                                if (v[h] > 0 && v[g] > 0) chunks++;
                            }
                        }

                        for (g = k + 1; g < s; g++)
                        {
                            if (v[g] > 0) chunks++;
                        }

                        strideOrb = chunks * M * M;

                        j = MapTT[strideTT[i] + strideOrb + q + l*M];

                        //v[k] -= 1;
                        //v[s] -= 1;
                        //v[q] += 1;
                        //v[l] += 1;
                        //j = FockToIndex(N, M, NCmat, v);

                        z += 4 * Hint[k + s*M + q*M2 + l*M3] * C[j] * sqrtOf;

                        //v[k] += 1;
                        //v[s] += 1;
                        //v[q] -= 1;
                        //v[l] -= 1;

                    }   // Finish l
                }       // Finish q
            }           // Finish s
        }               // Finish k

        out[i] = w + 0.5 * z;
    }

    free(v);

    } // end of parallel region

}



#endif
