#ifndef _twobodyMatrix_h
#define _twobodyMatrix_h

#include <math.h>
#include "configurationsMap.h"



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



void TBrho(int N, int M, Iarray NCmat, Iarray IF, Carray C, Carray rho)
{

    int
        i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        t,
        nc,
        M2,
        M3;

    Iarray
        v;
    
    double
        mod2,   // |Cj| ^ 2
        sqrtOf; // Factors from the action of creation/annihilation

    double complex
        RHO;





    // Auxiliar to memory access
    M2 = M * M;
    M3 = M * M * M;

    nc = NC(N,M);
    v = iarrDef(M);



    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {

        RHO = 0;

        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO  = RHO + mod2 * IF[k + i*M] * (IF[k + i*M] - 1);
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

            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[k + i*M] * IF[s + i*M];
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

            for (i = 0; i < nc; i++)
            {
                if (IF[i*M + k] < 2) continue;
                for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                sqrtOf = sqrt((double)(v[k]-1)*v[k]*(v[q]+1)*(v[q]+2));

                v[k] -= 2;
                v[q] += 2;
                j = FockToIndex(N, M, NCmat, v);

                RHO += conj(C[i]) * C[j] * sqrtOf;
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

            for (i = 0; i < nc; i++)
            {
                if (IF[i*M + k] < 2) continue;
                for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                sqrtOf = (v[k] - 1) * sqrt((double)v[k] * (v[l] + 1));

                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);

                RHO += conj(C[i]) * C[j] * sqrtOf;
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

            for (i = 0; i < nc; i++)
            {
                if (IF[i*M + k] < 1 || IF[i*M + s] < 1) continue;
                for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                sqrtOf = v[s] * sqrt((double)v[k]*(v[s]+1));

                v[k] -= 1;
                v[s] += 1;
                j = FockToIndex(N, M, NCmat, v);

                RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    if (IF[i*M + k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    if (IF[i*M + k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    if (IF[i*M + k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                    sqrtOf = sqrt((double)v[k]*(v[k]-1)*(v[q]+1)*(v[l]+1));

                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    if (IF[i*M + k] < 1 || IF[i*M + s] < 1) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                    sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));

                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    if (IF[i*M + k] < 1 || IF[i*M + s] < 1) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                    sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));

                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    if (IF[i*M + k] < 1 || IF[i*M + s] < 1) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                    sqrtOf = v[s] * sqrt((double)v[k] * (v[l] + 1));

                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                    if (l == k || l == s) continue;

                    RHO = 0;

                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i*M + k] < 1 || IF[i*M + s] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i*M + t];
                        sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));

                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);

                        RHO += conj(C[i]) * C[j] * sqrtOf;
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





void TBrho_X(int N, int M, Iarray Map, Iarray MapOT, Iarray strideOT,
     Iarray NCmat, Iarray IF, Carray C, Carray rho)
{

    int i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        h,
        nc,
        M2,
        M3,
        vk,
        vs,
        vl,
        vq,
        chunks,
        strideOrb;

    Iarray
        v;

    double
        mod2,   // |Cj| ^ 2
        sqrtOf; // Factors from the action of creation/annihilation

    double complex
        RHO;





    // Auxiliar to memory access
    M2 = M * M;
    M3 = M * M * M;

    nc = NC(N,M);
    v = iarrDef(M);



    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {

        RHO = 0;

        for (i = 0; i < nc; i++)
        {
            if (IF[i*M + k] < 2) continue;
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO  = RHO + mod2 * IF[k + i*M] * (IF[k + i*M] - 1);
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

            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[k + i*M] * IF[s + i*M];
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

            for (i = 0; i < nc; i++)
            {
                h = i * M;
                if (IF[h + k] < 2) continue;
                vk = IF[h + k];
                vq = IF[h + q];

                chunks = 0;
                for (j = 0; j < k; j++)
                {
                    if (IF[h + j] > 1) chunks++;
                }
                strideOrb = chunks * M * M;

                j = MapOT[strideOT[i] + strideOrb + q + q * M];

                sqrtOf = sqrt((double)(vk-1)*vk*(vq+1)*(vq+2));
                RHO += conj(C[i]) * C[j] * sqrtOf;
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

            for (i = 0; i < nc; i++)
            {
                vk = IF[i*M + k];
                vl = IF[i*M + l];
                if (vk < 2) continue;

                j = Map[i + k * nc + l * M * nc];
                sqrtOf = (vk - 1) * sqrt((double) vk * (vl + 1));
                RHO += conj(C[i]) * C[j] * sqrtOf;
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

            for (i = 0; i < nc; i++)
            {
                vk = IF[i*M + k];
                vs = IF[i*M + s];
                if (vk < 1 || vs < 1) continue;

                j = Map[i + k * nc + s * M * nc];
                sqrtOf = vs * sqrt((double) vk * (vs + 1));
                RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    h = i * M; // auxiliat stride for IF
                    if (IF[h + k] < 2) continue;
                    vk = IF[h + k];
                    vq = IF[h + q];
                    vl = IF[h + l];

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[h + j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + l + q * M];

                    sqrtOf = sqrt((double)vk*(vk-1)*(vq+1)*(vl+1));

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    h = i * M;
                    if (IF[h + k] < 2) continue;
                    vk = IF[h + k];
                    vq = IF[h + q];
                    vl = IF[h + l];

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[h + j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + l + q * M];

                    sqrtOf = sqrt((double)vk*(vk-1)*(vq+1)*(vl+1));

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    h = i * M; // auxiliar stride for IF
                    if (IF[h + k] < 2) continue;
                    vk = IF[h + k];
                    vq = IF[h + q];
                    vl = IF[h + l];

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[h + j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + l + q * M];

                    sqrtOf = sqrt((double)vk*(vk-1)*(vq+1)*(vl+1));

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    vs = IF[i*M + s];
                    vk = IF[i*M + k];
                    vl = IF[i*M + l];
                    if (vk < 1 || vs < 1) continue;
                    j = Map[i + k * nc + l * M * nc];
                    sqrtOf = vs * sqrt((double) vk * (vl + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    vs = IF[i*M + s];
                    vk = IF[i*M + k];
                    vl = IF[i*M + l];
                    if (vk < 1 || vs < 1) continue;
                    j = Map[i + k * nc + l * M * nc];
                    sqrtOf = vs * sqrt((double) vk * (vl + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    vs = IF[i*M + s];
                    vk = IF[i*M + k];
                    vl = IF[i*M + l];
                    if (vk < 1 || vs < 1) continue;
                    j = Map[i + k * nc + l * M * nc];
                    sqrtOf = vs * sqrt((double) vk * (vl + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                    if (l == k || l == s) continue;

                    RHO = 0;

                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i*M + k] < 1 || IF[i*M + s] < 1) continue;
                        for (h = 0; h < M; h++) v[h] = IF[i*M + h];

                        sqrtOf = sqrt((double)v[k]*v[s]*(v[q]+1)*(v[l]+1));
                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N,M,NCmat,v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }

                    rho[k + s * M + q * M2 + l * M3] = RHO;
                    rho[s + k * M + q * M2 + l * M3] = RHO;
                    rho[k + s * M + l * M2 + q * M3] = RHO;
                    rho[s + k * M + l * M2 + q * M3] = RHO;
                }   // Finish l loop
            }       // Finish q loop
        }           // Finish s loop
    }               // Finish k loop

    free(v);

    /*       ------------------- END OF ROUTINE -------------------       */
}





void TBrho_XX(int N, int M, Iarray Map, Iarray MapOT, Iarray MapTT,
     Iarray strideOT, Iarray strideTT, Iarray NCmat, Iarray IF,
     Carray C, Carray rho)
{

    int i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        h,
        g,
        nc,
        M2,
        M3,
        vk,
        vs,
        vl,
        vq,
        chunks,
        strideOrb;

    double
        mod2,   // |Cj| ^ 2
        sqrtOf; // Factors from the action of creation/annihilation

    double complex
        RHO;





    // Auxiliar to memory access
    M2 = M * M;
    M3 = M * M * M;

    nc = NC(N,M);



    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {

        RHO = 0;

        for (i = 0; i < nc; i++)
        {
            if (IF[k + i*M] < 2) continue;
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO  = RHO + mod2 * IF[k + i*M] * (IF[k + i*M] - 1);
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

            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[k + i*M] * IF[s + i*M];
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

            for (i = 0; i < nc; i++)
            {
                h = i * M; // auxiliar stride for IF
                if (IF[h + k] < 2) continue;
                vk = IF[h + k];
                vq = IF[h + q];

                chunks = 0;
                for (j = 0; j < k; j++)
                {
                    if (IF[h + j] > 1) chunks++;
                }
                strideOrb = chunks * M * M;

                j = MapOT[strideOT[i] + strideOrb + q + q * M];

                sqrtOf = sqrt((double)(vk-1)*vk*(vq+1)*(vq+2));
                RHO += conj(C[i]) * C[j] * sqrtOf;
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

            for (i = 0; i < nc; i++)
            {
                vk = IF[i*M + k];
                vl = IF[i*M + l];
                if (vk < 2) continue;

                j = Map[i + k * nc + l * M * nc];
                sqrtOf = (vk - 1) * sqrt((double) vk * (vl + 1));
                RHO += conj(C[i]) * C[j] * sqrtOf;
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

            for (i = 0; i < nc; i++)
            {
                vk = IF[i*M + k];
                vs = IF[i*M + s];
                if (vk < 1 || vs < 1) continue;

                j = Map[i + k * nc + s * M * nc];
                sqrtOf = vs * sqrt((double) vk * (vs + 1));
                RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    h = i * M; // auxiliat stride for IF
                    if (IF[h + k] < 2) continue;
                    vk = IF[h + k];
                    vq = IF[h + q];
                    vl = IF[h + l];

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[h + j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + l + q * M];

                    sqrtOf = sqrt((double)vk*(vk-1)*(vq+1)*(vl+1));

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    h = i * M;
                    if (IF[h + k] < 2) continue;
                    vk = IF[h + k];
                    vq = IF[h + q];
                    vl = IF[h + l];

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[h + j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + l + q * M];

                    sqrtOf = sqrt((double)vk*(vk-1)*(vq+1)*(vl+1));

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    h = i * M; // auxiliar stride for IF
                    if (IF[h + k] < 2) continue;
                    vk = IF[h + k];
                    vq = IF[h + q];
                    vl = IF[h + l];

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[h + j] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + l + q * M];

                    sqrtOf = sqrt((double)vk*(vk-1)*(vq+1)*(vl+1));

                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    vs = IF[i*M + s];
                    vk = IF[i*M + k];
                    vl = IF[i*M + l];
                    if (vk < 1 || vs < 1) continue;
                    j = Map[i + k * nc + l * M * nc];
                    sqrtOf = vs * sqrt((double) vk * (vl + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    vs = IF[i*M + s];
                    vk = IF[i*M + k];
                    vl = IF[i*M + l];
                    if (vk < 1 || vs < 1) continue;
                    j = Map[i + k * nc + l * M * nc];
                    sqrtOf = vs * sqrt((double) vk * (vl + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                for (i = 0; i < nc; i++)
                {
                    vs = IF[i*M + s];
                    vk = IF[i*M + k];
                    vl = IF[i*M + l];
                    if (vk < 1 || vs < 1) continue;
                    j = Map[i + k * nc + l * M * nc];
                    sqrtOf = vs * sqrt((double) vk * (vl + 1));
                    RHO += conj(C[i]) * C[j] * sqrtOf;
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

                    if (l == k || l == s) continue;

                    RHO = 0;

                    for (i = 0; i < nc; i++)
                    {
                        j = i * M; // auxiliar stride for IF
                        if (IF[j + k] < 1 || IF[j + s] < 1) continue;

                        vk = IF[j + k];
                        vl = IF[j + l];
                        vq = IF[j + q];
                        vs = IF[j + s];

                        chunks = 0;
                        for (h = 0; h < k; h++)
                        {
                            for (g = h + 1; g < M; g++)
                            {
                                if (IF[h+j] > 0 && IF[g+j] > 0) chunks++;
                            }
                        }

                        for (g = k + 1; g < s; g++)
                        {
                            if (IF[g+j] > 0) chunks++;
                        }

                        strideOrb = chunks * M * M;

                        sqrtOf = sqrt((double)vk*vs*(vq+1)*(vl+1));

                        j = MapTT[strideTT[i] + strideOrb + q + l*M];

                        RHO += conj(C[i]) * C[j] * sqrtOf;
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



#endif
