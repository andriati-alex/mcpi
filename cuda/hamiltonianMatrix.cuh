#ifndef _hamiltonianMatrix_cuh
#define _hamiltonianMatrix_cuh

#include "configurationsMap.cuh"



__device__
void d_NC(int N, int M, int * nc)
{

/** Number of Configurations(NC) of N particles in M states **/

    long
        j,
        i,
        n;

    n = 1;
    j = 2;

    if  (M > N)
    {
        for (i = N + M - 1; i > M - 1; i --)
        {
            n = n * i;
            if (n % j == 0 && j <= N)
            {
                n = n / j;
                j = j + 1;
            }
        }

        for (i = j; i < N + 1; i++) n = n / i;

        *nc = ((int) n);
    }

    else
    {

        for (i = N + M - 1; i > N; i --)
        {
            n = n * i;
            if (n % j == 0 && j <= M - 1)
            {
                n = n / j;
                j = j + 1;
            }
        }

        for (i = j; i < M; i++) n = n / i;

        *nc =  ((int) n);
    }
}



__global__
void applyHconf (int N, int M, Iarray Map, Iarray MapOT, Iarray MapTT,
     Iarray strideOT, Iarray strideTT, Iarray IF, Carray C, Carray Ho,
     Carray Hint, Carray out)
{
    // Apply the many-body hamiltonian in a state expressed in
    // number-occupation basis with coefficients defined by C.


    int
        i,
        j,
        nc,
        hIndex,
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
        cuStride,
        init,
        M2,
        M3,
        vk,
        vl,
        vs,
        vq;

    cuDoubleComplex
        z,
        w,
        aux,
        sqrtOf;

    cuStride = blockDim.x * gridDim.x;
    init = blockIdx.x * blockDim.x + threadIdx.x;


    d_NC(N,M,&nc);
    M2 = M * M;
    M3 = M * M2;

    for (i = init; i < nc; i += cuStride)
    {
        w = make_cuDoubleComplex(0,0);
        z = make_cuDoubleComplex(0,0);

        /* ================================================================ *
         *                                                                  *
         *                       One-body contribution                      *
         *                                                                  *
         * ================================================================ */

        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];

            if (vk < 1) continue;

            aux = cuCmul(Ho[k+M*k],
                    cuCmul(C[i],make_cuDoubleComplex(vk,0)));

            w = cuCadd(w,aux);

            for (l = 0; l < M; l++)
            {
                if (l == k) continue;

                vl = IF[l + i*M];

                sqrtOf = make_cuDoubleComplex(sqrt((double)vk*(vl+ 1)),0);

                j = Map[i + k * nc + l * M * nc];

                aux = cuCmul(Ho[k+M*l],cuCmul(sqrtOf,C[j]));
                w = cuCadd(w,aux);
            }
        }



        /* ================================================================ *
         *                                                                  *
         *                       Two-body contribution                      *
         *                                                                  *
         * ================================================================ */


        // ---------------------------------------------
        // Rule 1: Creation on k k / Annihilation on k k
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            hIndex = k + M * k + M2 * k + M3 * k;
            sqrtOf = make_cuDoubleComplex(vk * (vk - 1),0);
            aux = cuCmul(Hint[hIndex],cuCmul(C[i],sqrtOf));
            z = cuCadd(z,aux);
        }


        // ---------------------------------------------
        // Rule 2: Creation on k s / Annihilation on k s
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 1) continue;

            for (s = k + 1; s < M; s++)
            {
                vs = IF[s + i*M];
                hIndex = k + s * M + k * M2 + s * M3;
                sqrtOf = make_cuDoubleComplex(4 * vk * vs,0);
                aux = cuCmul(Hint[hIndex],cuCmul(sqrtOf,C[i]));
                z = cuCadd(z,aux);
                // z += 4 * Hint[k + s*M + k*M2 + s*M3] * sqrtOf * C[i];
            }
        }


        // ---------------------------------------------
        // Rule 3: Creation on k k / Annihilation on q q
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 2) continue;

            for (q = 0; q < M; q++)
            {
                if (q == k) continue;

                vq = IF[q + i*M];

                hIndex = k + k * M + q * M2 + q * M3;

                sqrtOf = make_cuDoubleComplex(
                        sqrt((double)(vk-1)*vk*(vq+1)*(vq+2)),0);

                chunks = 0;
                for (j = 0; j < k; j++) { if (IF[j+i*M] > 1) chunks++; }
                strideOrb = chunks * M * M;

                j = MapOT[strideOT[i] + strideOrb + q + q * M];

                aux = cuCmul(Hint[hIndex],cuCmul(sqrtOf,C[j]));
                z = cuCadd(z,aux);
            }
        }




        // ---------------------------------------------
        // Rule 4: Creation on k k / Annihilation on k l
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 2) continue;

            for (l = 0; l < M; l++)
            {

                if (l == k) continue;
                vl = IF[l + i*M];

                hIndex = k + k * M + k * M2 + l * M3;

                sqrtOf = make_cuDoubleComplex(
                        2*(vk-1)*sqrt((double)vk*(vl+1)),0);

                j = Map[i + k * nc + l * M * nc];

                aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));
                z = cuCadd(z,aux);
            }
        }


        // ---------------------------------------------
        // Rule 5: Creation on k s / Annihilation on s s
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 1) continue;

            for (s = 0; s < M; s++)
            {
                vs = IF[s + i*M];
                if (s == k || vs < 1) continue;

                hIndex = k + s * M + s * M2 + s * M3;

                sqrtOf = make_cuDoubleComplex(2*vs*sqrt((double)vk*(vs+1)),0);

                j = Map[i + k * nc + s * M * nc];

                aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));
                z = cuCadd(z,aux);
            }
        }


        // -----------------------------------------------------------
        // Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 2) continue;

            for (q = k + 1; q < M; q++)
            {
                vq = IF[q + i*M];

                for (l = q + 1; l < M; l++)
                {
                    vl = IF[l + i*M];
                    hIndex = k + k * M + q * M2 + l * M3;

                    sqrtOf = make_cuDoubleComplex(
                            2*sqrt((double)vk*(vk-1)*(vq+1)*(vl+1)),0);

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[j + i*M] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
        // -------------------------------------------------------------------
        for (q = 0; q < M; q++)
        {
            vq = IF[q + i*M];

            for (k = q + 1; k < M; k++)
            {
                vk = IF[k + i*M];
                if (vk < 2) continue;

                for (l = k + 1; l < M; l++)
                {
                    vl = IF[l + i*M];
                    hIndex = k + k * M + q * M2 + l * M3;

                    sqrtOf = make_cuDoubleComplex(
                            2*sqrt((double)vk*(vk-1)*(vq+1)*(vl+1)),0);

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[j + i*M] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
        // -------------------------------------------------------------------
        for (q = 0; q < M; q++)
        {
            vq = IF[q + i*M];

            for (l = q + 1; l < M; l++)
            {
                vl = IF[l + i*M];

                for (k = l + 1; k < M; k++)
                {
                    vk = IF[k + i*M];

                    if (vk < 2) continue;

                    hIndex = k + k * M + q * M2 + l * M3;

                    sqrtOf = make_cuDoubleComplex(
                            2*sqrt((double)vk*(vk-1)*(vq+1)*(vl+1)),0);

                    chunks = 0;
                    for (j = 0; j < k; j++)
                    {
                        if (IF[j + i*M] > 1) chunks++;
                    }
                    strideOrb = chunks * M * M;

                    j = MapOT[strideOT[i] + strideOrb + q + l * M];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }



        // -----------------------------------------------------------
        // Rule 7.0: Creation on k s / Annihilation on q q (q > k > s)
        // -------------------------------------------------------------------
        for (q = 0; q < M; q++)
        {
            vq = IF[q + i*M];

            for (k = q + 1; k < M; k++)
            {

                vk = IF[k + i*M];
                if (vk < 1) continue;

                for (s = k + 1; s < M; s++)
                {
                    vs = IF[s + i*M];
                    if (vs < 1) continue;

                    hIndex = k + s * M + q * M2 + q * M3;

                    sqrtOf = make_cuDoubleComplex(
                            2*sqrt((double)vk*vs*(vq+1)*(vq+2)),0);

                    chunks = 0;
                    for (h = 0; h < k; h++)
                    {
                        for (g = h + 1; g < M; g++)
                        {
                            if (IF[h+i*M] > 0 && IF[g+i*M] > 0) chunks++;
                        }
                    }

                    for (g = k + 1; g < s; g++)
                    {
                        if (IF[g + i*M] > 0) chunks++;
                    }

                    strideOrb = chunks * M * M;

                    j = MapTT[strideTT[i] + strideOrb + q + q*M];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 7.1: Creation on k s / Annihilation on q q (k > q > s)
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 1) continue;

            for (q = k + 1; q < M; q++)
            {
                vq = IF[q + i*M];
                for (s = q + 1; s < M; s++)
                {
                    vs = IF[s + i*M];
                    if (vs < 1) continue;

                    hIndex = k + s * M + q * M2 + q * M3;

                    sqrtOf = make_cuDoubleComplex(
                            2*sqrt((double)vk*vs*(vq+1)*(vq+2)),0);

                    chunks = 0;
                    for (h = 0; h < k; h++)
                    {
                        for (g = h + 1; g < M; g++)
                        {
                            if (IF[h+i*M] > 0 && IF[g+i*M] > 0) chunks++;
                        }
                    }

                    for (g = k + 1; g < s; g++)
                    {
                        if (IF[g + i*M] > 0) chunks++;
                    }

                    strideOrb = chunks * M * M;

                    j = MapTT[strideTT[i] + strideOrb + q + q*M];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }


        // -----------------------------------------------------------
        // Rule 7.2: Creation on k s / Annihilation on q q (k > s > q)
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 1) continue;

            for (s = k + 1; s < M; s++)
            {
                vs = IF[s + i*M];
                if (vs < 1) continue;

                for (q = s + 1; q < M; q++)
                {
                    vq = IF[q + i*M];

                    hIndex = k + s * M + q * M2 + q * M3;

                    sqrtOf = make_cuDoubleComplex(
                            2*sqrt((double)vk*vs*(vq+1)*(vq+2)),0);

                    chunks = 0;
                    for (h = 0; h < k; h++)
                    {
                        for (g = h + 1; g < M; g++)
                        {
                            if (IF[h+i*M] > 0 && IF[g+i*M] > 0) chunks++;
                        }
                    }

                    for (g = k + 1; g < s; g++)
                    {
                        if (IF[g + i*M] > 0) chunks++;
                    }

                    strideOrb = chunks * M * M;

                    j = MapTT[strideTT[i] + strideOrb + q + q*M];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }


        // ---------------------------------------------
        // Rule 8: Creation on k s / Annihilation on s l
        // -------------------------------------------------------------------
        for (s = 0; s < M; s++)
        {
            vs = IF[s + i*M];
            if (vs < 1) continue;

            for (k = 0; k < M; k++)
            {
                vk = IF[k + i*M];
                if (vk < 1 || k == s) continue;

                for (l = 0; l < M; l++)
                {
                    if (l == k || l == s) continue;

                    vl = IF[l + i*M];

                    hIndex = k + s * M + s * M2 + l * M3;

                    sqrtOf = make_cuDoubleComplex(
                            4*vs*sqrt((double)vk*(vl+1)),0);

                    j = Map[i + k * nc + l * M * nc];

                    aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                    z = cuCadd(z,aux);
                }
            }
        }


        // ---------------------------------------------
        // Rule 9: Creation on k s / Annihilation on q l
        // -------------------------------------------------------------------
        for (k = 0; k < M; k++)
        {
            vk = IF[k + i*M];
            if (vk < 1) continue;

            for (s = k + 1; s < M; s++)
            {
                vs = IF[s + i*M];
                if (vs < 1) continue;

                for (q = 0; q < M; q++)
                {
                    if (q == s || q == k) continue;
                    vq = IF[q + i*M];

                    for (l = q + 1; l < M; l ++)
                    {
                        if (l == k || l == s) continue;
                        vl = IF[l + i*M];

                        hIndex = k + s * M + q * M2 + l * M3;

                        sqrtOf = make_cuDoubleComplex(
                                4*sqrt((double)vk*vs*(vq+1)*(vl+1)),0);

                        chunks = 0;
                        for (h = 0; h < k; h++)
                        {
                            for (g = h + 1; g < M; g++)
                            {
                                if (IF[h+i*M] > 0 && IF[g+i*M] > 0) chunks++;
                            }
                        }

                        for (g = k + 1; g < s; g++)
                        {
                            if (IF[g + i*M] > 0) chunks++;
                        }

                        strideOrb = chunks * M * M;

                        j = MapTT[strideTT[i] + strideOrb + q + l*M];

                        aux = cuCmul(Hint[hIndex],cuCmul(C[j],sqrtOf));

                        z = cuCadd(z,aux);
                    }   // Finish l
                }       // Finish q
            }           // Finish s
        }               // Finish k



        aux = cuCmul(z,make_cuDoubleComplex(0.5,0));

        out[i] = cuCadd(w,aux);
    }

}



#endif
