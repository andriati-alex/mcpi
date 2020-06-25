



void setupH(int N, int lmax, int mcsize, Iarray * ht, Carray Ho, Carray Hint)
{



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

    M = 2 * lmax + 1;   // total number of IPS
    v = iarrDef(M);     // vector of occupation numbers

    for (i = 0; i < mcsize; i++)
    {
        w = 0;
        z = 0;

        for (k = 0; k < M; k++) v[k] = ht[i][M];

        /* ================================================================ *
         *                                                                  *
         *                       One-body contribution                      *
         *                                                                  *
         * ================================================================ */

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            w = w + Ho[k] * v[k] * C[i];
        }



        /* ================================================================ *
         *                                                                  *
         *                       Two-body contribution                      *
         *                                                                  *
         * ================================================================ */

        // Rule 1: Creation on k k / Annihilation on k k
        for (k = 0; k < M; k++)
        {
            sqrtOf = v[k] * (v[k] - 1);
            z = z + Hint[k + M * k + M2 * k + M3 * k] * C[i] * sqrtOf;
        }

        // Rule 2: Creation on k s / Annihilation on k s
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                sqrtOf = v[k] * v[s];
                z += 4 * Hint[k + s*M + k*M2 + s*M3] * sqrtOf * C[i];
                /*
                z += Hint[s + k*M + k*M2 + s*M3] * sqrtOf * C[i];
                z += Hint[s + k*M + s*M2 + k*M3] * sqrtOf * C[i];
                z += Hint[k + s*M + s*M2 + k*M3] * sqrtOf * C[i];
                */
            }
        }

        // Rule 3: Creation on k k / Annihilation on q q
        // EXCLUDED since it only conserves the total
        // angular momentum if k = q that is rule 1

        // Rule 4: Creation on k k / Annihilation on k l
        // EXCLUDED since it only conserves the total
        // angular momentum if l = k that is rule 1


        // Rule 5: Creation on k s / Annihilation on s s
        // EXCLUDED since it only conserves the total
        // angular momentum if s = k that is rule 1

        // Rule 6: Creation on k k / Annihilation on q l
        // ONLY IN CASE q + l = 2 * k
        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < M; q++)
            {
                if (q == k || 2 * k - q < 0 || 2 * k - q >= M) continue;
                l = 2 * k - q;
                if (l > q) continue;

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

        // Rule 7: Creation on k s / Annihilation on q q
        // ONLY IN CASE k + s = 2 * q
        for (q = 0; q < M; q++)
        {
            for (k = 0; k < M; k++)
            {
                if (q == k || 2 * q - k < 0 || 2 * q - k >= M) continue;

                s = 2 * q - k;
                if (v[k] < 1 || v[s] < 1 || s > k) continue;

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

        // Rule 8: Creation on k s / Annihilation on s l
        // ONLY IN CASE k = l, BUT this case is included
        // in rule 2

        // Rule 9: Creation on k s / Annihilation on q l
        // ONLY IN CASE k + s = q + l
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = 0; q < M; q++)
                {
                    if (q == s || q == k) continue;
                    if (k + s - q < 0 || k + s - q >= M ) continue;
                    l = k + s - q;
                    if (l <= q) continue;

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

                }       // Finish q
            }           // Finish s
        }               // Finish k

        out[i] = w + 0.5 * z;
    }

    free(v);

}
