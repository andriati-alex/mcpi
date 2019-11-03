#ifndef _onebodyMatrix_h
#define _onebodyMatrix_h

#include <math.h>
#include "configurationsMap.h"


/****   AUTHOR INFORMATION

 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : Oct/31/2019

****/


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



void OBrho(int N, int M, Iarray NCmat, Carray C, Cmatrix rho)
{

/** The most naive implementation, without even the Hashing table
  * of Fock states **/

    int i,
        j,
        k,
        l,
        nc;

    Iarray
        v;

    double
        mod2;

    double complex
        RHO;

    nc = NC(N,M);
    v = iarrDef(M);

    for (k = 0; k < M; k++)
    {

        // Diagonal elements

        RHO = 0;

        for (i = 0; i < nc; i++)
        {
            IndexToFock(i,N,M,v);
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO = RHO + mod2 * v[k];
        }

        rho[k][k] = RHO;



        // Off-diagonal elements. Take advantage of hermiticity
        // of the one-body density matrix

        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            for (i = 0; i < nc; i++)
            {
                IndexToFock(i,N,M,v);
                if (v[k] < 1) continue;

                v[k] -= 1;
                v[l] += 1;

                // convert to index
                j = FockToIndex(N, M, NCmat, v);

                v[k] += 1;
                v[l] -= 1;

                RHO += conj(C[i]) * C[j] * sqrt((double)(v[l]+1)*v[k]);
            }

            rho[k][l] = RHO;
            // Use hermiticity
            rho[l][k] = conj(RHO);
        }

    }

    free(v);

}





void OBrho_X(int N, int M, Iarray NCmat, Iarray IF, Carray C, Cmatrix rho)
{

/** First basic improvement given by a predefined Hashing table in IF
  * array, that means Index to Fock, see setupFocks routine in the
  * configurationsMap.h file **/

    int i,
        j,
        k,
        l,
        vk,
        vl,
        nc;

    double
        mod2;

    double complex
        RHO;

    nc = NC(N,M);

    for (k = 0; k < M; k++)
    {

        // Diagonal elements

        RHO = 0;

        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO = RHO + mod2 * IF[k + i*M];
        }

        rho[k][k] = RHO;



        // Off-diagonal elements

        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            for (i = 0; i < nc; i++)
            {
                vk = IF[k + i*M];
                if (vk < 1) continue;
                vl = IF[l + i*M];

                IF[k + i*M] -= 1;
                IF[l + i*M] += 1;
                j = FockToIndex(N,M,NCmat,&IF[M*i]);
                IF[k + i*M] += 1;
                IF[l + i*M] -= 1;

                RHO += conj(C[i]) * C[j] * sqrt((double)(vl+1)*vk);
            }

            rho[k][l] = RHO;
            rho[l][k] = conj(RHO);
        }

    }

}





void OBrho_XM(int N, int M, Iarray Map, Iarray NCmat, Iarray IF,
     Carray C, Cmatrix rho)
{

/** Second improvement given by a Mapping structure between  configurations
  * related by jump of a particle from one orbital to another, given in Map
  * array. See OneOneMap routine in configurationsMap.h **/

    int i,
        j,
        k,
        l,
        vk,
        vl,
        nc;

    double
        mod2;

    double complex
        RHO;

    nc = NC(N,M);

    for (k = 0; k < M; k++)
    {

        // Diagonal elements

        RHO = 0;

        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO = RHO + mod2 * IF[k + i*M];
        }

        rho[k][k] = RHO;

        // Off-diagonal elements

        for (l = k + 1; l < M; l++)
        {

            RHO = 0;

            for (i = 0; i < nc; i++)
            {
                vk = IF[k + M*i];
                if (vk < 1) continue;
                vl = IF[l + M*i];

                j = Map[i + k * nc + l * M * nc];
                RHO += conj(C[i]) * C[j] * sqrt((double)(vl+1) * vk);
            }

            rho[k][l] = RHO;
            rho[l][k] = conj(RHO);
        }
    }

}



#endif
