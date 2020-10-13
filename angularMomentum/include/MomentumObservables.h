#ifndef _MomentumObservables_h
#define _MomentumObservables_h

#include "LBoseBoseFockSpace.h"
#include "LBoseFermiFockSpace.h"



void single_avgocc(int N, int lmax, int nc, Iarray * ht, Carray C, Rarray occ)
{

/** COMPUTE SINGLE PARTICLE REDUCED DENSITY MATRIX - ONE COMPONENT BOSONS **/

    int
        m,
        i;
    double
        avgocc,
        cmod2;

    for (m = 0; m < 2*lmax+1; m++)
    {
        // Compute the average occupation in the state 'm'
        avgocc = 0.0;
        for (i = 0; i < nc; i++)
        {
            cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
            avgocc += cmod2*ht[i][m];
        }
        occ[m] = avgocc/N;
    }
}



void mixture_avgocc(CompoundSpace S, Carray C, char type, Rarray occ)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        m,
        i,
        N,
        lmax;
    double
        avgOcc,
        cmod2;
    Iarray
        vA,
        vB;

    // select one of the species in the system
    if (type == 'A' || type == 'a')
    {
        lmax = S->lmaxA;
        N = S->Na;
    }
    else
    {
        lmax = S->lmaxB;
        N = S->Nb;
    }

    vA = iarrDef(2*S->lmaxA+1);
    vB = iarrDef(2*S->lmaxB+1);

    for (m = 0; m < 2*lmax+1; m++)
    {
        // Compute the average occupation in the state 'm'
        avgOcc = 0.0;
        for (i = 0; i < S->size; i++)
        {
            // get occupation numbers in vectors vA and vB
            BBgetConfigs(i,S,vA,vB);
            cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
            if (type == 'A' || type == 'a') avgOcc += cmod2*vA[m];
            else                            avgOcc += cmod2*vB[m];
        }
        // compute contribution for angular momentum
        occ[m] = avgOcc/N;
    }

    free(vA);
    free(vB);
}



double mixture_avgmom(CompoundSpace S, Carray C, char type)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        m,
        i,
        N,
        lmax;
    double
        avgAngMom,
        avgOcc,
        cmod2;
    Iarray
        vA,
        vB;

    // select one of the species in the system
    if (type == 'A' || type == 'a')
    {
        lmax = S->lmaxA;
        N = S->Na;
    }
    else
    {
        lmax = S->lmaxB;
        N = S->Nb;
    }

    vA = iarrDef(2*S->lmaxA+1);
    vB = iarrDef(2*S->lmaxB+1);

    avgAngMom = 0;
    for (m = 0; m < 2*lmax+1; m++)
    {
        // Compute the average occupation in the state 'm'
        avgOcc = 0.0;
        for (i = 0; i < S->size; i++)
        {
            // get occupation numbers in vectors vA and vB
            BBgetConfigs(i,S,vA,vB);
            cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
            if (type == 'A' || type == 'a') avgOcc += cmod2*vA[m];
            else                            avgOcc += cmod2*vB[m];
        }
        // compute contribution for angular momentum
        avgAngMom += (m-lmax)*avgOcc;
    }

    free(vA);
    free(vB);

    return avgAngMom/N;
}



double mixture_momcov(CompoundSpace S, Carray C)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        n,
        m,
        i,
        Na,
        Nb,
        Nla,
        Nlb,
        lmaxA,
        lmaxB;
    double
        sum,
        cmod2,
        twobody_avg;
    Iarray
        vA,
        vB;

    // EXTRACT MULTI-CONFIG. SPACE DATA FROM STRUCTURE
    lmaxA = S->lmaxA;   // maximum single particle angular momentum
    Na = S->Na;         // number of particles of species A
    Nla = 2*lmaxA+1;    // total number of single particle states
    lmaxB = S->lmaxB;
    Nb = S->Nb;
    Nlb = 2*lmaxB+1;

    vA = iarrDef(2*lmaxA+1);
    vB = iarrDef(2*lmaxB+1);

    twobody_avg = 0;
    // particle number operator of type 'A' with momentum 'n'
    for (n = 0; n < Nla; n++)
    {
        // particle number operator of type 'B' with momentum 'm'
        for (m = 0; m < Nlb; m++)
        {
            sum = 0;
            for (i = 0; i < S->size; i++)
            {
                BBgetConfigs(i,S,vA,vB);
                cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                sum = sum + vA[n]*vB[m]*cmod2;
            }
            twobody_avg += sum*(n-lmaxA)*(m-lmaxB);
        }
    }

    free(vA);
    free(vB);

    return twobody_avg/Na/Nb - mixture_avgmom(S,C,'A')*mixture_avgmom(S,C,'B');
}



double mixture_momvariance(CompoundSpace S, Carray C, char type)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        n,
        m,
        i,
        N,
        lmax;
    double
        sum,
        cmod2,
        twobody_avg;
    Iarray
        vA,
        vB;

    vA = iarrDef(2*S->lmaxA+1);
    vB = iarrDef(2*S->lmaxB+1);

    if (type == 'A' || type == 'a')
    {
        lmax = S->lmaxA;
        N = S->Na;
    }
    else
    {
        lmax = S->lmaxB;
        N = S->Nb;
    }

    if (type == 'A' || type == 'a')
    {
        twobody_avg = 0;
        // remove particle of type 'A' with momentum 'n'
        for (n = 0; n < 2*lmax+1; n++)
        {
            // remove particle of type 'B' with momentum 'm'
            for (m = 0; m < 2*lmax+1; m++)
            {
                sum = 0.0;
                for (i = 0; i < S->size; i++)
                {
                    BBgetConfigs(i,S,vA,vB);
                    cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                    sum = sum + vA[n]*vA[m]*cmod2;
                }
                twobody_avg += sum*(n-lmax)*(m-lmax);
            }
        }
    }
    else
    {
        twobody_avg = 0;
        // remove particle of type 'A' with momentum 'n'
        for (n = 0; n < 2*lmax+1; n++)
        {
            // remove particle of type 'B' with momentum 'm'
            for (m = 0; m < 2*lmax+1; m++)
            {
                sum = 0.0;
                for (i = 0; i < S->size; i++)
                {
                    BBgetConfigs(i,S,vA,vB);
                    cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                    sum = sum + vB[n]*vB[m]*cmod2;
                }
                twobody_avg += sum*(n-lmax)*(m-lmax);
            }
        }
    }

    free(vA);
    free(vB);

    return twobody_avg/N/N - mixture_avgmom(S,C,type)*mixture_avgmom(S,C,type);
}



#endif
