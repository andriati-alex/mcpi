#ifndef _DensityCorrelation_h
#define _DensityCorrelation_h

#include "LBoseBoseFockSpace.h"
#include "LBoseFermiFockSpace.h"



void setSingleParticleStates(int lmax, Carray orbs)
{

/** SET THE PLANE WAVES CORRESPONDING TO THE SINGLE-PARTICLE BASIS **/

    int
        i,
        j,
        l,
        stride;
    Rarray
        x;

    x = setGrid(); // grid points linearly spaced between -PI and PI
    for (i = 0; i < 2*lmax+1; i++)
    {
        l = i-lmax; // dimensionless angular momentum of state 'i'
        for (j = 0; j < NGRID_POINTS; j++)
        {
            stride = NGRID_POINTS*i;
            orbs[stride+j] = cexp(I*l*x[i])/sqrt(2*PI);
        }
    }
    free(x);
}



void set1RDM(int N, int lmax, int mcsize, Iarray * ht, Carray C, Carray rho)
{

/** COMPUTE SINGLE PARTICLE REDUCED DENSITY MATRIX - ONE COMPONENT BOSONS **/

    int
        m,
        i,
        l;
    double
        cmod2;
    double complex
        avg;
    Rarray
        xpos,
        avgOcc;

    avgOcc = rarrDef(2*lmax+1); // vector of average occupation in SPS

    for (m = 0; m < 2*lmax+1; m++)
    {
        // Compute the average occupation in the state 'm'
        avg = 0.0;
        for (i = 0; i < mcsize; i++)
        {
            cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
            avg += cmod2*ht[i][m];
        }
        avgOcc[m] = creal(avg);
    }

    xpos = setDoubleGrid();
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        avg = 0.0;
        for (m = 0; m < 2*lmax+1; m++)
        {
            l = m-lmax; // dimensionless angular momentum of state 'm'
            avg += cexp(I*l*xpos[i])*avgOcc[m];
        }
        rho[i] = avg / (2*PI) / N;
    }
    free(avgOcc);
    free(xpos);
}



void mutualProb(int N, int lmax, int nc, Iarray * ht, Carray C, Carray rho)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        n,
        m,
        q,
        i,
        j,
        l,
        max_exchg;
    double
        bosef,
        cmod2;
    double complex
        sum;
    Iarray
        occ;
    Rarray
        xpos;
    Carray
        rhoL;

    max_exchg = 2*lmax; // MAXIMUM POSITIVE MOMENTUM THE ATOMS CAN EXCHANGE

    // vector of occupations numbers
    occ = iarrDef(2*lmax+1);
    // momentum exchanged from (-max_exchg) to (max_exchg)
    rhoL = carrDef(2*max_exchg+1);

    // TREAT FIRST THE ZERO MOMENTUM EXCHANGED CASE (q = 0)
    sum = 0.0;
    // remove particle of type 'A' with momentum 'n'
    for (n = 0; n < 2*lmax+1; n++)
    {
        // remove particle of type 'B' with momentum 'm'
        for (m = 0; m < 2*lmax+1; m++)
        {
            for (i = 0; i < nc; i++)
            {
                if (m != n) bosef = ht[i][n]*ht[i][m];
                else        bosef = ht[i][n]*(ht[i][n]-1);
                cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                sum = sum + bosef*cmod2;
            }
        }
    }
    rhoL[max_exchg] = sum;

    // CASES WITH MOMENTUM EXCHANGE (q != 0)
    for (q = -max_exchg; q <= max_exchg; q++)
    {
        sum = 0.0;
        if (q == 0) continue; // case computed above
        for (n = 0; n < 2*lmax+1; n++)
        {
            if (n+q < 0 || n+q > 2*lmax) continue;
            for (m = 0; m < 2*lmax+1; m++)
            {
                if (m-q < 0 || m-q > 2*lmax) continue;
                for (i = 0; i < nc; i++)
                {
                    for (l = 0; l < 2*lmax+1; l++) occ[l] = ht[i][l];
                    bosef = 1.0;
                    // destroy in state n+q
                    bosef *= occ[n+q];
                    occ[n+q] -= 1;
                    // destroy in state m-q
                    bosef *= occ[m-q];
                    occ[m-q] -= 1;
                    // create in state n
                    bosef *= (occ[n]+1);
                    occ[n] += 1;
                    // create in state m
                    bosef *= (occ[m]+1);
                    occ[m] += 1;
                    // check if some destruction operator
                    // acted on an empty state
                    if (fabs(bosef) < 1E-10) continue;
                    bosef = sqrt(bosef);
                    // find index of non-vanishing scalar product
                    j = BgetIndex(lmax,nc,ht,occ);
                    sum = sum + bosef*conj(C[i])*C[j];
                }
            }
        }
        rhoL[q+max_exchg] = sum;
    }

    xpos = setDoubleGrid(); // from -2*PI to 2*PI
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        sum = 0.0;
        for (m = 0; m < 2*max_exchg+1; m++)
        {
            l = m-max_exchg; // dimensionless ang. momentum exchanged
            sum += cexp(I*l*xpos[i])*rhoL[m];
        }
        rho[i] = sum/N/(N-1);
    }
    free(occ);
    free(xpos);
    free(rhoL);
}



void mixture_set1RDM(CompoundSpace S, Carray C, Carray rho, char type)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        m,
        i,
        l,
        N,
        lmax,
        mcsize;
    double
        cmod2;
    double complex
        avg;
    Iarray
        vA,
        vB;
    Rarray
        xpos,
        avgOcc;

    mcsize = S->size;
    // select the species below
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
    avgOcc = rarrDef(2*lmax+1); // vector of average occupation in SPS
    for (m = 0; m < 2*lmax+1; m++)
    {
        // Compute the average occupation in the state 'm'
        avg = 0.0;
        for (i = 0; i < mcsize; i++)
        {
            // get occupation numbers in vectors vA and vB
            BBgetConfigs(i,S,vA,vB);
            cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
            if (type == 'A' || type == 'a') avg += cmod2*vA[m];
            else                            avg += cmod2*vB[m];
        }
        avgOcc[m] = creal(avg);
    }

    xpos = setDoubleGrid(); // from -2*PI to 2*PI
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        avg = 0.0;
        for (m = 0; m < 2*lmax+1; m++)
        {
            l = m-lmax; // dimensionless angular momentum of state 'i'
            avg += cexp(I*l*xpos[i])*avgOcc[m];
        }
        rho[i] = avg / (2*PI) / N;
    }
    free(vA);
    free(vB);
    free(xpos);
    free(avgOcc);
}



void bosefermi_set1RDM(BFCompoundSpace S, Carray C, Carray rho, char type)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - MIXTURE OF BOSONS AND FERMIONS
    RETURN THE COMPONENT SELECTED BY 'type' VARIABLE                    **/

    int
        m,
        i,
        l,
        N,
        lmax,
        mcsize;
    double
        cmod2;
    double complex
        avg;
    Iarray
        vB;
    Farray
        vF;
    Rarray
        xpos,
        avgOcc;

    mcsize = S->size;
    // select the species below
    if (type == 'B' || type == 'b')
    {
        lmax = S->lmaxB;
        N = S->Nb;
    }
    else
    {
        lmax = S->lmaxF;
        N = S->Nf;
    }

    vB = iarrDef(2*S->lmaxB+1);
    vF = farrDef(2*S->lmaxF+1);
    avgOcc = rarrDef(2*lmax+1); // vector of average occupation in SPS
    for (m = 0; m < 2*lmax+1; m++)
    {
        // Compute the average occupation in the state 'm'
        avg = 0.0;
        for (i = 0; i < mcsize; i++)
        {
            // copy current configuration occupations
            BFgetConfigs(i,S,vB,vF);
            cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
            if (type == 'B' || type == 'b') avg += cmod2*vB[m];
            else                            avg += cmod2*vF[m];
        }
        avgOcc[m] = creal(avg);
    }

    xpos = setDoubleGrid(); // from -2*PI to 2*PI
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        avg = 0.0;
        for (m = 0; m < 2*lmax+1; m++)
        {
            l = m-lmax; // dimensionless angular momentum of state 'i'
            avg += cexp(I*l*xpos[i])*avgOcc[m];
        }
        rho[i] = avg / (2*PI) / N;
    }
    free(vB);
    free(vF);
    free(xpos);
    free(avgOcc);
}



void mixture_setCov(CompoundSpace S, Carray C, Carray rho)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        n,
        m,
        q,
        i,
        j,
        l,
        Na,
        Nb,
        Nla,
        Nlb,
        lmaxA,
        lmaxB,
        mcsize,
        max_exchg;
    double
        bosef,
        cmod2;
    double complex
        sum;
    Iarray
        vA,
        vB;
    Rarray
        xpos;
    Carray
        rhoL;

    // EXTRACT MULTI-CONFIG. SPACE DATA FROM STRUCTURE
    mcsize = S->size;
    lmaxA = S->lmaxA;   // maximum single particle angular momentum
    Na = S->Na;         // number of particles of species A
    Nla = 2*lmaxA+1;    // total number of single particle states
    lmaxB = S->lmaxB;
    Nb = S->Nb;
    Nlb = 2*lmaxB+1;

    vA = iarrDef(2*lmaxA+1);
    vB = iarrDef(2*lmaxB+1);

    // MAXIMUM POSITIVE MOMENTUM THE ATOMIC SPECIES CAN EXCHANGE
    // DUE TO INTER-SPECIES INTERACTIONS IS CONSTRAINED
    if (lmaxA > lmaxB) max_exchg = 2*lmaxB;
    else               max_exchg = 2*lmaxA;

    // momentum exchanged from (-max_exchg) to (max_exchg)
    rhoL = carrDef(2*max_exchg+1);

    // TREAT FIRST THE ZERO MOMENTUM EXCHANGED CASE (q = 0)
    sum = 0.0;
    // remove particle of type 'A' with momentum 'n'
    for (n = 0; n < Nla; n++)
    {
        // remove particle of type 'B' with momentum 'm'
        for (m = 0; m < Nlb; m++)
        {
            for (i = 0; i < mcsize; i++)
            {
                BBgetConfigs(i,S,vA,vB);
                bosef = vA[n]*vB[m]; // create/destroy in the same states
                cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                sum = sum + bosef * cmod2;
            }
        }
    }
    rhoL[max_exchg] = sum;

    for (q = -max_exchg; q <= max_exchg; q++)
    {
        sum = 0.0;
        if (q == 0) continue; // case computed above
        for (n = 0; n < Nla; n++)
        {
            if (n+q < 0 || n+q >= Nla) continue;
            for (m = 0; m < Nlb; m++)
            {
                if (m-q < 0 || m-q >= Nlb) continue;
                for (i = 0; i < mcsize; i++)
                {
                    BBgetConfigs(i,S,vA,vB);
                    if (vA[n+q]*vB[m-q] == 0) continue;
                    bosef = sqrt((double)vA[n+q]*vB[m-q]*(vA[n]+1)*(vB[m]+1));
                    // action of the operators to the left
                    vA[n+q] -= 1;
                    vA[n] += 1;
                    vB[m-q] -= 1;
                    vB[m] += 1;
                    // find index of non-vanishing scalar product
                    j = BBgetIndex(S,vA,vB);
                    sum = sum + bosef*conj(C[i])*C[j];
                    vA[n+q] += 1;
                    vA[n] -= 1;
                    vB[m-q] += 1;
                    vB[m] -= 1;
                }
            }
        }
        rhoL[q+max_exchg] = sum;
    }

    xpos = setDoubleGrid(); // from -2*PI to 2*PI
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        sum = 0.0;
        for (m = 0; m < 2*max_exchg+1; m++)
        {
            l = m-max_exchg; // dimensionless ang. momentum exchanged
            sum += cexp(I*l*xpos[i])*rhoL[m];
        }
        rho[i] = sum/(Na*Nb)-1.0;
    }
    free(vA);
    free(vB);
    free(xpos);
    free(rhoL);
}



void mixture_MutualProb(CompoundSpace S, Carray C, Carray rho, char type)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        n,
        m,
        q,
        i,
        j,
        l,
        Npar,
        lmax,
        mcsize,
        max_exchg;
    double
        bosef,
        cmod2;
    double complex
        sum;
    Iarray
        vA,
        vB;
    Rarray
        xpos;
    Carray
        rhoL;

    mcsize = S->size;

    vA = iarrDef(2*S->lmaxA+1);
    vB = iarrDef(2*S->lmaxB+1);

    // MAXIMUM POSITIVE MOMENTUM THE ATOMIC SPECIES CAN EXCHANGE
    // DUE TO INTER-SPECIES INTERACTIONS IS CONSTRAINED
    if (type == 'A' || type == 'a')
    {
        lmax = S->lmaxA;
        Npar = S->Na;
    }
    else
    {
        lmax = S->lmaxB;
        Npar = S->Nb;
    }
    max_exchg = 2*lmax; // modulus of maximum momentum exchange

    // momentum exchanged from (-max_exchg) to (max_exchg)
    rhoL = carrDef(2*max_exchg+1);

    // TREAT FIRST THE ZERO MOMENTUM EXCHANGED CASE (q = 0)
    sum = 0.0;
    // remove particle of type 'A' with momentum 'n'
    for (n = 0; n < 2*lmax+1; n++)
    {
        // remove particle of type 'B' with momentum 'm'
        for (m = 0; m < 2*lmax+1; m++)
        {
            for (i = 0; i < mcsize; i++)
            {
                BBgetConfigs(i,S,vA,vB);
                if (type == 'A' || type == 'a')
                {
                    if (m != n) bosef = vA[n]*vA[m];
                    else        bosef = vA[n]*(vA[n]-1);
                }
                else
                {
                    if (m != n) bosef = vB[n]*vB[m];
                    else        bosef = vB[n]*(vB[n]-1);
                }
                cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                sum = sum + bosef*cmod2;
            }
        }
    }
    rhoL[max_exchg] = sum;

    if (type == 'A' || type == 'a')
    {
        for (q = -max_exchg; q <= max_exchg; q++)
        {
            sum = 0.0;
            if (q == 0) continue; // case computed above
            for (n = 0; n < 2*lmax+1; n++)
            {
                if (n+q < 0 || n+q > 2*lmax) continue;
                for (m = 0; m < 2*lmax+1; m++)
                {
                    if (m-q < 0 || m-q > 2*lmax) continue;
                    for (i = 0; i < mcsize; i++)
                    {
                        BBgetConfigs(i,S,vA,vB);
                        bosef = 1.0;
                        // destroy in state n+q
                        bosef *= vA[n+q];
                        vA[n+q] -= 1;
                        // destroy in state m-q
                        bosef *= vA[m-q];
                        vA[m-q] -= 1;
                        // create in state n
                        bosef *= (vA[n]+1);
                        vA[n] += 1;
                        // create in state m
                        bosef *= (vA[m]+1);
                        vA[m] += 1;
                        // check if some destruction operator
                        // acted on an empty state
                        if (fabs(bosef) < 1E-10) continue;
                        bosef = sqrt(bosef);
                        // find index of non-vanishing scalar product
                        j = BBgetIndex(S,vA,vB);
                        sum = sum + bosef*conj(C[i])*C[j];
                    }
                }
            }
            rhoL[q+max_exchg] = sum;
        }
    }
    else
    {
        for (q = -max_exchg; q <= max_exchg; q++)
        {
            sum = 0.0;
            if (q == 0) continue; // case computed above
            for (n = 0; n < 2*lmax+1; n++)
            {
                if (n+q < 0 || n+q > 2*lmax) continue;
                for (m = 0; m < 2*lmax+1; m++)
                {
                    if (m-q < 0 || m-q > 2*lmax) continue;
                    for (i = 0; i < mcsize; i++)
                    {
                        BBgetConfigs(i,S,vA,vB);
                        bosef = 1.0;
                        // destroy in state n+q
                        bosef *= vB[n+q];
                        vB[n+q] -= 1;
                        // destoy in state m-q
                        bosef *= vB[m-q];
                        vB[m-q] -= 1;
                        // create in state n
                        bosef *= (vB[n]+1);
                        vB[n] += 1;
                        // create in state m
                        bosef *= (vB[m]+1);
                        vB[m] += 1;
                        if (fabs(bosef) < 1E-10) continue;
                        bosef = sqrt(bosef);
                        // find index of non-vanishing scalar product
                        j = BBgetIndex(S,vA,vB);
                        sum = sum + bosef*conj(C[i])*C[j];
                    }
                }
            }
            rhoL[q+max_exchg] = sum;
        }
    }

    xpos = setDoubleGrid(); // from -2*PI to 2*PI
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        sum = 0.0;
        for (m = 0; m < 2*max_exchg+1; m++)
        {
            l = m-max_exchg; // dimensionless ang. momentum exchanged
            sum += cexp(I*l*xpos[i])*rhoL[m];
        }
        rho[i] = sum/Npar/(Npar-1);
    }
    free(vA);
    free(vB);
    free(xpos);
    free(rhoL);
}



void bosefermi_setCov(BFCompoundSpace S, Carray C, Carray rho)
{

/** SINGLE PARTICLE REDUCED DENSITY MATRIX - SELECTION OF ONE 'type'
    PRESENT IN THE TWO-COMPONENT BOSONIC MIXTURE                 **/

    int
        n,
        m,
        q,
        i,
        j,
        l,
        Nb,
        Nf,
        Nlb,
        Nlf,
        lmaxB,
        lmaxF,
        mcsize,
        max_exchg;
    double
        fermif,
        bosef,
        cmod2;
    double complex
        sum;
    Iarray
        vB;
    Farray
        vF;
    Rarray
        xpos;
    Carray
        rhoL;

    // EXTRACT MULTI-CONFIG. SPACE DATA FROM STRUCTURE
    mcsize = S->size;
    lmaxB = S->lmaxB;   // maximum single particle angular momentum of bosons
    Nb = S->Nb;         // number of bosons
    Nlb = 2*lmaxB+1;    // total number of single particle states for bosons
    lmaxF = S->lmaxF;
    Nf = S->Nf;
    Nlf = 2*lmaxF+1;

    vB = iarrDef(2*lmaxB+1);
    vF = farrDef(2*lmaxF+1);

    // MAXIMUM POSITIVE MOMENTUM THE ATOMIC SPECIES CAN EXCHANGE
    // DUE TO INTER-SPECIES INTERACTIONS IS CONSTRAINED
    if (lmaxB > lmaxF) max_exchg = 2*lmaxF;
    else               max_exchg = 2*lmaxB;

    // momentum exchanged from (-max_exchg) to (max_exchg)
    rhoL = carrDef(2*max_exchg+1);

    // TREAT FIRST THE ZERO MOMENTUM EXCHANGED CASE (q = 0)
    sum = 0.0;
    // remove particle of type 'A' with momentum 'n'
    for (n = 0; n < Nlf; n++)
    {
        // remove particle of type 'B' with momentum 'm'
        for (m = 0; m < Nlb; m++)
        {
            for (i = 0; i < mcsize; i++)
            {
                BFgetConfigs(i,S,vB,vF);
                bosef = vF[n]*vB[m]; // create/destroy in the same states
                cmod2 = creal(C[i])*creal(C[i])+cimag(C[i])*cimag(C[i]);
                sum = sum + bosef * cmod2;
            }
        }
    }
    rhoL[max_exchg] = sum;

    for (q = -max_exchg; q <= max_exchg; q++)
    {
        sum = 0.0;
        if (q == 0) continue; // case computed above
        for (n = 0; n < Nlf; n++)
        {
            if (n+q < 0 || n+q >= Nlf) continue;
            for (m = 0; m < Nlb; m++)
            {
                if (m-q < 0 || m-q >= Nlb) continue;
                for (i = 0; i < mcsize; i++)
                {
                    BFgetConfigs(i,S,vB,vF);
                    if (vF[n+q]*vB[m-q] == 0 || vF[n] > 0) continue;
                    // action of bosonic operators to the left
                    bosef = sqrt((double)vB[m-q]*(vB[m]+1));
                    vB[m-q] -= 1;
                    vB[m]   += 1;
                    // fermi factor is taken into account moving the
                    // operators and inverting sign at each crossing
                    // among them until find the right place
                    fermif = 1;
                    for (l = 0; l < n+q; l++)
                    {
                        if (vF[l] == 1) fermif = (-1)*fermif;
                    }
                    vF[n+q] = 0;
                    for (l = 0; l < n; l++)
                    {
                        if (vF[l] == 1) fermif = (-1)*fermif;
                    }
                    vF[n] = 1;
                    // find index of non-vanishing scalar product
                    j = BFgetIndex(S,vB,vF);
                    sum = sum + bosef*conj(C[i])*C[j];
                    vF[n+q] = 1;
                    vF[n]   = 0;
                    vB[m-q] += 1;
                    vB[m]   -= 1;
                }
            }
        }
        rhoL[q+max_exchg] = sum;
    }

    xpos = setDoubleGrid(); // from -2*PI to 2*PI
    for (i = 0; i < 2*NGRID_POINTS+1; i++)
    {
        sum = 0.0;
        for (m = 0; m < 2*max_exchg+1; m++)
        {
            l = m-max_exchg; // dimensionless ang. momentum exchanged
            sum += cexp(I*l*xpos[i])*rhoL[m];
        }
        rho[i] = sum/(Nf*Nb)-1.0;
    }
    free(vF);
    free(vB);
    free(xpos);
    free(rhoL);
}



void SCANNING_ANALYSIS(char prefix [])
{
    int
        i,
        j,
        L,
        nc,
        Npar,
        lmax,
        Ncases;
    char
        fname[100],
        strnum[5];
    Iarray
        * ht;
    Rarray
        gridPoints;
    Carray
        C,
        RDM,
        mutprob;
    FILE
        * setup_file;

    // record grid points
    gridPoints = setDoubleGrid();
    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_grid.dat");
    rarr_txt(fname,2*NGRID_POINTS+1,gridPoints);
    free(gridPoints);

    RDM = carrDef(2*NGRID_POINTS+1);
    mutprob = carrDef(2*NGRID_POINTS+1);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_setup.dat");

    Ncases = NumberOfLines(fname)-1;
    setup_file = openFileRead(fname);
    ReachNewLine(setup_file);   // jump comment line started with #

    for (i = 0; i < Ncases; i++)
    {
        sprintf(strnum,"%d",i+1);
        printf("\n[%3d/%d] Working ...",i+1,Ncases);

        // READ PARAMETERS
        fscanf(setup_file,"%d",&Npar);  // Num. of particles type A
        fscanf(setup_file,"%d",&lmax);  // max. individual ang. momemtum
        fscanf(setup_file,"%d",&L);     // Total ang. momentum constraint
        ReachNewLine(setup_file);

        nc = BFixedMom_mcsize(Npar,lmax,L);
        ht = BAssembleHT(Npar,lmax,L,nc);
        C = carrDef(nc);
        // set file name with input data
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_job");
        strcat(fname,strnum);
        strcat(fname,".dat");
        // read the input data
        carr_input_txt(fname,nc,C);
        assert_norm(nc,C);

        set1RDM(Npar,lmax,nc,ht,C,RDM);
        mutualProb(Npar,lmax,nc,ht,C,mutprob);

        // set file name with 1-RDM
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_rdm");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,RDM);
        // set file name with mutual-probability
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_mut");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,mutprob);

        free(C);
        for (j = 0; j < nc; j++) free(ht[j]);
        free(ht);
    }

    free(RDM);
    free(mutprob);
    fclose(setup_file);

}



void SCANNING_MIXTURE_ANALYSIS(char prefix [])
{
    int
        i,
        L,
        NparA,
        NparB,
        lmaxA,
        lmaxB,
        Ncases;
    char
        fname[100],
        strnum[5];
    Rarray
        gridPoints;
    Carray
        C,
        RDM_A,
        RDM_B,
        denCov,
        mutprob_A,
        mutprob_B;
    CompoundSpace
        MixSpace;
    FILE
        * setup_file;

    // record grid points
    gridPoints = setDoubleGrid();
    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_grid.dat");
    rarr_txt(fname,2*NGRID_POINTS+1,gridPoints);
    free(gridPoints);

    RDM_A = carrDef(2*NGRID_POINTS+1);
    RDM_B = carrDef(2*NGRID_POINTS+1);
    denCov = carrDef(2*NGRID_POINTS+1);
    mutprob_A = carrDef(2*NGRID_POINTS+1);
    mutprob_B = carrDef(2*NGRID_POINTS+1);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_setup.dat");

    Ncases = NumberOfLines(fname)-1;
    setup_file = openFileRead(fname);
    ReachNewLine(setup_file);

    for (i = 0; i < Ncases; i++)
    {
        sprintf(strnum,"%d",i+1);
        printf("\n[%3d/%d] Working ...",i+1,Ncases);

        // READ PARAMETERS
        fscanf(setup_file,"%d",&NparA);  // Num. of particles type A
        fscanf(setup_file,"%d",&lmaxA);  // max. individual ang. momemtum
        fscanf(setup_file,"%d",&NparB);  // Num. of particles type B
        fscanf(setup_file,"%d",&lmaxB);  // max. individual ang. momentum
        fscanf(setup_file,"%d",&L);      // Total ang. momentum constraint
        ReachNewLine(setup_file);

        MixSpace = AllocCompBasis(NparA,NparB,lmaxA,lmaxB,L);
        PrintMixSpaceInfo(MixSpace);
        C = carrDef(MixSpace->size);
        // set file name with input data
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_job");
        strcat(fname,strnum);
        strcat(fname,".dat");
        // read the input data
        carr_input_txt(fname,MixSpace->size,C);
        assert_norm(MixSpace->size,C);

        mixture_setCov(MixSpace,C,denCov);
        mixture_set1RDM(MixSpace,C,RDM_A,'A');
        mixture_set1RDM(MixSpace,C,RDM_B,'B');
        mixture_MutualProb(MixSpace,C,mutprob_A,'A');
        mixture_MutualProb(MixSpace,C,mutprob_B,'B');

        // set file name with covariance function in grid points
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_cov");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,denCov);
        // set file name with 1-RDM of species A
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_rdmA");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,RDM_A);
        // set file name with 1-RDM of species B
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_rdmB");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,RDM_B);
        // set file name with mutual-probability of species A
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_mutA");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,mutprob_A);
        // set file name with mutual-probability of species B
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_mutB");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,mutprob_B);

        free(C);
        freeCompSpace(MixSpace);
    }

    free(RDM_A);
    free(RDM_B);
    free(denCov);
    free(mutprob_A);
    free(mutprob_B);
    fclose(setup_file);

}



void SCANNING_BOSEFERMI_ANALYSIS(char prefix [])
{
    int
        i,
        L,
        NparF,
        NparB,
        lmaxF,
        lmaxB,
        Ncases;
    double
        g[3];
    char
        fname[100],
        strnum[5];
    Rarray
        gridPoints;
    Carray
        C,
        RDM_B,
        RDM_F,
        denCov;
    BFCompoundSpace
        MixSpace;
    FILE
        * setup_file;

    // record grid points
    gridPoints = setDoubleGrid();
    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_grid.dat");
    rarr_txt(fname,2*NGRID_POINTS+1,gridPoints);
    free(gridPoints);

    RDM_F = carrDef(2*NGRID_POINTS+1);
    RDM_B = carrDef(2*NGRID_POINTS+1);
    denCov = carrDef(2*NGRID_POINTS+1);

    strcpy(fname,OUTPUT_PATH);
    strcat(fname,prefix);
    strcat(fname,"_setup.dat");

    Ncases = NumberOfLines(fname)-1;
    setup_file = openFileRead(fname);
    ReachNewLine(setup_file);

    for (i = 0; i < Ncases; i++)
    {
        sprintf(strnum,"%d",i+1);
        printf("\n[%3d/%d] Working ...",i+1,Ncases);

        // READ PARAMETERS
        fscanf(setup_file,"%d",&NparB);  // Num. of particles type A
        fscanf(setup_file,"%d",&lmaxB);  // max. individual ang. momemtum
        fscanf(setup_file,"%d",&NparF);  // Num. of particles type B
        fscanf(setup_file,"%d",&lmaxF);  // max. individual ang. momentum
        fscanf(setup_file,"%d",&L);      // Total ang. momentum constraint
        ReachNewLine(setup_file);

        MixSpace = BoseFermiBasis(NparB,NparF,lmaxB,lmaxF,L);
        PrintBFSpaceInfo(MixSpace);
        C = carrDef(MixSpace->size);
        // set file name with input data
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_job");
        strcat(fname,strnum);
        strcat(fname,".dat");
        // read the input data
        carr_input_txt(fname,MixSpace->size,C);
        assert_norm(MixSpace->size,C);

        bosefermi_setCov(MixSpace,C,denCov);
        bosefermi_set1RDM(MixSpace,C,RDM_B,'B');
        bosefermi_set1RDM(MixSpace,C,RDM_F,'F');

        // set file name with covariance function in grid points
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_cov");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,denCov);
        // set file name with 1-RDM of species A
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_rdmB");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,RDM_B);
        // set file name with 1-RDM of species B
        strcpy(fname,OUTPUT_PATH);
        strcat(fname,prefix);
        strcat(fname,"_rdmF");
        strcat(fname,strnum);
        strcat(fname,".dat");
        carr_txt(fname,2*NGRID_POINTS+1,RDM_F);

        free(C);
        freeBoseFermiSpace(MixSpace);
    }

    free(RDM_F);
    free(RDM_B);
    free(denCov);
    fclose(setup_file);

}



#endif
