#include "onebodyMatrix.h"
#include "twobodyMatrix.h"
#include "outTextFile.h"
#include <time.h>



int main(int argc, char * argv[])
{

    int
        i,
        nc,
        Npar,
        Morb;

    double
        sum,
        time_used;

    clock_t
        start,
        end;

    Iarray
        Map,
        MapOT,
        MapTT,
        IFmat,
        NCmat,
        strideOT,
        strideTT;

    Carray
        C,
        rho2,
        rho2_X,
        rho2_XM;

    Cmatrix
        rho1_XM,
        rho1_X,
        rho1;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    nc = NC(Npar,Morb);

    NCmat = MountNCmat(Npar,Morb);
    IFmat = MountFocks(Npar,Morb);

    printf("\nNumber of particles : %3d", Npar);
    printf("\nNumber of orbitals  : %3d", Morb);
    printf("\nNumber of configurations : %d", nc);

    printf("\n\n======================================\n\n");

    printf("MEMORY CONSUMPTION (in Mb)");

    printf("\n\nMemory for coefficients : %.1lf",
            ((double) nc*sizeof(double complex)) / 1E6);

    printf("\nMemory for Fock states : %.1lf",
            ((double) nc*Morb*sizeof(int)) / 1E6);

    printf("\nMemory for one to one Map : %.1lf",
            ((double) nc*Morb*Morb*sizeof(int)) / 1E6);



    rho1 = (double complex **) malloc(Morb * sizeof(double complex *));
    for (i = 0; i < Morb; i++)
    {
        rho1[i] = (double complex *) malloc(Morb * sizeof(double complex));
    }

    rho1_X = (double complex **) malloc(Morb * sizeof(double complex *));
    for (i = 0; i < Morb; i++)
    {
        rho1_X[i] = (double complex *) malloc(Morb * sizeof(double complex));
    }

    rho1_XM = (double complex **) malloc(Morb * sizeof(double complex *));
    for (i = 0; i < Morb; i++)
    {
        rho1_XM[i] = (double complex *) malloc(Morb * sizeof(double complex));
    }



    C = carrDef(nc);
    rho2 = carrDef(Morb * Morb * Morb * Morb);
    rho2_X = carrDef(Morb * Morb * Morb * Morb);
    rho2_XM = carrDef(Morb * Morb * Morb * Morb);

    sum = 0.0;
    for (i = 0; i < nc; i++)
    {
        C[i] = sin( 20 * ((double) i) / nc) * (i % 13) - (i % 8) * I;
        sum = sum + creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
    }

    // normalize to 1
    for (i = 0; i < nc; i++) C[i] = C[i] / sqrt(sum);



    strideTT = iarrDef(nc);
    strideOT = iarrDef(nc);
    Map = OneOneMap(Npar,Morb,NCmat,IFmat);
    MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
    MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

    printf("\nMemory for one-two Map : %.1lf",
            ((double) strideOT[nc-1]*sizeof(int))/1E6);
    printf("\nMemory for two-two Map : %.1lf",
            ((double) strideTT[nc-1]*sizeof(int))/1E6);



    printf("\n\n======================================\n\n");

    printf("TIME DEMANDED");
/*
    start = clock();
    for (i = 0; i < 50; i++) OBrho(Npar,Morb,NCmat,C,rho1);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 50;
    printf("\n\nTime to setup rho1 : %.3lfms", time_used * 1000);
*/
    start = clock();
    for (i = 0; i < 50; i++) OBrho_X(Npar,Morb,NCmat,IFmat,C,rho1_X);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 50;
    printf("\n\nTime to setup rho1 with Fock states");
    printf(" : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 50; i++) OBrho_XM(Npar,Morb,Map,NCmat,IFmat,C,rho1_XM);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 50;
    printf("\n\nTime to setup rho1 using Maps and Fock states");
    printf(" : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++) TBrho(Npar,Morb,NCmat,IFmat,C,rho2);
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to setup rho2 : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++)
    {
        TBrho_X(Npar,Morb,Map,MapOT,strideOT,NCmat,IFmat,C,rho2_X);
    }
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to setup rho2 with one-map : %.3lfms", time_used * 1000);

    start = clock();
    for (i = 0; i < 5; i++)
    {
        TBrho_XX(Npar,Morb,Map,MapOT,MapTT,strideOT,strideTT,
                NCmat,IFmat,C,rho2_XM);
    }
    end = clock();
    time_used = ((double) (end - start)) / CLOCKS_PER_SEC / 5;
    printf("\n\nTime to setup rho2 with two-map : %.3lfms", time_used * 1000);

    cmat_txt("rho1_X.dat",Morb,Morb,rho1_X);
    cmat_txt("rho1_XM.dat",Morb,Morb,rho1_XM);

    carr_txt("rho2.dat",Morb*Morb*Morb*Morb,rho2);
    carr_txt("rho2_X.dat",Morb*Morb*Morb*Morb,rho2_X);
    carr_txt("rho2_XM.dat",Morb*Morb*Morb*Morb,rho2_XM);



    free(C);
    free(Map);

    for(i = 0; i < Morb; i++) free(rho1[i]);
    free(rho1);

    for(i = 0; i < Morb; i++) free(rho1_X[i]);
    free(rho1_X);

    for(i = 0; i < Morb; i++) free(rho1_XM[i]);
    free(rho1_XM);

    free(IFmat);
    free(NCmat);

    free(strideOT);
    free(strideTT);
    free(MapOT);
    free(MapTT);
    free(rho2);
    free(rho2_X);

    printf("\n\nDone.\n\n");
    return 0;
}
