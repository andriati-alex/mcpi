
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of São Paulo - Brazil

 Last update : 08/13/2019

---------------------------------------------------------------------------

 ****  EXECUTABLE TO MEASURE MEMORY REQUIRED ****
 *
 * compilation :
 * -------------
 *
 * icc MemCost.c -lm -o exe (if available)
 * gcc MemCost.c -lm -o exe
 *
 * comments :
 * ----------
 *
 *  On executing this program it will compute the memory for all  required
 *  structures for each improvement done, sweeping the number of particles.
 *
 * ----------------------------------------------------------------------- */

#include "onebodyMatrix.h"
#include "twobodyMatrix.h"
#include "hamiltonianMatrix.h"
#include "outTextFile.h"
#include <time.h>



int main(int argc, char * argv[])
{

    int
        nc,
        Npar,
        Morb;

    FILE
        * output;

    double
        Cmem,
        IFmem,
        MAPmem,
        MAPOTmem,
        MAPTTmem;

    Iarray
        Map,
        MapOT,
        MapTT,
        IFmat,
        NCmat,
        strideOT,
        strideTT;

    output = fopen("mem_used.dat", "w");

    Morb = 3;

    for (Npar = 30; Npar < 1050; Npar += 50)
    {

        printf("\n\n");

        nc = NC(Npar,Morb);

        NCmat = setupNCmat(Npar,Morb);
        IFmat = setupFocks(Npar,Morb);

        printf("\nNumber of particles : %3d", Npar);
        printf("\nNumber of orbitals  : %3d", Morb);
        printf("\nNumber of configurations : %d", nc);

        strideTT = iarrDef(nc);
        strideOT = iarrDef(nc);
        Map = OneOneMap(Npar,Morb,NCmat,IFmat);
        MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
        MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

        Cmem = ((double) nc * sizeof(double complex) ) / 1E6;
        IFmem = ((double) nc * Morb * sizeof(int) ) / 1E6;
        MAPmem = ((double) nc * Morb * Morb * sizeof(int) ) / 1E6;
        MAPOTmem = ((double) (strideOT[nc-1] + Morb*Morb) * sizeof(int)) / 1E6;
        MAPTTmem = ((double) (strideTT[nc-1]) * sizeof(int)) / 1E6;

        printf("\n\n======================================\n\n");

        printf("MEMORY CONSUMPTION (in Mb)");

        printf("\n\nMemory for coefficients : %.1lf", Cmem);

        printf("\nMemory for Fock states : %.1lf", IFmem);

        printf("\nMemory for one to one Map : %.1lf", MAPmem);

        printf("\nMemory for one-two Map : %.1lf", MAPOTmem);

        printf("\nMemory for two-two Map : %.1lf", MAPTTmem);

        free(Map);

        free(IFmat);
        free(NCmat);

        free(strideOT);
        free(strideTT);
        free(MapOT);
        free(MapTT);

        fprintf(output,"%d %d ",Npar,Morb);
        fprintf(output,"%.3lf %.3lf %.3lf %.3lf %.3lf",Cmem,IFmem,MAPmem,
                MAPOTmem,MAPTTmem);

        fprintf(output,"\n");

    }

    fclose(output);

    printf("\n\nDone.\n\n");
    return 0;
}
