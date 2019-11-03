
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of São Paulo - Brazil

 Last update : November/02/2019

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
 * Compute the memory cost of the several data structures used  in
 * improvement of routines to compute the density matrices and the
 * Hamiltonian action on configuration  space.  Vary the number of
 * particles or the number of individual particle states and record
 * the results in 'mem_used.dat' file. The layout of the output
 * file has the following structure
 *
 * column 1 : Number of particles
 * column 2 : Number of orbitals
 * ------------------------------- Memory -------------------------------
 * column 3 : Coefficients storage
 * column 4 : Hashing table
 * column 5 : Mapping of single Jump from one orbital (1J1O)
 * column 6 : Mapping of double Jump from one orbital (2J1O)
 * column 7 : Mapping of double jump from two different orbitals (2J2O)
 *
 * ----------------------------------------------------------------------- */

#include "configurationsMap.h"
#include "outTextFile.h"



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

/** COLLECT TIME VARYING CONFIGURATIONAL PARAMETER, NUMBER OF  PARTICLES
  * OR NUMBER OF INDIVIDUAL PARTICLE STATES. TO SWITCH WICH PARAMETER IS
  * FIXED WITH THE ONE IS BEING VARIED, JUST SWAP THE PLACES  OF  'Npar'
  * AND 'Morb' VARIABLES IN THE NEXT TWO LINES                       **/

    Morb = 3;

    for (Npar = 900; Npar < 951; Npar += 50)
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

        printf("\n\nMemory for coefficients : %.3lf", Cmem);

        printf("\nMemory for Fock states : %.3lf", IFmem);

        printf("\nMemory single jump from one orbital map : %.3lf", MAPmem);

        printf("\nMemory double jump from one orbital map : %.3lf", MAPOTmem);

        printf("\nMemory double jump from two orbitals map : %.3lf", MAPTTmem);

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
