
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of São Paulo - Brazil

 Last update : 08/13/2019

---------------------------------------------------------------------------

 ****  PROGRAM TO DEMONSTRATE HOW THE MAPPING OF CONFIGURATIONS WORKS

 * Demonstrate the basics functions useful for computation of many-body
 * quantities in Fock state basis. Compile it :
 *
 * icc demonstrateFockMap.c -lm -o exe (if available)
 * gcc demonstrateFockMap.c -lm -o exe
 *
 * ./exe Nparticles Norbitals
 *
 * NOTE : this is suppose to be a simple demonstration, then do not use
 * large number of particles or orbitals because it would mess up the
 * output on the screen.
 *
 *------------------------------------------------------------------------ */



#include "configurationsMap.h"



int main(int argc, char * argv[])
{

    int
        i,
        j,
        k,
        s,
        q,
        l,
        g,
        h,
        nc,
        Npar, // number of particles
        Morb, // number of orbitals
        chunks,
        map2index;

    Iarray
        NCmat,
        IFmat,
        Map,
        MapTT,
        MapOT,
        stridesOT,
        stridesTT;

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

    printf("\n\nNumber of particles : %3d", Npar);
    printf(  "\nNumber of orbitals  : %3d", Morb);
    printf(  "\nNumber of configurations : %d", nc);

    if (nc > 5000 || Morb > 7)
    {
        printf("\n\nWARNING: lARGE SYSTEM.\n\n");
    }

    stridesTT = iarrDef(nc);
    stridesOT = iarrDef(nc);

    NCmat = setupNCmat(Npar,Morb);
    IFmat = setupFocks(Npar,Morb);

    Map = OneOneMap(Npar,Morb,NCmat,IFmat);
    MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,stridesTT);
    MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,stridesOT);

    printf("\nSize of one to one jump map : %d", nc * Morb * Morb);
    printf("\nSize of two to two jump map : %d", stridesTT[nc-1]);
    printf("\nSize of one to two jump map : %d", stridesOT[nc-1] + Morb*Morb);

    printf("\n\n\n");



    if (nc < 10000)
    {

        printf("Configuration  [");

        for (i = 0; i < Morb - 1; i++) printf(" Orb%d ,",i + 1);
        printf(" Orb%d ]", Morb);

        printf("    Transitions : from -> to\n");

        printf("=============================================");
        printf("================================");

        for (i = 0; i < nc; i++)
        {
            printf("\n%8d       [", i);

            for (j = 0; j < Morb - 1; j++) printf(" %3d  ,",IFmat[j+Morb*i]);
            printf(" %3d  ]", IFmat[Morb - 1 + Morb*i]);

            k = FockToIndex(Npar,Morb,NCmat,&IFmat[Morb*i]);
            if (k != i)
            {
                printf("\n\nERROR: Wrong map from FockToIndex\n\n");
            }

            k = i % Morb;
            j = (2 * i / 3 + 1) % Morb;
            printf("    %d -> %d", k + 1, j + 1);
            printf(" index = %4d",Map[i+k*nc+j*nc*Morb]);
        }

        printf("\n\n\nDouble transitions from different orbitals.\n\n");
        printf("=============================================");
        printf("================================\n\n");

        for (i = 0; i < nc; i++)
        {
            printf("\n%8d       [", i);

            for (j = 0; j < Morb - 1; j++) printf(" %3d  ,",IFmat[j+Morb*i]);
            printf(" %3d  ]", IFmat[Morb - 1 + Morb*i]);

            k = FockToIndex(Npar,Morb,NCmat,&IFmat[Morb*i]);
            if (k != i)
            {
                printf("\n\nERROR: Wrong map from FockToIndex\n\n");
            }

            k = i % (Morb - 1);
            s = (3 * (i + 1) / 2) % (Morb - k - 1) + k + 1;
            q = (11 * (i + 1) / 3 - 2) % Morb;
            l = (17 * i / 8) % Morb;

            if (IFmat[k+i*Morb] < 1 || IFmat[s+i*Morb] < 1) continue;

            chunks = 0;
            printf("    (%d,%d) -> (%d,%d)", k+1,s+1,q+1,l+1);

            for (h = 0; h < k; h++)
            {
                for (g = h + 1; g < Morb; g++)
                {
                    if (IFmat[h+i*Morb] > 0 && IFmat[g+i*Morb] > 0) chunks++;
                }
            }

            for (g = k + 1; g < s; g++)
            {
                if (IFmat[k+i*Morb] > 0 && IFmat[g+i*Morb] > 0) chunks++;
            }


            map2index = chunks*Morb*Morb + q + l * Morb + stridesTT[i];
            printf(" index = %4d",MapTT[map2index]);
        }



        printf("\n\n\nDouble transitions from the same orbital.\n\n");
        printf("=============================================");
        printf("================================\n\n");

        for (i = 0; i < nc; i++)
        {
            printf("\n%8d       [", i);

            for (j = 0; j < Morb - 1; j++) printf(" %3d  ,",IFmat[j+Morb*i]);
            printf(" %3d  ]", IFmat[Morb - 1 + Morb*i]);

            k = FockToIndex(Npar,Morb,NCmat,&IFmat[Morb*i]);
            if (k != i)
            {
                printf("\n\nERROR: Wrong map from FockToIndex\n\n");
            }

            k = i % Morb;
            q = (11 * (i + 1) / 3 - 2) % Morb;
            l = (17 * i / 8) % Morb;

            if (IFmat[k+i*Morb] < 2) continue;

            chunks = 0;
            printf("    (%d) -> (%d,%d)", k+1,q+1,l+1);

            for (h = 0; h < k; h++)
            {
                if (IFmat[h+i*Morb] > 1) chunks++;
            }

            map2index = chunks*Morb*Morb + q + l * Morb + stridesOT[i];
            printf(" index = %4d",MapOT[map2index]);
        }

    }

    else
    {
        printf("\n\nWARNING : Not printing, too large system with ");
        printf("total number of combinations higher than 10000\n");
    }



    free(IFmat);
    free(NCmat);
    free(Map);
    free(MapTT);
    free(stridesTT);
    free(MapOT);
    free(stridesOT);

    printf("\n\nDone.\n");
    return 0;
}
