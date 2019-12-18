
/****   AUTHOR INFORMATION
 
 NAME : Alex Valerio Andriati
 AFFILIATION : University of Sao Paulo - Brazil

 Last update : November/02/2019

 -------------------------------------------------------------------------

 ****  TEST ROUTINES TO COMPUTE PHYSICAL OPERATORS
 *
 * COMPILE :
 *
 * nvcc performanceTest.c -o exe
 *
 * HOW TO EXECUTE :
 *
 * ./exe Nparticles Norbitals
 *
 * where 'Nparticles' and 'Morbitals' are command line arguments for the
 * number of particles and individual particle states respectively.
 *
 * In order to measure time one may use the command nvproof before ./exe
 * For more information check out CUDA developer guide
 *
 * ----------------------------------------------------------------------- */

#include "hamiltonianMatrix.cuh"



int main(int argc, char * argv[])
{

    int
        i,
        j,
        q,
        l,
        nc,
        Npar,
        Morb,
        blocks;

    size_t
        nbytes;

    double
        sum,
        realPart,
        imagPart;

    cudaError_t
        err;

    cuDoubleComplex
        z;

    Iarray
        Map,
        d_Map,
        MapOT,
        d_MapOT,
        MapTT,
        d_MapTT,
        IFmat,
        d_IFmat,
        NCmat,
        strideOT,
        d_strideOT,
        strideTT,
        d_strideTT;

    Carray
        C,
        d_C,
        out,
        d_out,
        Ho,
        d_Ho,
        Hint,
        d_Hint;



    if (argc != 3)
    {
        printf("\n\nERROR: Need two integer numbers from command line ");
        printf("the first number of particles and second the number of ");
        printf("orbitals.\n\n");
        exit(EXIT_FAILURE);
    }

    err = cudaSuccess;

    sscanf(argv[1],"%d",&Npar);
    sscanf(argv[2],"%d",&Morb);
    nc = NC(Npar,Morb);

    NCmat = setupNCmat(Npar,Morb);
    IFmat = setupFocks(Npar,Morb);

    // alloc NCmat and IFmat on device

    cuda_iarrDef( nc * Morb , &d_IFmat );

    strideTT = iarrDef(nc);
    strideOT = iarrDef(nc);

    // alloc Map's strides on device

    cuda_iarrDef( nc , &d_strideTT );
    cuda_iarrDef( nc , &d_strideOT );

    Map = OneOneMap(Npar,Morb,NCmat,IFmat);
    MapTT = TwoTwoMap(Npar,Morb,NCmat,IFmat,strideTT);
    MapOT = OneTwoMap(Npar,Morb,NCmat,IFmat,strideOT);

    // alloc Maps on device

    cuda_iarrDef( nc * Morb * Morb , &d_Map );
    cuda_iarrDef( strideOT[nc-1] + Morb*Morb , &d_MapOT );
    cuda_iarrDef( strideTT[nc-1] , &d_MapTT );

    // Copy data from host to device

    nbytes = nc * Morb * sizeof(int);
    err = cudaMemcpy(d_IFmat,IFmat,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    nbytes = nc * sizeof(int);
    err = cudaMemcpy(d_strideOT,strideOT,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    nbytes = nc * sizeof(int);
    err = cudaMemcpy(d_strideTT,strideTT,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    nbytes = nc * Morb * Morb * sizeof(int);
    err = cudaMemcpy(d_Map,Map,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    nbytes = (strideOT[nc-1] + Morb*Morb) * sizeof(int);
    err = cudaMemcpy(d_MapOT,MapOT,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    nbytes = strideTT[nc-1] * sizeof(int);
    err = cudaMemcpy(d_MapTT,MapTT,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    printf("\nNumber of particles : %3d", Npar);
    printf("\nNumber of orbitals  : %3d", Morb);
    printf("\nNumber of configurations : %d", nc);

    printf("\n\n======================================\n\n");

    printf("MEMORY CONSUMPTION (in Mb)");

    printf("\n\nMemory for coefficients : %.1lf",
            ((double) 2*nc*sizeof(double)) / 1E6);

    printf("\nMemory for Fock states : %.1lf",
            ((double) nc*Morb*sizeof(int)) / 1E6);

    printf("\nMemory for single jump from 1 orbital Map : %.1lf",
            ((double) nc*Morb*Morb*sizeof(int)) / 1E6);

    printf("\nMemory for double jump from 1 orbital Map : %.1lf",
            ((double) strideOT[nc-1]*sizeof(int))/1E6);
    printf("\nMemory for double jump from 2 orbitals Map : %.1lf",
            ((double) strideTT[nc-1]*sizeof(int))/1E6);



    Ho = carrDef(Morb*Morb);

    Hint = carrDef(Morb*Morb*Morb*Morb);

    C = carrDef(nc);

    out = carrDef(nc);

    // alloc device arrays
    cuda_carrDef(nc,&d_out);
    cuda_carrDef(nc,&d_C);
    cuda_carrDef(Morb*Morb,&d_Ho);
    cuda_carrDef(Morb*Morb*Morb*Morb,&d_Hint);

    sum = 0.0;
    for (i = 0; i < nc; i++)
    {
        realPart = sin( 20 * ((double) i) / nc) * (i % 13);
        imagPart = (i % 8) - (i % 3);
        C[i] = make_cuDoubleComplex(realPart,imagPart);
        sum = sum + realPart * realPart + imagPart * imagPart;
    }

    // normalize to 1
    for (i = 0; i < nc; i++)
    {
        C[i] = cuCdiv(C[i],make_cuDoubleComplex(sqrt(sum),0));
    }

    // Done in the Host copy to device
    nbytes = nc * sizeof(cuDoubleComplex);
    err = cudaMemcpy(d_C,C,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }



    for (i = 0; i < Morb; i++)
    {
        Ho[i + Morb*i] = make_cuDoubleComplex((i % 4) - 1,0);
        for (j = i + 1; j < Morb; j++)
        {
            realPart = i * (j % 3) - (i % 4) + 5 * (j % 2);
            imagPart = -4.123 * i / (j + 1);
            Ho[i + Morb*j] = make_cuDoubleComplex(realPart,imagPart);
            Ho[j + Morb*i] = cuConj(Ho[i + Morb*j]);
        }
    }

    for (i = 0; i < Morb*Morb*Morb*Morb; i++)
    {
        Hint[i] = make_cuDoubleComplex(1.234,0);
    }

    for (i = 0; i < Morb; i++)
    {
        for (j = i + 1; j < Morb; j++)
        {
            for (q = 0; q < Morb; q++)
            {
                if (q == i || q == j) continue;
                for (l = q + 1; l < Morb; l++)
                {
                    if (l == i || l == j) continue;
                    // real part
                    realPart = i - 2 + 10 * (j % (i+1)) - q * l;
                    // imag part
                    imagPart = ((double) i * q - j * l) / Morb;

                    z = make_cuDoubleComplex(realPart,imagPart);
                    Hint[i+j*Morb+q*Morb*Morb+l*Morb*Morb*Morb] = z;
                    Hint[i+j*Morb+l*Morb*Morb+q*Morb*Morb*Morb] = z;
                    Hint[j+i*Morb+l*Morb*Morb+q*Morb*Morb*Morb] = z;
                    Hint[j+i*Morb+q*Morb*Morb+l*Morb*Morb*Morb] = z;
                    Hint[q+l*Morb+i*Morb*Morb+j*Morb*Morb*Morb] = cuConj(z);
                    Hint[q+l*Morb+j*Morb*Morb+i*Morb*Morb*Morb] = cuConj(z);
                    Hint[l+q*Morb+i*Morb*Morb+j*Morb*Morb*Morb] = cuConj(z);
                    Hint[l+q*Morb+j*Morb*Morb+i*Morb*Morb*Morb] = cuConj(z);
                }
            }
        }
    }

    // Copy values of matrix elements to device

    nbytes = Morb * Morb * sizeof(cuDoubleComplex);
    err = cudaMemcpy(d_Ho,Ho,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    nbytes = Morb * Morb * Morb * Morb * sizeof(cuDoubleComplex);
    err = cudaMemcpy(d_Hint,Hint,nbytes,cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Host to Device - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }



    printf("\n\n======================================\n\n");

    // optimally choose the number of blocks
    blocks = (nc + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    for (i = 0; i < 5; i++)
    {
        applyHconf<<<blocks,THREADS_PER_BLOCK>>>(Npar,Morb,d_Map,d_MapOT,d_MapTT,
            d_strideOT,d_strideTT,d_IFmat,d_C,d_Ho,d_Hint,d_out);

        cudaDeviceSynchronize();
    }

    nbytes = nc * sizeof(cuDoubleComplex);
    err = cudaMemcpy(out,d_out,nbytes,cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("\n\nFailed to copy data from Device to Host - ");
        printf(" error code : %s!\n\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // free Host memory
    free(C);
    free(out);
    free(Ho);
    free(Hint);
    free(IFmat);
    free(NCmat);
    free(strideOT);
    free(strideTT);
    free(MapOT);
    free(MapTT);
    free(Map);

    // free device memory

    cudaFree(d_C);
    cudaFree(d_out);
    cudaFree(d_Ho);
    cudaFree(d_Hint);
    cudaFree(d_IFmat);
    cudaFree(d_strideOT);
    cudaFree(d_strideTT);
    cudaFree(d_MapOT);
    cudaFree(d_MapTT);
    cudaFree(d_Map);

    printf("\n\nDone.\n\n");
    return 0;
}
