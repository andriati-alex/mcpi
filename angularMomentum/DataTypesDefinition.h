#ifndef _DataTypesDefinition_h
#define _DataTypesDefinition_h

#include <complex.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define PI 3.141592653589793
#define MEMORY_TOL 1E9      // SAFEGUARD TO DO NOT EXCEED CAPABILITY
#define LNCZS_STOP_TOL 1E-8 // MAX. VARIATION IN EIGENVALUE TO STOP ALGORITHM
#define LNCZS_TIME_IT 4     // LANCZOS ITERATIONS FOR TIME EVOLUTION
#define MAX_LNCZS_IT 500    // MAX. NUMBER OF LANCZOS ITERATIONS



typedef int * Iarray;               // Vector of integer numbers
typedef char * Farray;              // Integer numbers for Fermi occ.
typedef double * Rarray;            // Vector of real numbers
typedef double ** Rmatrix;          // Matrix of real numbers
typedef double complex * Carray;    // Vector of complex numbers
typedef double complex ** Cmatrix;  // Matrix of complex numbers



struct _HConfMat
{
/** Sparse matrix structure for Hamiltonian matrix in Config. space
    It is commonly feasible for  single  species case for which the
    configurational space is not too large. More information in
    https://en.wikipedia.org/wiki/Sparse_matrix
    Look for the Yale sparse matrix storage                     **/
    int
        nnze;
    Iarray
        rows,
        cols;
    Carray
        vals;
};
typedef struct _HConfMat * HConfMat;



struct _TensorProd_ht
{
/** The Tensor product of two spaces with specific momentum La
    for species A and Lb for species B is a  subspace  of  the
    general combined space with momentum L = La + Lb       **/
    int
        La,     // momentum of config. space A
        Lb,     // momentum of config. space B
        sizeA,  // multiconfig. space size A
        sizeB;  // multiconfig. space size B

    Iarray
        * hta,  // hashing table for config. in space A
        * htb;  // hashing table for config. in space B
};
typedef struct _TensorProd_ht * TensorProd_ht;



struct _CompoundSpace
{
/** The compound space is design by selecting the momenta in  integer
    units separately for species A and B, such that  La + Lb = L.  In
    this way, it is possible to define an array of subspaces for each
    possible combination of angular momenta of each species       **/
    int
        L,     // total angular momentum (from both species)
        Na,
        Nb,
        size,  // total number of elements in compound basis
        n_sub, // number of subspaces - size of arrays below
        lmaxA, // max. momentum for IPS for species A
        lmaxB; // max. momentum for IPS for species B
    // In a contiguous enumeration of the product space,
    // the strides mark in which index  a  new  subspace
    // begins, which are related to momenta of  A and  B
    Iarray
        strides;
    struct _TensorProd_ht
        * sub; // subspaces all possible Tensor product spaces
};
typedef struct _CompoundSpace * CompoundSpace;



struct _BFTensorProd_ht
{
/** The Tensor product of two spaces with specific momentum La
    for species A and Lb for species B is a  subspace  of  the
    general combined space with momentum L = Lb + Lf       **/
    int
        Lb,     // momentum of config. space Bosons
        Lf,     // momentum of config. space Fermions
        sizeB,  // multiconfig. space size Bosons
        sizeF;  // multiconfig. space size Fermions
    Iarray
        * htb;  // hashing table for config. of Bosons
    Farray
        * htf;  // hashing table for config. of Fermions
};
typedef struct _BFTensorProd_ht * BFTensorProd_ht;



struct _BFCompoundSpace
{
/** The compound space is design by selecting the momenta in  integer
    units separately for species B and F, such that  Lb + Lf = L.  In
    this way, it is possible to define an array of subspaces for each
    possible combination of angular momenta of each species       **/
    int
        L,     // total angular momentum (from both species)
        Nb,
        Nf,
        size,  // total number of elements in compound basis
        n_sub, // number of subspaces - size of arrays below
        lmaxB, // max. momentum for Bosons   single particle states
        lmaxF; // max. momentum for Fermions single particle states
    // In a contiguous enumeration of the product space,
    // the strides mark in which index  a  new  subspace
    // begins, which are related to momenta of  B and  F
    Iarray
        strides;
    struct _BFTensorProd_ht
        * sub; // subspaces all possible Tensor product spaces
};
typedef struct _BFCompoundSpace * BFCompoundSpace;



#endif
