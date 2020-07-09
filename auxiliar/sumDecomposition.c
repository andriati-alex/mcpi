/** PROGRAM TO FIND ALL POSSIBLE ADDITIONS OF POSITIVE INTEGER NUMBERS WHICH
    YEILD THE SAME RESULT
    ========================================================================

    Developed by : Alex Andriati
    email        : andriati@if.usp.br
    University of Sao Paulo, Physics Institute
    jun/2020
    ======================================================================== **/



#include <stdio.h>
#include <stdlib.h>



void printPoss(int n, int * vec)
{

/** PRINT ONE OF MANY POSSIBLE ADDITIONS THAT YIELD THE SAME RESULT **/

    int
        k,
        i,
        j;

    for (k = n; k > 0; k--)
    {
        if (vec[k-1] > 0) break;
    }

    printf("%d",k);
    for (j = 1; j < vec[k-1]; j++) printf(" + %d",k);

    for (i = k-1; i > 0; i--)
    {
        for (j = 0; j < vec[i-1]; j++) printf(" + %d",i);
    }
    printf("\n");
}



void fun_rec(int L, int l, int * rep, int rep_size, int * countAdd)
{

/** RECURSION TO BUILD THE ADDITIONS
    --------------------------------
    The first argument 'L' is how much need to be added yet.  The
    second 'l' is the last number added, and used as upper bound,
    in order to skip permutations of previous results.   rep[i-1]
    counts the repetitions of the number 'i' in the decomposition **/

    int
        x,
        init;

    if (L == 0)
    {
        *countAdd = *countAdd + 1;
        printPoss(rep_size,rep);
        return;
    }

    if (L > l) init = l;
    else       init = L;

    for (x = init; x > 0; x--)
    {
        rep[x-1] = rep[x-1] + 1;
        fun_rec(L-x,x,rep,rep_size,countAdd);
        rep[x-1] = rep[x-1] - 1;
    }
}



int printDecomp(unsigned int L, unsigned int cutoff)
{

/** PRINT ALL POSSIBLE SUMS OF POSITIVE INTEGERS THAT YIELD THE SAME RESULT
    -----------------------------------------------------------------------
    Given a positive integer number 'L' the basic ideia is to decompose as
    L = x + (L - x)
    for any x = L, ..., 1,
    L = (L - x) + (x - y) + y
    with y = x-1, ..., 1, and so on. We can do it recursively. In this way,
    we can build all possible sums that yield 'L' as result.  The  'cutoff'
    applies an upper bound for the largest number in the sum,  thus it only
    has effect if 'cutoff' < 'L' **/

    int
        i,
        l,
        Npos,
        largestNum;

    int
        * rep;

    printf("\n");

    // cutoff > L does not make sense, since the largest number
    // in the sum cannot exceed 'L'. Filter this possibility
    if (L > cutoff) largestNum = cutoff;
    else            largestNum = L;

    // Each index of the 'rep' is a positive integer number that  can
    // be used in the decomposition, thus its size is constrained  by
    // the smallest between 'L' and 'cutoff'. The entry v[i-1] counts
    // how many times the same number 'i' appears in the addition, i.e
    // repetitions in the terms being added.
    rep = (int *) malloc(largestNum*sizeof(int));
    Npos = 0;

    // From the largest possible number that can be used,  recursively
    // search for possible additions that yield the remaining quantity
    // yet to be added, in order to obtain 'L'
    for (l = largestNum; l > 0; l--)
    {
        for (i = 0; i < largestNum; i++) rep[i] = 0;
        rep[l-1] = 1; // (L - l) is how much we need to add yet
        fun_rec(L-l,l,rep,largestNum,&Npos); // recursive function
    }

    free(rep);

    return Npos;
}



int main(int argc, char * argv[])
{

/** MAIN FUNCTION - PRINT ALL POSSIBLE ADDITIONS THAT YEILD THE SAME RESULT **/

    unsigned int
        cutoff,
        sumValue,
        possibilities;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two command line arguments : \n");
        printf("\t1. Positive Integer number - Additions result\n");
        printf("\t2. Max. integer number allowed in decomposition\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&sumValue);
    sscanf(argv[2],"%d",&cutoff);
    possibilities = printDecomp(sumValue,cutoff);

    printf("\nNumber of possible decompositions : %d",possibilities);

    printf("\n\n");
    return 0;
}
