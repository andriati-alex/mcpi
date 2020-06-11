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



void fun_rec(int L, int l, int * vec, int vec_size, int * countAdd)
{

/** RECURSION TO BUILD THE ADDITIONS
    --------------------------------
    The first argument 'L' is how much need to be added yet.  The
    second 'l' is the last number added. 'vec' stack up the terms
    already used. **/

    int
        x,
        init;

    if (L == 0)
    {
        *countAdd = *countAdd + 1;
        printPoss(vec_size,vec);
        return;
    }

    if (L > l) init = l;
    else       init = L;

    for (x = init; x > 0; x--)
    {
        vec[x-1] = vec[x-1] + 1;
        fun_rec(L-x,x,vec,vec_size,countAdd);
        vec[x-1] = vec[x-1] - 1;
    }
}



int printDecomp(int L, int cutoff)
{

/** PRINT ALL POSSIBLE SUMS OF INTEGERS THAT YIELD THE SAME RESULT
    --------------------------------------------------------------
    Given an integer number 'L' the basic ideia is to decompose as
    L = x + (L - x)
    for any x = L, ..., 1,
    L = (L - x) + (x - y) + y
    with y = x-1, ..., 1, and so on. We can do it recursively.  In
    this way,  we can build all possible sums that  yield  'L'  as
    result.  There is a possibility to add a  maximum  value for x
    different of 'L', given by the 'cutoff'.                   **/

    int
        i,
        l,
        Npos,
        startNum;

    int
        * vec;

    printf("\n");

    // Each index of the vector is a possible integer number that can
    // be used in the decomposition. The entry v[i-1] counts how many
    // times the same number 'i' appears in the addition.
    // 'Npos' counts how many possible decomposition there are.
    vec = (int *) malloc(L*sizeof(int));
    Npos = 0;

    // cutoff > L does not make sense. Filter this possible entry
    if (L > cutoff) startNum = cutoff;
    else            startNum = L;

    for (l = startNum; l > 0; l--)
    {
        for (i = 0; i < L; i++) vec[i] = 0;
        vec[l-1] = 1;
        fun_rec(L-l,l,vec,L,&Npos); // recursion
    }

    free(vec);

    return Npos;
}



int main(int argc, char * argv[])
{

/** MAIN FUNCTION - PRINT ALL POSSIBLE ADDITIONS THAT YEILD THE SAME RESULT **/

    int
        cutoff,
        sumValue,
        possibilities;

    if (argc != 3)
    {
        printf("\n\nERROR: Need two command line arguments : \n");
        printf("\t1. Integer number\n");
        printf("\t2. Max. integer number allowed in decomposition\n\n");
        exit(EXIT_FAILURE);
    }

    sscanf(argv[1],"%d",&sumValue);
    sscanf(argv[2],"%d",&cutoff);
    possibilities = printDecomp(sumValue,cutoff);

    printf("\nTotal number of possibilities : %d",possibilities);

    printf("\n\n");
    return 0;
}
