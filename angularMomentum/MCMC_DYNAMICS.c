#include <time.h>
#include "TimePropagation.h"

/** MAIN PROGRAM TO COMPUTE GROUND STATE(S) FROM INPUT FILES (.inp)

        **************************************************************
        **                                                          **
        **  DEVELOPER                                               **
        **      Alex Andriati   andriati@if.usp.br                  **
        **  FILIATION                                               **
        **      University of Sao Paulo - Physics Institute         **
        **  PAGES                                                   **
        **      https://www.researchgate.net/profile/Alex_Andriati  **
        **      https://github.com/andriati-alex                    **
        **                                                          **
        **************************************************************

Use external files to set up the problem of  computing ground state
and properly record the output data. It requires a .inp file  whose
name must be provide in the job configuration  file  'job.conf'  as
well as the species type and quantity (1 or 2). The file 'job.conf'
provides comments in sentences preceded by # explaining  the  input
keywords and numbers. By default a scanning is evaluate  using  all
the lines of .inp file as a different problem set up

COMPILATION
icc main.c -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -qopenmp -O3 -o exe

**/



int main(int argc, char * argv[])
{

    char
        c,
        typeA,  // Species A is aways bosons
        typeB;  // Can be either (b)osons or (f)ermions

    char
        TimeInfo[20],
        spec_name[20],
        inp_fname[100],
        fname_prefix[100];

    int
        i,
        Njobs,  // number of lines as input parameters to evaluate
        Nspec,  // number of species must be 1 or 2
        nthreads;

    double
        time_used,
        record_interval;

    FILE
        * job_file;



    // set up parallel info. based on workstation capabilities
    nthreads = omp_get_max_threads()/2;
    omp_set_num_threads(nthreads);



    /* CONFIGURE TYPE OF SYSTEM
    ==================================================================== */
    job_file = openFileRead("dynamics-job.conf");
    i = 1;
    while ((c = getc(job_file)) != EOF)
    {
        // jump comment line. They are start with #
        if (c == '#') { ReachNewLine(job_file); continue; }
        else          { fseek(job_file,-1,SEEK_CUR);    }

        switch (i)
        {
            case 1:
                fscanf(job_file,"%d%c%c",&Nspec,&typeA,&typeB);
                i = i + 1;
                break;
            case 2:
                fscanf(job_file,"%s",fname_prefix);
                i = i + 1;
                break;
            case 3:
                fscanf(job_file,"%lf",&record_interval);
                i = i + 1;
                break;
            case 4:
                fscanf(job_file,"%s",TimeInfo);
                break;
        }
        ReachNewLine(job_file);
    }
    fclose(job_file);
    /* FINISHED READING OF CONFIGURATION PARAMETERS
    =================================================================== */

    // Assess input parameter
    if (i != 4)
    {
        printf("\n\nINPUT ERROR : not enough input lines in job.conf\n\n");
        exit(EXIT_FAILURE);
    }
    if (Nspec != 2 && Nspec != 1)
    {
        printf("\n\nINPUT ERROR : The number of species must be 1 or 2\n\n");
        exit(EXIT_FAILURE);
    }

    // set up file name for output data
    strcpy(inp_fname,fname_prefix);
    strcat(inp_fname,".inp");
    Njobs = NumberOfLines(inp_fname);

    printf("\n\n");
    printf("\t******************************************************\n");
    printf("\t*                                                    *\n");
    printf("\t*     FIXED MOMENTUM MULTICONFIGURATIONAL METHOD     *\n");
    printf("\t*                  TIME PROPAGATION                  *\n");
    printf("\t*                                                    *\n");
    printf("\t******************************************************\n");

    printf("\nJOB SET UP INFORMATION");
    if (Nspec == 1) printf("\n\tSingle species of bosons");
    else
    {
        if (typeB == 'f') strcpy(spec_name,"fermions");
        else              strcpy(spec_name,"bosons");
        printf("\n\tMixture : bosons + %s",spec_name);
    }
    printf("\n\tFrom input file %s %d jobs are requested",inp_fname,Njobs);

    // CALLING DRIVER ROUTINE TO COMPUTE GROUND-STATE
    time_used = omp_get_wtime(); // trigger to measure time
    switch (Nspec)
    {
        case 1:
            TIME_SCANNING(Njobs,fname_prefix,record_interval,TimeInfo[0]);
            break;
        case 2:
            if (typeB == 'F' || typeB == 'f')
            {
                TIME_BOSEFERMI_SCANNING(Njobs,fname_prefix,
                    record_interval,TimeInfo[0]);
            }
            else
            {
                TIME_MIXTURE_SCANNING(Njobs,fname_prefix,
                    record_interval,TimeInfo[0]);
            }
            break;
    }
    time_used = omp_get_wtime() - time_used; // finish time measure

    printf("\n\n\n");
    sepline();
    printf("\nTotal Time Elapsed (with %d threads) : ",nthreads);
    printf("%.1lf(s) = ",time_used);
    TimePrint(time_used);
    if (Njobs > 1)
    {
        printf("\nAverage time elapsed per ground state computed : ");
        printf("%.1lf(s) = ",time_used/Njobs);
        TimePrint(time_used);
    }

    printf("\n\nDone !\n\n");

    return 0;
}
