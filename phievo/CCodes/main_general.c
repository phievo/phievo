/*
General main.c, should not need any others, see comments in fitness_template for
the functions that must be supplied in that file for the main to work,
*/

int main()  {
    srand( SEED );
    int i,k,l;
    double score = 0;
    int len_hist = SIZE*NCELLTOT*NSTEP;
    double *hptr = &history[0][0][0];
    /* following incase one wants to average history before doing fitness */
    static double history2[SIZE][NSTEP][NCELLTOT];
    double *h2ptr = &history2[0][0][0];//table for averaging output (used for multicell problems)

    /* dummy return when no outputs for fitness function */
    if(NOUTPUT <= 0) {
        printf("%s","no output variables? terminating without integration" );
    }

    for(i=0; i<len_hist; i++)  {
        *(h2ptr + i) = 0;
    }

    for (k=0; k<NTRIES; k++){
        integrator(k);
        fitness(history, trackout,k);
        if( PRINT_BUF )  {
            print_history(k);
        }

        for(i=0; i<len_hist; i++)  {
            *(h2ptr+i) += *(hptr +i);
        }
    }
    treatment_fitness(history2,trackout);
}
