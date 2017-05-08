
/* compute the RHS of equations and run over NSTEP's with 1st order Euler method
   The arugment kk, is an index passed to inputs that records which iteration of
   initial or boundary conditions the same system of equs is being integrated for

   This version of code used with ligand at cell computed from sum neighbors
*/

void integrator(int kk){

    double s[SIZE];
    double ds[SIZE];
    double sumligands[SIZE];
    double memory[SIZE];
    int index,n,pas,ncell;

    for (index=0;index<SIZE;index++){
	s[index] = 0;
        ds[index]=0;
        memory[index]=0;
    }

    /* initialize geometry here, incase cells move  */
    init_geometry();
    init_history(kk);

    /* loop over time steps, then over each cell etc */
    for (pas=0;pas<NSTEP-1;pas++)  {
	for (ncell=0;ncell<NCELLTOT;ncell++)  {
            inputs(pas,ncell,kk);
            for (index=0;index<SIZE;index++) {
	        s[index]=history[index][pas][ncell];
            }
            derivC(s,history,pas,ds,memory,ncell);  //local integration
            sum_concentration(ncell,pas,sumligands);  //perform sum of ligands concentrations for non external ligands
	    diffusion(ncell,pas,ds,history,geometry);//computes diffusion of external ligands
            LRinC(s,ds,sumligands);

            for (index=0;index<SIZE;index++) {
	 	 history[index][pas+1][ncell] = s[index] + DT*ds[index];
		 if (history[index][pas+1][ncell]<0)//might happen for langevin
		   history[index][pas+1][ncell]=0;
	    }
	}
    }

    /* fill in inputs for last time.  */
    for (ncell=0;ncell<NCELLTOT;ncell++)  {
      inputs(NSTEP-1,ncell,kk);
    }
}
