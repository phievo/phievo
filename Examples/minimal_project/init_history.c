// Sets gene concentration at t=0, this function is used before starting the integration.

void init_history(int trial)  {
  int ncell,n_gene;
  int input_index;
  int output_index;
  // Everything to 0
  for (ncell=0;ncell<NCELLTOT;ncell++)
    {
      for (n_gene=0;n_gene<SIZE;n_gene++)
	{
	  history[n_gene][0][ncell] = 0;
        }
    }

  // Set input Species to 1 (for the sake of the example)
  for (ncell=0;ncell<NCELLTOT;ncell++)
    {
      for (n_gene=0;n_gene<NINPUT;n_gene++)
	{
	  input_index = trackin[n_gene]; // trackin is a static array(accessible from every function) 
	  history[input_index][0][ncell] = 0;
        }
    }

  // Set output Species to 0.5 (for the sake of the example)
  for (ncell=0;ncell<NCELLTOT;ncell++)
    {
      for (n_gene=0;n_gene<NOUTPUT;n_gene++)
	{
	  output_index = trackout[n_gene]; // trackout is a static array
	  history[output_index][0][ncell] = 0.5;
        }
    }
}
