// Use an input to set input for t>2.0
// else input=1
void inputs(int pas,int ncell,int trial)
{
  int n_gene;
  int input_index;
   
  for (n_gene=0;n_gene<NINPUT;n_gene++)
    {
      input_index = trackin[n_gene];
      if(pas<2.0)
	{
	  history[input_index][pas][ncell] = 1.0;
	}
      else
	{
	  history[input_index][pas][ncell] = 0.0;
	}
    }
}


