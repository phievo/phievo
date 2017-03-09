/* define input variables as a function of time.  n-attempts = loop [0,ntries)  */

void inputs(int pas,int ncell, int ntry){
  history[trackin[0]][pas][ncell]= (double) 10.0*ntry/9;//input at the beginning
}


