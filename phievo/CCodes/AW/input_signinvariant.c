/* define input variables as a function of time.  n-attempts = loop [0,ntries)  */

void inputs(int pas,int ncell, int ntry){
  if( pas <=NSTEP/2)  history[trackin[0]][pas][ncell]= (double) 2.0;//input at the beginning
  else if( pas > NSTEP/2 && ntry == 0)  history[trackin[0]][pas][ncell]= (double) 4.0;//input at the end -- try1
  else if( pas > NSTEP/2 && ntry == 1)  history[trackin[0]][pas][ncell]= (double) 1.0;//input at the end -- try2
}


