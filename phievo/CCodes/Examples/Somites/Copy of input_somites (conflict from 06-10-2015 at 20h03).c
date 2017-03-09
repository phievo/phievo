/* define input variables as a function of time.  n-attempts = loop [0,ntries) */

void inputs(int pas,int ncell,int n_attempts){
  int trackg1=trackin[0];
  double ss;
 
     ss=0.04*exp(0.75*0.33*0.33*(ncell)-0.6*pas*DT);
     if (ss>1) 
     	ss=1;
     else 
     	ss=0;
     if ((pas>10000)&&(ss>0))
     	ss=0.6;
     history[trackg1][pas][ncell]=ss;

}


