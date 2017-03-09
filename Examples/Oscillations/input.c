/* define input variables as a function of step, which cell, and n-attempts = loop [0,ntries)
set history[input# ][pas][ncell] to correct value.
*/

void inputs(int pas,int ncell,int n_attempts){
    int track0 = trackin[0];
    int track1 = trackin[1];
    history[track0][pas][ncell] = isignal[pas][ncell][0];
    history[track1][pas][ncell] = isignal[pas][ncell][1];
}
