/* Define geometry[NCELLTOT][NNEIGHBOR] array defining the neighbors of each cell
   Convention geometry < 0 ie invalid index -> no neighbor, ie used for end cells
   NB even for NCELLTOT=1, should have NNEIGHBOR == 3 to avoid odd seg faults since
   3 neighbors accessed in various loops, eventhough does nothing.
*/
   
void init_geometry() {

  int index;

  for (index=0;index<NCELLTOT;index++){
    geometry[index][0]=index-1;//NB no left neighbour for cell 0 (-1)
    geometry[index][1]=index;
    geometry[index][2]=index+1;
  }

  geometry[NCELLTOT-1][2]=-1;//NB no right neighbour for cell ncell-1 (-1)

}


