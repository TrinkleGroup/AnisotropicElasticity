#ifndef __NNPAIR_H
#define __NNPAIR_H

/*
  Program: nnpair.H
  Author:  D. Trinkle
  Date:    Jan. 23, 2003
  Purpose: The header file for ident.C; contains structures and
           functions for handling the nearest neighbor grid construction,
	   
           This header file contains the following subroutines:

	   calc_grid(cart[9], Rcut, Ngrid[3])
	     --determines the grid sizing in each direction, given
	       the cell dimensions in cart, and the cutoff distance

	   make_grid(Ngrid[3], grid_list[])
	     --constructs the grid elements (including connectivity
	       information).

	   populate_grid(Ngrid[3], grid_list[], Natoms, u[])
	     --bins up all of the atoms into the grid element list
	       (so that we can construct the nearest neighbor list).

	   free_grid(Ngrid[3], grid_list[])
	     --frees up the grid list

	   nn_grid (cart[9], Ngridelem, grid_list[], Natoms, u[],
	            Rcut, NNpairs, nn_list[])
	     --constructs the list of nearest neighbor pairs (basically,
	       a bond list).  The list is symmetric--so both the i-j
	       and j-i bonds appear in the list.

	   nn_raw (cart[9], Natoms, u[], Rcut, NNpairs, nn_list[])
	     --constructs the list of nearest neighbor pairs, but
	       does it the old-fashioned way.  Useful only for small
	       lattices where Rcut is large compared to cart[i].  This
	       will search beyond just the neighboring cells for neighbors.

	   sort_nn_list(NNpairs, nn_pair_list[], Natoms, nn_list[][])
	     --takes the list of bonds, and makes a sorted matrix
	       style list out of them, where for each atom i:
        	 nn_list[i][0]    : number of neighbors
		 nn_list[i][1..?] : index of the bond in nn_pair_list
	       We just give the index to the bond elem in nn_pair_list.
	       The list is sorted from shortest to longest bond length.

	   free_nn_list(Natoms, nn_list[])
	     --frees up the matrix style list of nearest neighbors.
*/

//****************************** STRUCTURES *****************************
// Done this way to make the linked list kinda structure work:
typedef struct
{
  int Nneigh;       // Number of neighboring grid_elem's
  int* neighlist;   // List of neighboring grid_elem's
  int Natoms;       // Number of atoms in grid block
  int* atomlist;    // List of atoms in grid block
} grid_elem_type;


typedef struct 
{
  int i, j;       // The two atoms in the nn_pair
  double r;       // Length of bond
  double v_ij[3]; // Unit length vector pointing from i to j (r_j-r_i)
} nn_pair_type;


//****************************** SUBROUTINES ****************************
// Determine the number of grid elements in each direction:
void calc_grid(double cart[9], double Rcut, int Ngrid[3]);

// Make the grid (connect the grid blocks)
void make_grid(int Ngrid[3], grid_elem_type* &grid_list);

// Delete the grid
void free_grid(int Ngrid[3], grid_elem_type* &grid_list);

// Put the atoms in the grid blocks (according to grid_elem function)
void populate_grid(int Ngrid[3], grid_elem_type* grid_list, 
		   int Natoms, double** u);

// Using the grid list, construct the nn list
void nn_grid (double cart[9], int Ngridelem, grid_elem_type* grid_list, 
	      int Natoms, double **u, double Rcut,
	      int &NNpairs, nn_pair_type* &nn_list);

// Without a grid list, construct the nn list.  Also searches beyond
// the neighboring periodic image cells to make bonds.
void nn_raw (double cart[9], int Natoms, double **u, double Rcut,
	     int &NNpairs, nn_pair_type* &nn_list);


// Rearranging the pairing list into a sorted "matrix" style list
void sort_nn_list(int NNpairs, nn_pair_type* nn_pair_list,
		  int Natoms, int** &nn_list);

void free_nn_list(int Natoms, int** &nn_list);



//******************************* calc_grid *****************************
// Determine the number of grid elements in each direction:
void calc_grid(double cart[9], double Rcut, int Ngrid[3]) 
{
  int i, j, k, l;
  double det_cart;
  double len[3];
  double ajXak[3];
  double ajXak_len;

  // a1.(a2 x a3)
  det_cart = cart[0]*(cart[4]*cart[8] - cart[7]*cart[5]) -
    cart[1]*(cart[3]*cart[8] - cart[6]*cart[5]) +
    cart[2]*(cart[3]*cart[7] - cart[6]*cart[4]);
  if (det_cart < 0) det_cart = -det_cart;

  // Calculate aj x ak foreach i.
  for (i=0; i<3; ++i) {
    j = (i+1)%3;
    k = (i+2)%3;
    for (l=0; l<3; ++l)
      ajXak[l] = cart[3*j + (l+1)%3] * cart[3*k + (l+2)%3]
	- cart[3*j + (l+2)%3] * cart[3*k + (l+1)%3];
    ajXak_len = 0;
    for (l=0; l<3; ++l) ajXak_len += ajXak[l]*ajXak[l];
    len[i] = det_cart / sqrt(ajXak_len);
  }
  
  for (i=0; i<3; ++i) {
    Ngrid[i] = (int) (len[i]/Rcut);
    if (Ngrid[i] == 0) Ngrid[i] = 1;
  }
}

//******************************* make_grid ****************************
// Make the grid (connect the grid blocks)

// Convert from a three integer triplet to a single index
inline int trip2int (int Ngrid[3], int i[3]) 
{
  return i[0] + Ngrid[0]*(i[1] + Ngrid[1]*i[2]);
}

// Returns which grid a given u_vect lies inside (used by populate)
inline int grid_elem(int Ngrid[3], double u_vect[3]) 
{
  return (int)(u_vect[0]*Ngrid[0]) 
    + Ngrid[0]*((int)(u_vect[1]*Ngrid[1])
		+ Ngrid[1]*((int)(u_vect[2]*Ngrid[2])));
}


void make_grid(int Ngrid[3], grid_elem_type* &grid_list) 
{
  int i[3], j[3], jtemp[3], k;
  int index, n_count;
  int Nelem;
  int nneigh;

  Nelem = Ngrid[0]*Ngrid[1]*Ngrid[2];
  grid_list = new grid_elem_type[Nelem];

  // Cute little piece of code to handle Ngrid = 1, 2, or greater.
  int start[3], end[3];

  nneigh = 1;
  for (k=0; k<3; ++k) {
    if (Ngrid[k] == 1) { start[k] =  0; end[k] = 0; }
    if (Ngrid[k] == 2) { start[k] =  0; end[k] = 1; nneigh *= 2;}
    if (Ngrid[k] >= 3) { start[k] = -1; end[k] = 1; nneigh *= 3;}
  }

  for (i[0]=0; i[0]<Ngrid[0]; ++(i[0]))
    for (i[1]=0; i[1]<Ngrid[1]; ++(i[1]))
      for (i[2]=0; i[2]<Ngrid[2]; ++(i[2])) {
	index = trip2int(Ngrid, i);
	grid_list[index].Nneigh = nneigh;
	grid_list[index].neighlist = new int[nneigh];
	grid_list[index].Natoms = 0;
	grid_list[index].atomlist = NULL;
	// Now, do the connections:
	n_count = 0;
	for (j[0]=(i[0]+start[0]); j[0]<=(i[0]+end[0]); ++(j[0]))
	  for (j[1]=(i[1]+start[1]); j[1]<=(i[1]+end[1]); ++(j[1]))
	    for (j[2]=(i[2]+start[2]); j[2]<=(i[2]+end[2]); ++(j[2])) {
	      for (k=0; k<3; ++k) jtemp[k] = (j[k]+Ngrid[k]) % (Ngrid[k]);
	      grid_list[index].neighlist[n_count] = trip2int(Ngrid, jtemp);
	      ++n_count;
	    }
      }
}


//***************************** populate_grid **************************
// Put the atoms in the grid blocks (according to grid_elem function)
void populate_grid(int Ngrid[3], grid_elem_type* grid_list, 
		   int Natoms, double** u) 
{
  int i, k, n;
  double approx_density;
  int Ngridelem;
  int* Nalloc, Nalloc0; // Temporary variable.
  int* t_list;
  grid_elem_type* g;

  Ngridelem = Ngrid[0]*Ngrid[1]*Ngrid[2];
  approx_density = (double)Natoms / (double)Ngridelem;
  // Let's assume that the density at best gets to be twice as big.
  Nalloc0 = (int)(2.*approx_density);
  if (Nalloc0 < 4) Nalloc0 = 4;
  Nalloc = new int[Ngridelem];
  for (i=0; i<Ngridelem; ++i) {
    grid_list[i].Natoms = 0;
    grid_list[i].atomlist = new int[Nalloc0];
    Nalloc[i] = Nalloc0;
  }
  
  for (n=0; n<Natoms; ++n) {
    // Where does this atom belong?
    i = grid_elem(Ngrid, u[n]);
    g = grid_list + i;
    if (g->Natoms == Nalloc[i]) {
      // We need to reallocate, so let's do it:
      Nalloc[i] = 2*Nalloc[i]; // Double memory
      t_list = g->atomlist;
      g->atomlist = new int[Nalloc[i]];
      // Copy elements
      for (k=0; k<(g->Natoms); ++k) g->atomlist[k] = t_list[k];
      delete[] t_list; // Delete old one
    }
    // Now, we can insert it:
    g->atomlist[g->Natoms] = n;
    ++(g->Natoms);
  }
  // Quick sanity check for populate:
  for (n=0, i=0; i<Ngridelem; ++i) n += grid_list[i].Natoms;
  if (n!= Natoms) {
    fprintf(stderr, "Error occured in populate grid: found %d atoms, not %d atoms\n", n, Natoms);
  }
  
  delete[] Nalloc;
}

//***************************** free_grid ******************************
// Delete the grid
void free_grid(int Ngrid[3], grid_elem_type* &grid_list) 
{
  int i;
  int Ngridelem;

  if (grid_list == NULL) return;

  Ngridelem = Ngrid[0]*Ngrid[1]*Ngrid[2];
  for (i=0; i<Ngridelem; ++i) {
    delete[] grid_list[i].atomlist;
    delete[] grid_list[i].neighlist;
  }
  delete[] grid_list;
  grid_list = NULL;
}


//****************************** nn_grid *******************************
// Using the grid list, construct the nn list

// Keeps the difference between two numbers between -0.5 and 0.5.
inline double diff (double x) {return x + 2 - (int)(x+2.5);}

void nn_grid (double cart[9], int Ngridelem, grid_elem_type* grid_list, 
	      int Natoms, double **u, double Rcut,
	      int& NNpairs, nn_pair_type* &nn_list) 
{
  int ng, ngn, i, ii, j, jj, k;
  grid_elem_type *g, *gn;
  double vect[3];
  double r2, Rcut2;
  
  Rcut2 = Rcut*Rcut;

  // Determine the maximum number of nn pairs allowed:
  int natoms, npair, nneigh;
  natoms = 0; nneigh = 0; // Find the maximum number of atoms and grid neigh's
  for (ng=0; ng<Ngridelem; ++ng) {
    if (grid_list[ng].Natoms > natoms) natoms = grid_list[ng].Natoms;
    if (grid_list[ng].Nneigh > nneigh) nneigh = grid_list[ng].Nneigh;
  }
  // Assume that every atom in each of the boxes is a neighbor:
  npair = natoms*nneigh*Natoms;
  nn_list = new nn_pair_type[npair];

  // Now, go through each grid element, and find the pairs:
  npair = 0;
  for (ng=0; ng<Ngridelem; ++ng) {
    g = grid_list + ng;
    // Loop through the atoms in this grid:
    for (ii=0; ii<(g->Natoms); ++ii) {
      i = (g->atomlist)[ii];
      // Go through the neighboring boxes:
      for (ngn=0; ngn<(g->Nneigh); ++ngn) {
	gn = grid_list + g->neighlist[ngn];
	for (jj=0; jj<(gn->Natoms); ++jj) {
	  j = (gn->atomlist)[jj];
	  if (i==j) continue;
	  // Convert into cartesian coord.:
	  r2 = 0;
	  for (k=0; k<3; ++k) {
	    vect[k] = cart[k]*(diff(u[j][0] - u[i][0]))
	      + cart[3+k]*(diff(u[j][1] - u[i][1]))
	      +	cart[6+k]*(diff(u[j][2] - u[i][2]));
	    r2 += vect[k]*vect[k];
	  }
	  if (r2 > Rcut2) continue;
	  // Now, we've got a pair... let's add it to the list:
	  nn_list[npair].i = i;
	  nn_list[npair].j = j;
	  nn_list[npair].r = sqrt(r2);
	  r2 = 1./sqrt(r2);
	  for (k=0; k<3; ++k) nn_list[npair].v_ij[k] = r2*vect[k];
	  ++npair;
	}
      }
    }
  }
  NNpairs = npair;
}


//******************************** nn_raw ******************************
// Without a grid list, construct the nn list.  Also searches beyond
// the neighboring periodic image cells to make bonds.
void nn_raw (double cart[9], int Natoms, double **u, double Rcut,
	     int &NNpairs, nn_pair_type* &nn_list) 
{
  int i, j, k;
  double vect[3];
  double r2, Rcut2;

  Rcut2 = Rcut*Rcut;

  // First, let's find the smallest vector we can make with our
  // cartesian cell, and then we'll determine how many pairs, etc.
  double min_dist;
  int n[3], maxn;
  min_dist = 1.e20;
  const int RANGE = 3;
  for (n[0]=-RANGE; n[0]<=RANGE; ++(n[0]))
    for (n[1]=-RANGE; n[1]<=RANGE; ++(n[1]))
      for (n[2]=-RANGE; n[2]<=RANGE; ++(n[2])) {
	if ( (n[0] == 0) && (n[1] == 0) && (n[2] == 0) ) continue;
	for (k=0; k<3; ++k)
	  vect[k] = n[0]*cart[k] + n[1]*cart[3+k] + n[2]*cart[6+k];
	r2 = vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2];
	if (r2 < min_dist) min_dist = r2;
      }
  maxn = (int) (Rcut / sqrt(min_dist) + 1.);
  
  int npair;
  npair = (2*maxn + 1);
  npair = Natoms*Natoms*npair*npair*npair;
  nn_list = new nn_pair_type[npair];
  
  // Let's pair 'em up!
  npair = 0;
  for (i=0; i<Natoms; ++i) {
    // n[] is the matrix of shifts.
    for (n[0]=-maxn; n[0]<=maxn; ++(n[0]))
      for (n[1]=-maxn; n[1]<=maxn; ++(n[1]))
	for (n[2]=-maxn; n[2]<=maxn; ++(n[2]))
	  for (j=0; j<Natoms; ++j) {
	    if ( (i==j) && (n[0] == 0) && (n[1] == 0) && (n[2] == 0) )
	      continue;
	    // Convert into cartesian coord.:
	    r2 = 0;
	    for (k=0; k<3; ++k) {
	      vect[k] = cart[k]*(u[j][0] - u[i][0] + n[0])
		+ cart[3+k]*(u[j][1] - u[i][1] + n[1])
		+ cart[6+k]*(u[j][2] - u[i][2] + n[2]);
	      r2 += vect[k]*vect[k];
	    }
	    if (r2 > Rcut2) continue;
	    // Now, we've got a pair... let's add it to the list:
	    nn_list[npair].i = i;
	    nn_list[npair].j = j;
	    nn_list[npair].r = sqrt(r2);
	    r2 = 1./sqrt(r2);
	    for (k=0; k<3; ++k) nn_list[npair].v_ij[k] = r2*vect[k];
	    ++npair;
	  }
  }
  NNpairs = npair;
}


//****************************** sort_nn_list **************************
// Rearranging the pairing list into a sorted "matrix" style list:
// nn_list[i][0]    : number of neighbors
// nn_list[i][1..?] : index of the bond in nn_pair_list
//   Our goal is to sort this list from shortest to longest bond for
//   each i.

void sort_pair_list (nn_pair_type* nn_pair_list, int Npairs, int* list) 
{
  // Do a quick type insertion sort on list:
  int i, j, t;
  for (i=0; i<Npairs; ++i)
    for (j=(Npairs-1); j>i; --j)
      if (nn_pair_list[list[j]].r < nn_pair_list[list[j-1]].r) {
	t = list[j];
	list[j] = list[j-1];
	list[j-1] = t;
      }
}

void sort_nn_list(int NNpairs, nn_pair_type* nn_pair_list,
		  int Natoms, int** &nn_list) 
{
  int i, inew, j, p, count;
  int* tlist;

  // Make a "temporary" list
  i = (10*NNpairs) / Natoms;  if (i < 10) i = 10;
  tlist = new int[i];
  
  if (nn_list == NULL) nn_list = new int *[Natoms];

  // Now, let's make those pairings:
  // Start us off 
  p = 0; count = 0;
  inew = nn_pair_list[p].i;
  i = inew;
  // Note: the exit for this loop is inside of it (inew == -1)
  while (-1) {
    if (inew == i) {
      tlist[count] = p;
      ++count;
      ++p;
      if (p != NNpairs) inew = nn_pair_list[p].i;
      else              inew = -1;
    }
    else {
      // Sort our list:
      sort_pair_list(nn_pair_list, count, tlist);
      // Add to our list:
      nn_list[i] = new int[count+1];
      nn_list[i][0] = count;
      for (j=0; j<count; ++j) nn_list[i][j+1] = tlist[j];
      i = inew;
      count = 0;
      if (inew == -1) break; // Our flag to get out...
    }
  }
  // Garbage collection:
  delete[] tlist;
}


//****************************** free_nn_list **************************
void free_nn_list(int Natoms, int** &nn_list) 
{
  int i;
  for (i=0; i<Natoms; ++i) delete[] nn_list[i];
  delete[] nn_list;
  nn_list = NULL;
}


#endif
