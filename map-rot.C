/*
  Program: map.C
  Author:  D. Trinkle
  Date:    August 12, 2003
  Purpose: Read in two XYZ files and make the differential displacement
           map.  Output to a fig file for plotting.  We rotate our
	   coordinates according to 

  Param.:  <disl> <perfect-xtal.file> <dislocated-xtal.file> <Rcut> <Rmax>
           perfect-xtal:     XYZ file for perfect crystal
	   dislocated-xtal:  XYZ file with dislocation
	   Rcut:             cutoff for 2d NN search
	   Rmax:             radius around COM; if 0, not used.
	   -e    show edge components with direction and center changes
	   -n    write text amount of burgers vector for each pair
	   -b    write text for mini-burgers loops on triads


  Flags:   MEMORY:  the amount of space allocated; not used.
	   VERBOSE: verbose?
	   TESTING: output practically everything as we do it.

  Algo.:   Read in two XYZ files, and calculate...

  Output:  We output a differential displacement map in fig format
           for plotting and what-not.  Maybe at some point we'll figure
	   out way to output the numbers in an intelligent form...

	   Now with edge components!
	   If this switch is off, we do a simplified DD map; namely,
	   we put the atom between the *undisplaced* atoms, and point
	   it between them.
	   With this switch on, we do the true DD map, using the
	   *displaced* atom positions to calculate the connecting
	   vector.

	   A user might not want the TRUE one if there is some
	   unexpected weirdness in the two files.
*/

// ************************** COMPILIATION OPTIONS *********************

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "io.H"   // All of our "read in file", etc.
#include "drawfig.H"
#include "nnpair.H"

// ****************************** SUBROUTINES **************************

// Given a set of atoms, determines the ideal scale factor to maximize
// space on the page, ideal atom size (can be recalculated later), and
// portrait vs. landscape
void auto_scale (int N, double** r, 
		 double &a0, double &x0, double &y0, double &r0, 
		 int &portrait);

// Determine the 2d NN list for our disc.
void plane_nn_pair (int Natoms, double** r, double Rcut, int &NNpairs, 
		    nn_pair_type* &nn_pair_list, int** &nn_list);

// Construct the set of all right-handed triads:
void make_triads (int Natoms, int NNpairs, nn_pair_type* nn_pair_list, 
		  int** nn_list, int &Ntriad, int** &triad);


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<perfect-xtal.file> <dislocated-xtal.file> <Rcut> <Rmax>";

const int NFLAGS = 3;
const char USERFLAGLIST[NFLAGS] = {'e', 'n', 'b'};

const char* ARGEXPL = 
"  perfect-xtal:     XYZ file for perfect crystal\n\
  dislocated-xtal:  XYZ file with dislocation\n\
  Rcut:             cutoff for 2d NN search\n\
  Rmax:             radius around COM; if 0, not used.\n\
  -e     show edge components with direction and center changes\n\
  -n     write text amount of burgers vector for each pair\n\
  -b     write text for mini-burgers loops on triads";

int main ( int argc, char **argv ) 
{
  int i, j, k; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 65536; // 2^16, default.
  int FLAGON[NFLAGS]; // We use this to determine which flags are on.
  int EDGE_COMP;
  int NUMBERS;
  int TRIADS;

  char* args[NUMARGS];
  for (i=0; i<NFLAGS; ++i) FLAGON[i] = 0;

  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, FLAGON);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) )
      print_long_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // Okay, now we can actually read in our input file.
  FILE* perfect_file;
  FILE* disloc_file;

  int Natoms;
  double** pos;
  double** pos_d;
  double Rcut, Rmax;
  double z_thick;

  EDGE_COMP = FLAGON[0];
  NUMBERS = FLAGON[1];
  TRIADS = FLAGON[2];

  sscanf(args[2], "%lf", &Rcut);
  sscanf(args[3], "%lf", &Rmax);
  if (Rmax <= 0.0) Rmax = 1.e30;
  if (Rcut <= 0.0) {
    fprintf(stderr, "Rcut = %.5lf needs to be larger than 0.\n", Rcut);
    exit(-1);
  }  

  // Read each file!
  perfect_file = myopenr(args[0]);


  if (perfect_file == NULL) {
    fprintf(stderr, "File %s could not be opened.\n", args[0]);
    exit(ERROR_NOFILE);
  }

  char dump[512];
  // ==== Perfect file ====
  // N  # number of atoms.
  fgets(dump, sizeof(dump), perfect_file);
  sscanf(dump, "%d", &Natoms);
  // # comment line: has to contain z_thickness.
  fgets(dump, sizeof(dump), perfect_file);
  sscanf(dump, "%lf", &z_thick);
  if (Natoms <= 0) {
    fprintf(stderr, "Number of atoms less than 1?  N = %d\n", Natoms);
    ERROR = ERROR_BADFILE;
  }
  if (Natoms > MEMORY) {
    fprintf(stderr, "Number of atoms too big: N = %d; increase -m\n", Natoms);
    ERROR = ERROR_BADFILE;
  }

  // Do all the reading now!
  if (!ERROR) {
    char null[512];
    pos = new (double*)[Natoms];
    // ==== Perfect file ====
    for (i=0; i<Natoms; ++i) {
      pos[i] = new double[3];
      // <name> x y z
      fgets(dump, sizeof(dump), perfect_file);
      sscanf(dump, "%s %lf %lf %lf", null, pos[i], pos[i]+1, pos[i]+2);
    }
  }
  myclose(perfect_file);

  // ==== Dislocated file ====
  disloc_file = myopenr(args[1]);
  if (disloc_file == NULL) {
    fprintf(stderr, "File %s could not be opened.\n", args[1]);
    exit(ERROR_NOFILE);
  }
  fgets(dump, sizeof(dump), disloc_file);
  sscanf(dump, "%d", &i);
  // # comment line: has to contain z_thickness; check for compatibility.
  double ztest;
  fgets(dump, sizeof(dump), disloc_file);
  sscanf(dump, "%lf", &ztest);
  if (Natoms != i) {
    fprintf(stderr, "Number of atoms in files don't match?  %d != %d\n",
	    Natoms, i);
    ERROR = ERROR_BADFILE;
  }
  if (! dcomp(z_thick, ztest)) {
    fprintf(stderr, "thickness in undislocated = %.5lf, while in dislocated = %.5lf\n", z_thick, ztest);
    fprintf(stderr, "check to make sure files are compatible.\n");
    ERROR = ERROR_BADFILE;
  }
  if ( dcomp(z_thick, 0.) || (z_thick < 0.) ) {
    fprintf(stderr, "thickness = %.5lf is too small\n", z_thick);
    fprintf(stderr, "check to make sure that second line of XYZ file has the z thickness.\n");
    ERROR = ERROR_BADFILE;
  }

  // Do all the reading now!
  if (!ERROR) {
    char null[512];
    pos_d = new (double*)[Natoms];
    // ==== Dislocated file ====
    for (i=0; i<Natoms; ++i) {
      pos_d[i] = new double[3];
      // <name> x y z
      fgets(dump, sizeof(dump), disloc_file);
      sscanf(dump, "%s %lf %lf %lf", null, pos_d[i], pos_d[i]+1, pos_d[i]+2);
    }

    // Calc. COM shift, and set to 0.
    double COM[3];
    for (j=0; j<3; ++j) COM[j] = 0.;
    for (i=0; i<Natoms; ++i) for(j=0; j<3; ++j) COM[j] += pos[i][j];
    for (j=0; j<3; ++j) COM[j] /= Natoms;
    for (i=0; i<Natoms; ++i) for(j=0; j<3; ++j) pos[i][j] -= COM[j];

    for (j=0; j<3; ++j) COM[j] = 0.;
    for (i=0; i<Natoms; ++i) for(j=0; j<3; ++j) COM[j] += pos_d[i][j];
    for (j=0; j<3; ++j) COM[j] /= Natoms;
    for (i=0; i<Natoms; ++i) for(j=0; j<3; ++j) pos_d[i][j] -= COM[j];
  }

  myclose(disloc_file);

  // Bail now if we encountered some error along the way.
  if (ERROR) exit(ERROR);

  // Now, let's do a sweep through our list of atoms, and only keep
  // those in the maximum:
  double Rmax2;
  Rmax2 = Rmax*Rmax;
  MEMORY = Natoms; // Keep track of total amount allocated...
  Natoms = 0;
  for (i=0; i<MEMORY; ++i) {
    if ( (pos[i][0]*pos[i][0]+pos[i][1]*pos[i][1]) <= Rmax2 ) {
      // Keep this atom!
      if (i == Natoms)
	++Natoms; // No copying to do...
      else {
	// Shift atom i to position Natoms:
	for (j=0; j<3; ++j) {
	  pos[Natoms][j] = pos[i][j];
	  pos_d[Natoms][j] = pos_d[i][j];
	}
	++Natoms;
      }
    }
  }

  // *********************** DIFF DISP ANALYSIS **********************
  
  double* disp_z;
  disp_z = new double[Natoms];
  for (i=0; i<Natoms; ++i)
    disp_z[i] = pos_d[i][2] - pos[i][2];


  // *************************** NN ANALYSIS *************************
  // We need to analyze our disc to determine what all of the neighbors
  // are.  For this, we do a lot of weird things that are hidden away
  // in another routine, and for good reason.
  nn_pair_type* nn_pair_list;
  int NNpairs;
  int** nn_list = NULL;

  plane_nn_pair (Natoms, pos, Rcut, NNpairs, nn_pair_list, nn_list);

  // ****************************** OUTPUT ***************************
  // Autoscale that sucker!
  double a0, x0, y0, r_atom;
  int portrait;
  
  auto_scale(Natoms, pos, a0, x0, y0, r_atom, portrait);

  // Declare a figure drawing object.
  drawfig draw(stdout, portrait, a0, x0, y0);

  // Output all of the atoms:
  for (i=0; i<Natoms; ++i) {
    // Set fillstyle based on depth:
    draw.fillstyle(GREEN, (int)(WHITEFILL*insidecell(pos[i][2]/z_thick)) );
    draw.circle(pos[i][0], pos[i][1], r_atom);
  }
  

  // Now, let's do the differential displacements.
  // For now, this is going to do each vector *twice*, but we'll
  // deal with that later.
  // Run over all of the "bonds" in our list:
  draw.depth(draw.depth()+1);
  double zdisp;
  nn_pair_type* nn_pair;
  double burgers;
  double x, y;
  double vx, vy;

  burgers = z_thick;

  // Replace the "dscale" with the nn. dist scale (i.e., 1/3 b == nn dist)
  // Set text size:
  double dscale = 1./3.; // Scale factor for burgers vector
  
  draw.textstyle(FONT_COURIER, 12.);
  for (i=0; i<NNpairs; ++i) {
    nn_pair = nn_pair_list + i;
    // Only do the pairs where j > i (ensures that we don't double count):
    if (nn_pair->j > nn_pair->i) {
      zdisp = disp_z[nn_pair->j] - disp_z[nn_pair->i];
      
      // Make sure that it's between -1/2 burgers and 1/2 burgers:
      for ( ; (2.*zdisp) > (burgers+TOLER); zdisp -= burgers) ;
      for ( ; (2.*zdisp) <= -(burgers+TOLER); zdisp += burgers) ;
      // Now, scale zdisp by burgers vector
      zdisp *= 1./burgers;

      if (!EDGE_COMP) {
	x = 0.5*(pos[nn_pair->i][0] + pos[nn_pair->j][0]);
	y = 0.5*(pos[nn_pair->i][1] + pos[nn_pair->j][1]);
	vx = nn_pair->v_ij[0] * nn_pair->r;
	vy = nn_pair->v_ij[1] * nn_pair->r;
      }
      else {
	x = 0.5*(pos_d[nn_pair->i][0] + pos_d[nn_pair->j][0]);
	y = 0.5*(pos_d[nn_pair->i][1] + pos_d[nn_pair->j][1]);
	vx = pos_d[nn_pair->j][0] - pos_d[nn_pair->i][0];
	vy = pos_d[nn_pair->j][1] - pos_d[nn_pair->i][1];
      }
      // dscale determines what length of b is equal to the nn dist:
      draw.cvector(x, y, zdisp*vx/dscale, zdisp*vy/dscale);
      if (NUMBERS) {
	// Output a number there too.
	if (fabs(zdisp) >= 0.01) {
	  sprintf(dump, "%.0lf%%", fabs(zdisp)*100.);
	  draw.text(x, y, dump);
	}
      }
    }
  }

  draw.textstyle(FONT_COURIER, 14.);
  if (TRIADS) {
    int Ntriad;
    int** triad;
    make_triads (Natoms, NNpairs, nn_pair_list, nn_list, Ntriad, triad);
    // Loop through the triads, and calculate the Burgers loop :)
    double zd[3];
    int t;
    for (t=0; t<Ntriad; ++t) {
      // Go through our triad:
      for (j=0; j<3; ++j) {
	zd[j] = disp_z[triad[(j+1)%3][t]] - disp_z[triad[j][t]];
      
	// Make sure that it's between -1/2 burgers and 1/2 burgers:
	for ( ; (2.*zd[j]) > (burgers+TOLER); zd[j] -= burgers) ;
	for ( ; (2.*zd[j]) <= -(burgers+TOLER); zd[j] += burgers) ;
      }
      zdisp = (zd[0] + zd[1] + zd[2]) / burgers;
      i = triad[0][t];
      j = triad[1][t];
      k = triad[2][t];

      if ( fabs(zdisp) >= 0.01 ) {
	sprintf(dump, "%.0lf%%", zdisp*100.);
	draw.text( (pos[i][0]+pos[j][0]+pos[k][0])/3.,
		   (pos[i][1]+pos[j][1]+pos[k][1])/3.,
		   dump);
      }
    }
    
    // Garbage collection
    for (i=0; i<3; ++i)
      delete[] triad[i];
    delete[] triad;
  }

  // ************************* GARBAGE COLLECTION ********************
  delete[] disp_z;
  free_nn_list(Natoms, nn_list);
  delete[] nn_pair_list;

  for (i=0; i<MEMORY; ++i) { // Remember, we read more than we needed.
    delete[] pos[i];
    delete[] pos_d[i];
  }
  delete[] pos;
  delete[] pos_d;
    
  return 0;
}



// ============================= auto_scale ============================
// Given a set of atoms, determines the ideal scale factor to maximize
// space on the page, and portrait vs. landscape
// Assumes a 1200 dpi resolution, and a 1" margin on each side

const double DPI = 1200;
const double MARGIN_WIDTH = 6.5 * DPI;
const double MARGIN_HEIGHT = 9 * DPI;
const double R0_FRAC = 0.1;

const double MINSCALE = 1.e-7;
const double MAXSCALE = 1.e100;

void auto_scale (int N, double** r, 
		 double &a0, double &x0, double &y0, double &r0, 
		 int &portrait) 
{
  int i;
  double xmin, xmax, ymin, ymax;
  double *rp;
  double xs, ys;     // Scale factors in each direction
  double a0_p, a0_l; // Scale factor for portrait vs. landscape
  
  
  xmin = r[0][0]; xmax = r[0][0];
  ymin = r[0][1]; ymax = r[0][1];
  
  for (i=0; i<N; ++i) {
    rp = r[i];
    if (rp[0] < xmin) xmin = rp[0];
    if (rp[0] > xmax) xmax = rp[0];
    if (rp[1] < ymin) ymin = rp[1];
    if (rp[1] > ymax) ymax = rp[1];
  }
  
  x0 = 0.5*(xmax+xmin);
  y0 = 0.5*(ymax+ymin);

  // Portrait?
  if ( (xmax-xmin) < MINSCALE ) xs = MAXSCALE;
  else xs = MARGIN_WIDTH / (xmax-xmin);
  if ( (ymax-ymin) < MINSCALE ) ys = MAXSCALE;
  else ys = MARGIN_HEIGHT / (ymax-ymin);
  a0_p = ( xs < ys ) ? xs : ys;
  if (a0_p == MAXSCALE) a0_p = 1.;
  
  // Landscape?
  if ( (xmax-xmin) < MINSCALE ) xs = MAXSCALE;
  else xs = MARGIN_HEIGHT / (xmax-xmin);
  if ( (ymax-ymin) < MINSCALE ) ys = MAXSCALE;
  else ys = MARGIN_WIDTH / (ymax-ymin);
  a0_l = ( xs < ys ) ? xs : ys;
  if (a0_l == MAXSCALE) a0_l = 1.;
  
  if (a0_p > a0_l) {
    a0 = a0_p;
    portrait = -1;
  }
  else {
    a0 = a0_l;
    portrait = 0;
  }

  r0 = R0_FRAC * sqrt(MARGIN_HEIGHT*MARGIN_WIDTH / (M_PI * N) ) / a0;
}


// ============================ plane_nn_pair ==========================
// Determine the 2d NN list for our disc.
// There are a lot of bizarre little tricks to make this work, so pay
// attention.  At no time will my hands leave my wrists.
// Basically, we make a bounding box that's just a little too big,
// project out the z component, convert into "unit" coordinates, and
// feed to the NN engine.
// Most efficient method?  Probably not.
// Quick to code?  Relatively speaking, yes.
// My apologies in advance...
void plane_nn_pair (int N, double** r, double Rcut, int &NNpairs, 
		    nn_pair_type* &nn_pair_list, int** &nn_list) 
{
  int i;
  double xmin, xmax, ymin, ymax;
  double* rp;

  xmin = r[0][0]; xmax = r[0][0];
  ymin = r[0][1]; ymax = r[0][1];
  
  for (i=0; i<N; ++i) {
    rp = r[i];
    if (rp[0] < xmin) xmin = rp[0];
    if (rp[0] > xmax) xmax = rp[0];
    if (rp[1] < ymin) ymin = rp[1];
    if (rp[1] > ymax) ymax = rp[1];
  }

  // Now, make the "supercell"
  double Lx, Ly, Lz;
  double Lx1, Ly1;
  double cart[9];
  
  Lx = (xmax-xmin) + Rcut + 1.e-5;  Lx1 = 1./Lx;
  Ly = (ymax-ymin) + Rcut + 1.e-5;  Ly1 = 1./Ly;
  Lz =               Rcut + 1.e-5;

  for (i=0; i<9; ++i) cart[i] = 0.;
  cart[0] = Lx;
  cart[4] = Ly;
  cart[8] = Lz;

  // "Populate" it.
  double** u;
  u = new (double*)[N];
  for (i=0; i<N; ++i) {
    u[i] = new double[3];
    u[i][0] = insidecell( r[i][0] * Lx1 );
    u[i][1] = insidecell( r[i][1] * Ly1 );
    u[i][2] = 0.;
  }

  // Now, do the NN analysis!
  int Ngrid[3];      // Number of grid elements
  calc_grid(cart, Rcut, Ngrid);
  grid_elem_type* grid_list;
  make_grid(Ngrid, grid_list);
  populate_grid(Ngrid, grid_list, N, u);
  nn_grid(cart, Ngrid[0]*Ngrid[1]*Ngrid[2], grid_list, N, u, Rcut,
	  NNpairs, nn_pair_list);
  // Garbage collection
  free_grid(Ngrid, grid_list);
  // NOTE: we ASSUME that nn_list == NULL here!
  sort_nn_list(NNpairs, nn_pair_list, N, nn_list);

  // Garbage collection:
  for (i=0; i<N; ++i)
    delete[] u[i];
  delete[] u;    
}


// =============================== triads ==============================
// Construct the set of all right-handed triads:
void make_triads (int Natoms, int NNpairs, nn_pair_type* nn_pair_list, 
		  int** nn_list, int& Ntriad, int** &triad) 
{
  int i, j, k;
  int Napprox;
  
  Napprox = Natoms*6;
  triad = new (int*)[3];
  for (i=0; i<3; ++i) triad[i] = new (int)[Napprox];
  
  Ntriad = 0;
  int ni, nj, nk;
  nn_pair_type *pij, *pjk;
  int found;
  for (i=0; i<Natoms; ++i) {
    for (ni=1; ni<=nn_list[i][0]; ++ni) {
      pij = nn_pair_list + nn_list[i][ni];
      j = pij->j;
      if (i < j) {
	for (nj=1; nj<=nn_list[j][0]; ++nj) {
	  pjk = nn_pair_list + nn_list[j][nj];
	  k = pjk->j;
	  if (j < k) {
	    // Now, see if k has i for a neighbor:
	    found = 0;
	    for (nk=1; (nk<=nn_list[k][0]) && (!found); ++nk)
	      found = (i == nn_pair_list[nn_list[k][nk]].j);
	    if (found) {
	      // Add our triad
	      if (Ntriad >= Napprox) {
		fprintf(stderr, "Too many triads; set Rcut smaller.\n");
		return;
	      }
	      triad[0][Ntriad] = i;
	      // Make sure we're right handed.
	      if ( (pij->v_ij[0] * pjk->v_ij[1]) >
		   (pij->v_ij[1] * pjk->v_ij[0]) ) {
		triad[1][Ntriad] = j;
		triad[2][Ntriad] = k;
	      }
	      else {
		triad[1][Ntriad] = k;
		triad[2][Ntriad] = j;
	      }
	      ++Ntriad;
	    }
	  }
	}
      }
    }
  }
}
	 
