/*
  Program: anisotropic-xyz-ref.C
  Author:  D. Trinkle (slight modification by Anne Marie Tan)
  Date:    August 14, 2003 (modified August 21, 2015)
  Purpose: Calculate the anisotropic elastic solution for a general
           dislocation given:
	   1. dislocation line vector t = |t| (t1, t2, t3)
	   2. burgers vector          b =     (b1, b2, b3)
	   3. dislocation cut vector  m =     (m1, m2, m3) (normalized)
	   4. elastic constants Cmn, crystal class c

	   Fixed to use the correct slip plane definition (important
	   for edge and mixed dislocations); that is, the slip plane
	   vector normal is:

	     n0 = t x b

	   unless it's zero; then n0 = t x m0.  Note: m0 is to be 
	   perpendicular to t and in the plane of n0.

	   This new version reads in 2 XYZ files, where x=m, y=n, z=t.
       The first file is the undislocated geometry, the second file is the 
       "reference" geometry. The displacement field is evaluated according 
       to the "reference" geometry and then applied onto the undislocated 
       geometry, centered at (0,0). R_disloc = R_undisloc + u(R_reference).
 
       Ideally, you would want to run this code multiple times until self-
       consistency is achieved. The first time this code is called, 
       R_reference = R_undisloc. Subsequently, R_reference = R_disloc from
       the previous time. This is repeated until the new R_disloc = R_reference.

	   You need to be VERY CAREFUL to be consistent about
	   what information you feed this routine--it does next to no
	   checks on its own, so you can easily get nonsense out.

	   This code does NOT shift the origin in any way--it puts the
	   displacement field right at 0,0.  It needs a good "makeslab"
	   type code to construct the undislocated slab first, with the
	   appropriate center.  But this has the advantage that you can
	   stovepipe multiple anisotropic calls to, e.g., create a pair of
	   partials.

  Param.:  <cell> <infile> <undisloc> <reference> <outputstrainfile>
           cell:      cell file (see below for format)
           infile:    input file (see below for format)
	       undisloc:  undislocated crystal input XYZ file
           reference: input XYZ file to be used as reference for evaluating 
                      the displacement field (undislocated/dislocated crystal)
           outputstrainfile: file to output the strain tensor

	   ==== cell ====
       a0                               # Scale factor for unit cell
       a1.x a1.y a1.z                   # Cartesian coord of unit cell
       a2.x a2.y a2.z
       a3.x a3.y a3.z
       crystal-class <C_11> ... <C_66>  # Crystal class and elastic const.
       Natoms                           # Number of atoms in first unit cell
       u1.1 u1.2 u1.3                   # Atom locations, in direct coord.
       ...
       uN.1 uN.2 uN.3
	   ==== cell ====
	   
	   ==== infile ====
	   t1 t2 t3     # dislocation line direction (unit cell coord.)
	   b1 b2 b3 bd  # burgers vector (unit cell coord.)/bd
	   m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)
	   ==== infile ====

	   ==== undisloc ====
	   N               # standard xyz format
	   comment         # this *should* be the threading length
	   atomtype x y z
	   ...
	   ==== undisloc ====
 
       ==== reference ====
       N               # standard xyz format
       comment         # this *should* be the threading length
       atomtype x y z
       ...
       ==== reference ====
 
       ==== outputstrainfile ====
       N               # standard xyz format
       comment
       atomtype xx xy xz yx yy yz zx zy zz
       ...
       ==== outputstrainfile ====


  Flags:   MEMORY:  our setting for step size
	   VERBOSE: output the displacement fields too
	   TESTING: output practically everything as we do it.

  Algo.:   Read in everything, and just go *nuts* a calculatin'.

           First, we make sure that m0 is perp. to t and to b, and 
	   normalized.  We also construct n0 = t x m0.

	   We then define the vectors m(theta) and n(theta) as:

	   m(theta) =  m0*cos(theta) + n0*sin(theta)
	   n(theta) = -m0*sin(theta) + n0*cos(theta)

	   and the matrices (ab)_ij as

	   (ab)_ij = sum  a_k C_ikjl b_l
	              kl

           We have to do four integrals and store two of them as
	   functions of theta, namely, the two constant matrices:

	            1   Pi      -1
	   S_ij = - -  Int  (nn)  (nm)   dtheta
	            Pi  0       ik    kj

	            1     Pi                    -1
	   B_ij = -----  Int  (mm)   - (mn) (nn)  (nm)   dtheta
	          4Pi^2   0       ij      ik    kl    lj

	     (Note: B_ij = B_ji)

  	   and the two matrices as a function of theta:

	                    theta     -1
	   N_ij(theta) = 4Pi Int  (nn)   dtheta
	                      0       ij

	                theta     -1
	   L_ij(theta) = Int  (nn)  (nm)   dtheta
	                  0       ik    kj

	     (Note: N_ij = N_ji)

	   Also, due to the theta periodicity, we only have to evaluate
	   these from 0..Pi, since 

	     N_ij(theta+Pi) = N_ij(theta) + N_ij(Pi)

	   and similarly for L_ij.

	   Also, S_ij = -1/Pi * L_ij(Pi), so we have only three integrals
	   to do.

	   *Then*, to turn these all into our displacement using the 
	   equation:

	   u (|x|, theta) = [-S  ln |x| + N  B   + L  S  ] b  / 2 Pi
	    i                  is          ik ks    ik ks   s

	   Voila! (whew)

	   We integrate using a stepping scheme based on Simpson's
	   extended rule; basically, to integrate f(x) from 0 to x,
	   we calculate f(x) at a grid of points, and if F(N-1) is
	   the integral to x-h, and f[i] = f(x-ih), then:

   	     F(N) = F(N-1) + h*SUM(i=0..3, int_weight[i]*f[i])

	   which works very well.  To get the first two points, we 
	   actually have to evaluate f at -h and -2h, and then
	   start with F(0) = 0.  It works (go figure).

	   We do 16384 integration steps (2^14... woohoo!) to make
	   sure that we have something reasonable :)

  Output:  If we're verbose, we'll output the theta dependence of u_i, and
           also the ln |x| prefactor.

	   For fun, and profit, we can output the energy prefactor:

	   E = b_i B_ij b_j  : self-energy prefactor per length

	   The real meat of the code, though, is to actually put these
	   displacements to work by outputting the XYZ files for a
	   cylindrical slab material.  We do this by adding in the 
	   displacements... not too hard.
*/

// ************************** COMPILIATION OPTIONS ***********************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include "dcomp.H"
#include "io-short.H"
#include "matrix.H"
#include "elastic.H"
#include "cell.H"
#include "integrate.H"

// This is the permutation matrix; eps[i][j][k] =
//  1: if ijk is an even permutation of (012)
// -1: if ijk is an odd permutation of (012)
//  0: otherwise
const int eps[3][3][3] = {
  {{0,0,0}, {0,0,1}, {0,-1,0}},
  {{0,0,-1}, {0,0,0}, {1,0,0}},
  {{0,1,0}, {-1,0,0}, {0,0,0}}
};

// ****************************** SUBROUTINES ****************************

inline double dot(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


void m_theta(double theta, double m[3], double n[3], double mt[3]) 
{
  mt[0] = cos(theta)*m[0] + sin(theta)*n[0];
  mt[1] = cos(theta)*m[1] + sin(theta)*n[1];
  mt[2] = cos(theta)*m[2] + sin(theta)*n[2];
}

void n_theta(double theta, double m[3], double n[3], double nt[3]) 
{
  nt[0] = -sin(theta)*m[0] + cos(theta)*n[0];
  nt[1] = -sin(theta)*m[1] + cos(theta)*n[1];
  nt[2] = -sin(theta)*m[2] + cos(theta)*n[2];
}

void a_mult_b (double a[3], double b[3], double Cijkl[9][9],
	       double ab[9]) 
{
  int i, j, k, l;
  for (i=0; i<3; ++i)
    for (j=0; j<3; ++j) {
      ab[index(i,j)] = 0.;
      for (k=0; k<3; ++k)
	for (l=0; l<3; ++l)
	  ab[index(i,j)] += a[k]*Cijkl[index(k,i)][index(j,l)]*b[l];
    }
}

void a_mult_a (double a[3], double Cijkl[9][9], double aa[9]) 
{
  int i, j, k, l;
  for (i=0; i<3; ++i) {
    for (j=0; j<i; ++j)
      aa[index(i,j)] = aa[index(j,i)];
    for (   ; j<3; ++j) {
      aa[index(i,j)] = 0.;
      for (k=0; k<3; ++k)
	for (l=0; l<3; ++l)
	  aa[index(i,j)] += a[k]*Cijkl[index(k,i)][index(j,l)]*a[l];
    }
  }
}

void print_mat (double a[9]) 
{
  int i, j;
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j)
      printf("  %10.5lf", a[index(i,j)]);
    printf("\n");
  }
}

void calcstrain(double x, double theta, double m0[3], double n0[3], double Cijkl[9][9], double Sint[9], double b0[3], double Bint[9], double strain[9])
{
    double mt[3], nt[3], nnt[9], nmt[9], nni[9], nninm[9];
    
    //calc theta pieces
    m_theta(theta, m0, n0, mt);
    n_theta(theta, m0, n0, nt);
    a_mult_a(nt, Cijkl, nnt);
    a_mult_b(nt, mt, Cijkl, nmt);
    double detnn = 1./inverse(nnt, nni);
    for (int i=0; i<9; ++i) nni[i] *= detnn;
    mult(nni, nmt, nninm);
    
    //strain
    double Sb[3], NBb[3], NB[9], LS[9], sum[9],strain_tot[9];
    mult_vect(Sint, b0, Sb);
    mult(nni, Bint, NB);
    mult(nninm, Sint, LS);
    for (int i=0; i<9; ++i) sum[i] = 4.0*M_PI*NB[i] + LS[i];
    mult_vect(sum, b0, NBb);
    
    for(int idx = 0; idx < 3; ++idx) {
        for(int jdx = 0; jdx < 3; ++jdx) {
            strain_tot[idx + 3*jdx] = 0.5*M_1_PI*(-mt[jdx]*Sb[idx] + nt[jdx]*NBb[idx])/x;
        }
    }
    // symmetrize strain
    for(int idx = 0; idx < 3; ++idx) {
        for(int jdx = 0; jdx < 3; ++jdx) {
            strain[idx + 3*jdx] = 0.5*(strain_tot[idx + 3*jdx] + strain_tot[jdx + 3*idx]);
        }
    }
    
}


void calcstrain_fd(double mcoord, double ncoord, double aln, double inv_dtheta, double xyz0[3], double** u_xyz, double strain[9])
{
    // Consider a point slightly in the +m direction:
    double dist_posdisp_m = sqrt((mcoord+0.0001)*(mcoord+0.0001) + ncoord*ncoord);
    double theta_posdisp_m = atan2(ncoord, (mcoord+0.0001));
    if (theta_posdisp_m < 0.) theta_posdisp_m += (2.*M_PI);
    
    // Compute displacements at this point:
    // xyz0*(ln|x| - ln(a0)) + u_xyz(theta)
    double lnr_posdisp_m = log(dist_posdisp_m) + aln;
    // Now, linearly interpolate for theta:
    double kreal_posdisp_m = theta_posdisp_m * inv_dtheta;
    int k_posdisp_m = (int) kreal_posdisp_m;
    double alpha_posdisp_m = kreal_posdisp_m - k_posdisp_m, beta_posdisp_m = 1. - alpha_posdisp_m;
    
    double disp_posdisp_m[3];
    for (int d=0; d<3; ++d) {
        disp_posdisp_m[d] = xyz0[d]*lnr_posdisp_m + beta_posdisp_m*u_xyz[k_posdisp_m][d] + alpha_posdisp_m*u_xyz[k_posdisp_m+1][d];
    }
    
    // Consider a point slightly in the -m direction:
    double dist_negdisp_m = sqrt((mcoord-0.0001)*(mcoord-0.0001) + ncoord*ncoord);
    double theta_negdisp_m = atan2(ncoord, (mcoord-0.0001));
    if (theta_negdisp_m < 0.) theta_negdisp_m += (2.*M_PI);
    
    // Compute displacements at this point:
    // xyz0*(ln|x| - ln(a0)) + u_xyz(theta)
    double lnr_negdisp_m = log(dist_negdisp_m) + aln;
    // Now, linearly interpolate for theta:
    double kreal_negdisp_m = theta_negdisp_m * inv_dtheta;
    int k_negdisp_m = (int) kreal_negdisp_m;
    double alpha_negdisp_m = kreal_negdisp_m - k_negdisp_m, beta_negdisp_m = 1. - alpha_negdisp_m;

    double disp_negdisp_m[3];
    for (int d=0; d<3; ++d) {
        disp_negdisp_m[d] = xyz0[d]*lnr_negdisp_m + beta_negdisp_m*u_xyz[k_negdisp_m][d] + alpha_negdisp_m*u_xyz[k_negdisp_m+1][d];
    }

    // Consider a point slightly in the +n direction:
    double dist_posdisp_n = sqrt(mcoord*mcoord + (ncoord+0.0001)*(ncoord+0.0001));
    double theta_posdisp_n = atan2((ncoord+0.0001), mcoord);
    if (theta_posdisp_n  < 0.) theta_posdisp_n += (2.*M_PI);
    
    // Compute displacements at this point:
    // xyz0*(ln|x| - ln(a0)) + u_xyz(theta)
    double lnr_posdisp_n = log(dist_posdisp_n) + aln;
    // Now, linearly interpolate for theta:
    double kreal_posdisp_n = theta_posdisp_n * inv_dtheta;
    int k_posdisp_n = (int) kreal_posdisp_n;
    double alpha_posdisp_n = kreal_posdisp_n - k_posdisp_n, beta_posdisp_n = 1. - alpha_posdisp_n;
    
    double disp_posdisp_n[3];
    for (int d=0; d<3; ++d) {
        disp_posdisp_n[d] = xyz0[d]*lnr_posdisp_n + beta_posdisp_n*u_xyz[k_posdisp_n][d] + alpha_posdisp_n*u_xyz[k_posdisp_n+1][d];
    }
    
    // Consider a point slightly in the -n direction:
    double dist_negdisp_n = sqrt(mcoord*mcoord + (ncoord-0.0001)*(ncoord-0.0001));
    double theta_negdisp_n = atan2((ncoord-0.0001), mcoord);
    if (theta_negdisp_n  < 0.) theta_negdisp_n += (2.*M_PI);

    // Compute displacements at this point:
    // xyz0*(ln|x| - ln(a0)) + u_xyz(theta)
    double lnr_negdisp_n = log(dist_negdisp_n) + aln;
    // Now, linearly interpolate for theta:
    double kreal_negdisp_n = theta_negdisp_n * inv_dtheta;
    int k_negdisp_n = (int) kreal_negdisp_n;
    double alpha_negdisp_n = kreal_negdisp_n - k_negdisp_n, beta_negdisp_n = 1. - alpha_negdisp_n;

    double disp_negdisp_n[3];
    for (int d=0; d<3; ++d) {
        disp_negdisp_n[d] = xyz0[d]*lnr_negdisp_n + beta_negdisp_n*u_xyz[k_negdisp_n][d] + alpha_negdisp_n*u_xyz[k_negdisp_n+1][d];
    }
    
    // compute strains using finite differences (central differences)
    strain[0] = (disp_posdisp_m[0] - disp_negdisp_m[0])/0.0002; //eps_mm
    strain[1] = 0.5*((disp_posdisp_m[1] - disp_negdisp_m[1])/0.0002 + (disp_posdisp_n[0] - disp_negdisp_n[0])/0.0002); //eps_mn
    strain[2] = 0.5*((disp_posdisp_m[2] - disp_negdisp_m[2])/0.0002 + 0.0); //eps_mt -> should be 0
    strain[3] = strain[1]; //eps_nm
    strain[4] = (disp_posdisp_n[1] - disp_negdisp_n[1])/0.0002; //eps_nn
    strain[5] = 0.5*((disp_posdisp_n[2] - disp_negdisp_n[2])/0.0002 + 0.0); //eps_nt -> should be 0
    strain[6] = strain[2]; //eps_tm -> should be 0
    strain[7] = strain[5]; //eps_tn -> should be 0
    strain[8] = 0.0; //eps_tt -> should be 0
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "[-hvt] [-s STEPS] cell infile undisloc reference";

const char* ARGEXPL =
" cell:      cell file (-h for format)\n\
  infile:    input file (-h for format)\n\
  undisloc:  undislocated crystal input XYZ file\n\
  reference: reference crystal input XYZ file\n\
\n\
  -s STEPS  number of integration steps\n\
  -v        verbosity\n\
  -t        testing\n\
  -h        help";

const char* FILEEXPL =
"==== cell ====\n\
a0                              # Scale factor for unit cell\n\
a1.x a1.y a1.z                  # Cartesian coord of unit cell\n\
a2.x a2.y a2.z\n\
a3.x a3.y a3.z\n\
crystal-class <C_11> ... <C_66> # Crystal class and elastic const.\n\
Natoms                          # Number of atoms in first unit cell\n\
u1.1 u1.2 u1.3                  # Atom locations, in direct coord.\n\
...\n\
uN.1 uN.2 uN.3\n\
==== cell ====\n\
\n\
==== infile ====\n\
t1 t2 t3     # dislocation line direction (unit cell coord.)\n\
b1 b2 b3 bd  # burgers vector (unit cell coord.)/bd \n\
m1 m2 m3     # dislocation cut vector (perp. to t, in slip plane)\n\
==== infile ====\n\
\n\
==== undisloc ====\n\
N               # standard xyz format\n\
comment         # this *should* be the threading length\n\
atomtype x y z\n\
...\n\
==== undisloc ====\n\
\n\
==== reference ====\n\
N               # standard xyz format\n\
comment         # this *should* be the threading length\n\
atomtype x y z\n\
...\n\
==== reference ====\n";

int main ( int argc, char **argv ) 
{
  int i, j, k; // General counting variables.

  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int Nsteps = 16384; // 2^14, default

  char ch;
  while ((ch = getopt(argc, argv, "vths:")) != -1) {
    switch (ch) {
    case 's':
      Nsteps = (int)strtol(optarg, (char**)NULL, 10);
      break;
    case 'v':
      VERBOSE = 1;
      break;
    case 't':
      TESTING = 1;
      VERBOSE = 1;
      break;
    case 'h':
    case '?':
    default:
      ERROR = 1;
    }
  }
     
  argc -= optind; if (argc<NUMARGS && !ERROR) ERROR = 2;
  argv += optind;

  if (TESTING) {
    printf("# Nsteps=%d\n", Nsteps);
  }
  // We're going to use the number of steps according to our preferred
  // amount of memory allocation.
  if (Nsteps < 4) {
    fprintf(stderr, "Nsteps (%d) must be 4 or larger.\n", Nsteps);
    ERROR = 2;
  }  

  // All hell broken loose yet?
  if (ERROR != 0) {
    fprintf(stderr, "%s %s\n%s\n", progname, ARGLIST, ARGEXPL);
    if (ERROR == 1) {
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "Crystal classes:\n%s\n", CRYSTAL_CLASS);
      fprintf(stderr, "\nElastic constants ordering:\n");
      for (k=0; k<NCLASSES; ++k) {
	    fprintf(stderr, "  Class %2d (%d):", k, class_len[k]);
	    for (i=0; i<class_len[k]; ++i) {
	      fprintf(stderr, " C_%2d", class_Cij[k][i]);
        }
	    fprintf(stderr, "\n");
      }
    }
    exit(ERROR);
  }


  // ****************************** INPUT ****************************
  char dump[512];
  char *cell_name = argv[0];
  char *infile_name = argv[1];
  char *undisloc_name = argv[2];
  char *reference_name = argv[3];
  char *strainfile_name = argv[4];
  FILE* infile;
  FILE* infile_ref;

  double cart[9];
  int crystal; // crystal class
  double* Cmn_list; // elastic constant input
  double Cijkl[9][9];

  // disl. line, burgers vect, cut, center of dislocation (all in unit coord)
  int tu0[3], bu0[3], mu0[3];		// all in unit cell coord; must be int.
  int bu_denom;				// denominator for burgers vector (partials)
  double t0[3], b0[3], m0[3], n0[3];	// n0 = t0 x m0, in cart. coord. 

  // First, read in the cell.
  infile = myopenr(cell_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", cell_name);
    exit(1);
  }
  {
    int Natoms=NO_ATOMS;
    double **u_atoms=NULL;
    ERROR = read_cell(infile, cart, crystal, Cmn_list, u_atoms, Natoms);
  }
  myclose(infile);
  
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_ZEROVOL) ) 
      fprintf(stderr, "Cell had zero volume.\n");
    if ( has_error(ERROR, ERROR_LEFTHANDED) )
      fprintf(stderr, "Left-handed cell.\n");
    exit(ERROR);
  }

  if (TESTING)
    verbose_output_cell(cart, crystal, Cmn_list, NULL, 0);

  // Now, read in the dislocation information
  infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", infile_name);
    exit(1);
  }

  // **** NOTE: all input in unit cell coord, so first three vect. are int.
  //  t1 t2 t3            # dislocation line
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d", &tu0[0], &tu0[1], &tu0[2]);

  //  b1 b2 b3            # burgers vector
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d %d", &bu0[0], &bu0[1], &bu0[2], &bu_denom);
  // For backwards compatibility...
  if (bu_denom == 0) bu_denom = 1;

  //  m1 m2 m3            # dislocation cut vector (perp. to t)
  nextnoncomment(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d", &mu0[0], &mu0[1], &mu0[2]);

  myclose(infile);

  // Now, convert vectors from unit cell to cartesian coord.:
  mult_vect(cart, tu0, t0);
  mult_vect(cart, bu0, b0); for (i=0; i<3; ++i) b0[i] *= 1./bu_denom;
  mult_vect(cart, mu0, m0);
  // Sanity check on vectors:
  if ( dot(t0, t0) < 1e-8 ) {
    fprintf(stderr, "Bad t vector.\n");
    ERROR = ERROR_BADFILE;
  }
  if ( dot(b0, b0) < 1e-8 ) {
    fprintf(stderr, "Bad b vector.\n");
    ERROR = ERROR_BADFILE;
  }
  // We also need to project out any t components of m, and place
  // it in the slip plane (provided t x b isn't 0):
  for (i=0; i<3; ++i)
    m0[i] -= dot(m0, t0)/dot(t0,t0) * t0[i];
  // Now, calculate n0 (we'll recalc it later, correctly)
  for (i=0; i<3; ++i) {
    n0[i] = 0.;
    for (j=0; j<3; ++j)
      for (k=0; k<3; ++k)
	n0[i] += eps[i][j][k]*t0[j]*b0[k];
  }
  if (! dcomp(dot(n0,n0), 0.) ) 
    // We have a non-screw dislocation...
    for (i=0; i<3; ++i)
      m0[i] -= dot(m0, n0)/dot(n0,n0) * n0[i];

  if ( dcomp(dot(m0, m0), 0.) ) {
    fprintf(stderr, "Bad m0 vector (parallel to t or out of the t x b slip plane).\n");
    ERROR = ERROR_BADFILE;
  }

  // Now, normalize:
  double magn;
  magn = 1./sqrt(dot(m0,m0));
  for (i=0; i<3; ++i) m0[i] *= magn;

  
  if (VERBOSE) {
    printf("# Run dislocation along (%.5lf %.5lf %.5lf)\n",t0[0],t0[1],t0[2]);
    printf("# Burgers vector        (%.5lf %.5lf %.5lf), magn = %.5lf\n", 
	   b0[0],b0[1],b0[2], sqrt(dot(b0,b0)));
    printf("# Cut direction         (%.5lf %.5lf %.5lf)\n",m0[0],m0[1],m0[2]);
  }

  // Calculate elastic constant matrix:

  make_Cijkl(crystal, Cmn_list, Cijkl);

  // ***************************** ANALYSIS **************************

  if (VERBOSE) {
    double comp;
    comp = fabs(dot(b0,t0)/sqrt(dot(b0,b0)*dot(t0,t0)));
    printf("# Screw component: %5.2lf%%  Edge component: %5.2lf%%\n",
	   comp*100.0, (1.-comp)*100.0);
  }
  
  // Now, compute n0 = t0 x m0:
  for (i=0; i<3; ++i) {
    n0[i] = 0;
    for (j=0; j<3; ++j)
      for (k=0; k<3; ++k)
	n0[i] += eps[i][j][k] * t0[j] * m0[k];
  }
  // Normalize:
  magn = 1./sqrt(dot(n0,n0));
  for (i=0; i<3; ++i) n0[i] *= magn;
  
  if (TESTING) {
    printf("##\n## Normalized vectors:\n");
    printf("## Run dislocation along (%.5lf %.5lf %.5lf)\n", t0[0],t0[1],t0[2]);
    printf("## Cut direction         (%.5lf %.5lf %.5lf)\n", m0[0],m0[1],m0[2]);
    printf("## Perp direction        (%.5lf %.5lf %.5lf)\n", n0[0],n0[1],n0[2]);
  }
  if (VERBOSE) {
     printf("# %17.12lf %17.12lf %17.12lf : normalized x axis\n", m0[0], m0[1], m0[2]);
     printf("# %17.12lf %17.12lf %17.12lf : normalized y axis\n", n0[0], n0[1], n0[2]);
     printf("# %17.12lf %17.12lf %17.12lf : normalized z axis\n",
            t0[0]/sqrt(dot(t0,t0)), t0[1]/sqrt(dot(t0,t0)), t0[2]/sqrt(dot(t0,t0))); 
  }


  // Now some evaluating of integrals :)
  double theta;
  double dtheta;
  dtheta = M_PI / Nsteps;
  double mt[3], nt[3];
  double nnt[9], mmt[9], nmt[9], mnt[9], nnti[9];
  double detnn;

  // We have to integrate three functions.
  double **Nint, **Lint;
  double Bint[9];

  Nint = new double*[Nsteps+1];
  Lint = new double*[Nsteps+1];
  for (i=0; i<=Nsteps; ++i) {
    Nint[i] = new double[9];
    Lint[i] = new double[9];
  }
  
  // Function evaluations, stored for integration purposes.
  double nn_old[4][9], nnnm_old[4][9], mnnnnm_old[4][9];

  // First, prime the integration pump:
  for (k=1; k<=3; ++k) {
    theta = -(k-1)*dtheta;
    // Eval (nn), (nm), (mn), (mm), and (nn)^-1
    m_theta(theta, m0, n0, mt);
    n_theta(theta, m0, n0, nt);
    a_mult_a(mt, Cijkl, mmt);
    a_mult_a(nt, Cijkl, nnt);
    a_mult_b(nt, mt, Cijkl, nmt);
    transpose(nmt, mnt);
    detnn = 1./inverse(nnt, nnti);
    for (i=0; i<9; ++i) nnti[i] *= detnn;
    // Now, put into the function evaluations:
    for (i=0; i<9; ++i) nn_old[k][i] = nnti[i];
    mult(nnti, nmt, nnnm_old[k]);
    mult(mnt, nnnm_old[k], mnnnnm_old[k]);
    for (i=0; i<9; ++i)
      mnnnnm_old[k][i] = mmt[i] - mnnnnm_old[k][i];
    // And HERE's where we'd integrate, if we wanted to... :)
  }

  // Now we've got EVERYTHING, let's integrate!
  // theta = 0 is easy...
  for (i=0; i<9; ++i) {
    Nint[0][i] = 0.;
    Lint[0][i] = 0.;
    Bint[i] = 0.;
  }

  for (k=1; k<=Nsteps; ++k) {
    theta = k*dtheta;
    // Eval (nn), (nm), (mn), (mm), and (nn)^-1
    m_theta(theta, m0, n0, mt);
    n_theta(theta, m0, n0, nt);
    a_mult_a(mt, Cijkl, mmt);
    a_mult_a(nt, Cijkl, nnt);
    a_mult_b(nt, mt, Cijkl, nmt);
    transpose(nmt, mnt);
    detnn = 1./inverse(nnt, nnti);
    for (i=0; i<9; ++i) nnti[i] *= detnn;

    // Now, put into the function evaluations:
    for (i=0; i<9; ++i) nn_old[0][i] = nnti[i];
    mult(nnti, nmt, nnnm_old[0]);
    mult(mnt, nnnm_old[0], mnnnnm_old[0]);
    for (i=0; i<9; ++i)
      mnnnnm_old[0][i] = mmt[i] - mnnnnm_old[0][i];
    // Now, we can integrate!
    for (i=0; i<9; ++i) {
      Nint[k][i] = Nint[k-1][i];
      Lint[k][i] = Lint[k-1][i];
      for (j=0; j<4; ++j) {
	Nint[k][i] += dtheta*int_weight[j]*nn_old[j][i];
	Lint[k][i] += dtheta*int_weight[j]*nnnm_old[j][i];
	Bint[i]    += dtheta*int_weight[j]*mnnnnm_old[j][i];
      }
    }
    // Now, we slide down all of our "old" values:
    for (j=3; j>0; --j)
      for (i=0; i<9; ++i) {
	nn_old[j][i] = nn_old[j-1][i];
	nnnm_old[j][i] = nnnm_old[j-1][i];
	mnnnnm_old[j][i] = mnnnnm_old[j-1][i];
      }
    // And do it all again!
  }
  // Finally, define S, and scale everything appropriately.
  double Sint[9];

  for (k=0; k<=Nsteps; ++k)
    for (i=0; i<9; ++i) 
      Nint[k][i] *= (4.*M_PI);
  
  for (i=0; i<9; ++i) {
    Sint[i] = -Lint[Nsteps][i] * M_1_PI;
    Bint[i] *= 0.25*M_1_PI*M_1_PI;
  }      

  // Displacement!
  double** u;
  double NB[9], LS[9], sum[9];

  u = new double*[Nsteps+1];
  for (k=0; k<=Nsteps; ++k) {
    u[k] = new double[3];
    theta = k*dtheta;
    // Eval. the theta part of u_i:
    mult(Nint[k], Bint, NB);
    mult(Lint[k], Sint, LS);
    for (i=0; i<9; ++i) sum[i] = NB[i] + LS[i];
    mult_vect(sum, b0, u[k]);
    for (i=0; i<3; ++i) u[k][i] *= 0.5*M_1_PI;
  }

  // Now, let's put those displacements into cylindrical coordinates:
  double** u_xyz;
  double tmagn;
  tmagn = 1./sqrt(dot(t0,t0));
  u_xyz = new double*[2*Nsteps+1];
  for (k=0; k<=Nsteps; ++k) {
    u_xyz[k] = new double[3];
    u_xyz[k][0] = dot(u[k], m0);
    u_xyz[k][1] = dot(u[k], n0);
    u_xyz[k][2] = dot(u[k], t0) * tmagn;
  }
  for ( ; k<=(2*Nsteps); ++k) {
    u_xyz[k] = new double[3];
    u_xyz[k][0] = dot(u[k-Nsteps], m0) + u_xyz[Nsteps][0];
    u_xyz[k][1] = dot(u[k-Nsteps], n0) + u_xyz[Nsteps][1];
    u_xyz[k][2] = dot(u[k-Nsteps], t0) * tmagn + u_xyz[Nsteps][2];
  }

    
  
  // ****************************** OUTPUT ***************************

  // Human readable (sorta) first:

  if (VERBOSE) {
    // Let's give the energy per-length prefactor:
    double u0[3];
    double tnorm[3];
    for (i=0; i<3; ++i) tnorm[i] = t0[i] *tmagn;

    mult_vect(Bint, b0, u0);
    printf("# Energy per unit length prefector = %.15lf\n", dot(b0, u0));
    
    if (TESTING) {
      // First, let's dump out the radial part (ln |x| prefactor):
      mult_vect(Sint, b0, u0);
      for (i=0; i<3; ++i) u0[i] *= -0.5*M_1_PI;
      printf("# radial prefactor: u.t, u.m(0), u.n(0) =\n");
      printf("# %.15lf %.15lf %.15lf\n", 
	     dot(u0, tnorm), dot(u0, m0), dot(u0, n0));
      printf("# \n");
      
      // Now, let's output it; next, just the angular part.
      printf("# theta  u.t  u.m(theta)  u.n(theta)\n");
      // 0..Pi
      for (k=0; k<=Nsteps; ++k) {
	theta = k*dtheta;
	// For output in "dislocation coordinates":
	m_theta(theta, m0, n0, mt);
	n_theta(theta, m0, n0, nt);
	printf("%10.7lf %.15lf %.15lf %.15lf\n", theta,
	       dot(u[k], tnorm), dot(u[k], mt), dot(u[k], nt));
      }
      // Pi .. 2Pi
      // We handle this simply adding in the u0 = u[Nsteps]
      for (i=0; i<3; ++i) u0[i] = u[Nsteps][i];
      for (k=1; k<=Nsteps; ++k) {
	theta = k*dtheta + M_PI;
	// For output in "dislocation coordinates":
	m_theta(theta, m0, n0, mt);
	n_theta(theta, m0, n0, nt);
	printf("%10.7lf %.15lf %.15lf %.15lf\n", theta,
	       dot(u[k], tnorm)+dot(u0,tnorm), 
	       dot(u[k], mt)+dot(u0,mt),
	       dot(u[k], nt)+dot(u0,nt));
      }
    }
  }

  if (ERROR) {
    fprintf(stderr, "An error occured, and we're getting out now.\n");
  }
  else {
    // Now, let's treat the logarithmic part:
    double u0[3], xyz0[3];
    mult_vect(Sint, b0, u0);
    for (i=0; i<3; ++i) u0[i] *= -0.5*M_1_PI;
    xyz0[0] = dot(u0, m0);
    xyz0[1] = dot(u0, n0);
    xyz0[2] = dot(u0, t0)*tmagn;
    // Our scaling factor:
    double aln;
    // double a0;
    // a0 = exp( log(det(cart)/Natoms) / 3.);
    // aln = -log(a0);
    aln = - log(det(cart)) / 3.;
    // for interpolation purposes:
    double inv_dtheta = 1./dtheta;

    // Output XYZ files!!
    infile = myopenr(undisloc_name);
    infile_ref = myopenr(reference_name);
    int Nslab;
    // Natoms
    nextnoncomment(dump, sizeof(dump), infile);
    printf("%s\n", dump);
    sscanf(dump, "%d", &Nslab);
    nextnoncomment(dump, sizeof(dump), infile_ref); // dummy readline
    // comment
    nextnoncomment(dump, sizeof(dump), infile);
    printf("%s\n", dump);
    nextnoncomment(dump, sizeof(dump), infile_ref); // dummy readline
      
    FILE *strainfile = myopenw(strainfile_name);
    fprintf(strainfile, "%d\n", Nslab);
    fprintf(strainfile, "%s", dump);
      
    for (int n=0; n<Nslab; ++n) {
      char atomname[512];
      double xyz[3];
      double xyz_ref[3];
      // undislocated atom x y z
      nextnoncomment(dump, sizeof(dump), infile);
      sscanf(dump, "%s %lf %lf %lf", atomname, xyz, xyz+1, xyz+2);
      // reference atom x y z
      nextnoncomment(dump, sizeof(dump), infile_ref);
      sscanf(dump, "%*s %lf %lf %lf", xyz_ref, xyz_ref+1, xyz_ref+2);

      // Now, we need to do some analysis on our displacements; first,
      // we need to calculate the distance from the dislocation,
      // and the magical angle theta for each:
      double dist_ref = sqrt(xyz_ref[0]*xyz_ref[0] + xyz_ref[1]*xyz_ref[1]);
      ERROR = dcomp(dist_ref, 0.);
      double theta_ref = atan2(xyz_ref[1], xyz_ref[0]);
      if (theta_ref < 0.) theta_ref += (2.*M_PI);
      if (ERROR) {
	    fprintf(stderr, "You managed to center your dislocation right on an atom... that's not so good.\n");
	    break;
      }
      
      // Let's displace all of the atoms accordingly:
      // xyz0*(ln|x| - ln(a0)) + u_xyz(theta)
      double lnr = log(dist_ref) + aln;
      // Now, linearly interpolate for theta:
      double kreal = theta_ref * inv_dtheta;
      int k = (int) kreal;
      double alpha = kreal - k, beta = 1. - alpha;
      //double disp[3];
      for (int d=0; d<3; ++d) {
          xyz[d] += xyz0[d]*lnr + beta*u_xyz[k][d] + alpha*u_xyz[k+1][d];
      }
      // output
      printf("%s %20.15lf %20.15lf %20.15lf\n", atomname, xyz[0], xyz[1], xyz[2]);
        
      double strain_xyz[9];
      //alcstrain(dist_ref, theta_ref, m0, n0, Cijkl, Sint, b0, Bint, strain_xyz);
      calcstrain_fd(xyz_ref[0], xyz_ref[1], aln, inv_dtheta, xyz0, u_xyz, strain_xyz);

      fprintf(strainfile, "%s % 20.15lf % 20.15lf % 20.15lf % 20.15lf % 20.15lf % 20.15lf % 20.15lf % 20.15lf % 20.15lf %20.15lf\n", atomname, strain_xyz[0], strain_xyz[1], strain_xyz[2], strain_xyz[3], strain_xyz[4], strain_xyz[5], strain_xyz[6], strain_xyz[7], strain_xyz[8], dist_ref);
    }
    myclose(infile);
    myclose(infile_ref);
  }

  // ************************* GARBAGE COLLECTION ********************
  for (i=0; i<=(2*Nsteps); ++i)
    delete[] u_xyz[i];
  delete[] u_xyz;

  for (i=0; i<=Nsteps; ++i) {
    delete[] Nint[i];
    delete[] Lint[i];
    delete[] u[i];
  }
  delete[] u;
  delete[] Nint;
  delete[] Lint;

  delete[] Cmn_list;

  return 0;
}
