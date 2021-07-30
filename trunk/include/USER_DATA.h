#include "mpi.h"
#include "mycomplex.h"

typedef struct {
  int iRows, iCols, memory;
  int **a;
} IntMatrix;

typedef struct {
  int iRows, iCols, memory;
  double **a;
} DoubleMatrix;

typedef struct {
  int iRows, iCols, memory;
  Complex **a;
} ComplexMatrix;

typedef struct {

  int comp1; //!< comp1, comp2, comp3 components of integer vector
  int comp2; 
  int comp3;

} VECTOR_INT;

typedef struct {

  double comp1; //!< comp1, comp2, comp3 components of double vector
  double comp2;
  double comp3;

} VECTOR_DOUBLE;

typedef struct {

  int memory;         //!< Memory allocated to this array
  VECTOR_INT oblique; //!< Components of k vectors in fractions of reciprocal lattice vectors expressed as integers
  VECTOR_DOUBLE cart; //!< Components of k vectors in Cartesians

} VECTOR_KNET; 

typedef struct {

  int monk_pack;        //!< not used
  int verbosity ;       //!< Printing level during runs. Default value 1
  int print_pairs ;     //!< Print pair_p, pair_c or pair_e when generated
  int check_pairs ;     //!< Check pair_p, pair_c or pair_e when generated
  int numtasks ;        //!< Number of tasks generated under MPI_WORLD_COMM
  int taskid ;          //!< Number of current task
  int type ;            //!< Job type  0 = scf, 1 = bse
  int *Fermi ;          //!< 
  int xc_num ;          //!< Number of functionals for xclib
  int xc_typ[2] ;       //!< XC_FAMILY parameter for xclib
  int xc_hfx ;          //!< Hamiltonian requires Hartree-Fock exchange
  int xc_lmx ;          //!< Max L value in Lededev grid for angular integration in DFT
  int xc_sph ;          //!< Number of angular integration points in DFT
  int xc_rad ;          //!< Number of radial integration points in DFT
  int xc_grd ;          //!< Grid type in DFT
  int ham_type ;        //!< Hamiltonian type: Set to 0 for DFT, 1 for hybrid DFT, 2 for Hartree-Fock
  int ham_dft  ;        //!< DFT Hamiltonian type: Set to 0 for DFT, 1 for DFT-GGA
  int ham_dft_num ;     //!< DFT Hamiltonian type: Number of density functionals needed - max value 2
  int ham_dft_exc[2] ;  //!< DFT Hamiltonian type: Exchange functional list
  int ham_dft_cor[2] ;  //!< DFT Hamiltonian type: Correlation functional list
  int ham_hyb  ;        //!< Hybrid DFT Hamiltonian type: Set to 1 for hybrid DFT
  int coul_int ;        //!< Calculate Coulomb integrals without multipole approximation
  int exch_int ;        //!< Calculate Exchange integrals without multipole approximation
  int spin_dim ;        //!< Used for summation over spin in spin-polarised(= 2)/unpolarised(= 1) cases
  int spin_fac ;        //!< Spin factor in spin-polarised ( =2)/unpolarised (= 1) cases
  int spin_tot ;        //!< Total spin for spin polarised cases: singlet = 0, doublet = 1, triplet = 2
  int spin_index ;      //!< Used to look over spin up and down
  int spin_pol ;
  int spin_orb ;        //!< Spin orbit coupling switch - on iff set to 1  && spin_pol == 1 - off by default
  int spin_polarisation;//!< Switch (= 0) for spin-unpolarised cases (= 1) for spin-polarised cases
  int mixing_type ;     //!< Parameter which determines type of mixing in SCF: 
  int max_cycle ;       //!< Maximum number of cycles in SCF calculation
  int iter ;            //!< Current iteration in SCF calculation
  int sym_adapt ;       //!< Symmetry adapt basis (= 1) Default is symmetry adaptation on
  int guess_type ;      //!< Type of guess used to start/restart SCF calculation
  //int mpi_io ;          //!< IO switch (= 0) use multiple local files (= 1) use single MPI files
  int mpp ;             //!< MPP switch (= 0) write eigenvectors to disk for density matrix, 1 node max (= 1) no eigenvector write use for >1 node 
  int scf_direct ;      //!< Set to 1 if SCF integrals are recalculated each cycle, set to 0 if they are stored on disk
  int scf_denfit ;      //!< Set to 1 if SCF integrals are calculated by density fitting
  int scf_dencou ;      //!< Set to 1 if density fitting coefficients are from 3-centre coulomb integrals, 0 if 3-centre overlap integrals
  int scf_coulomb ;     //!< Set to 1 if SCF coulomb integrals are calculated, set to 0 if they are not calculated
  int scf_exchange ;    //!< Set to 1 if SCF exchange integrals are calculated, set to 0 if they are not calculated
  int scf_trans ;       //!< Set to 0 canonical transformation set to 1 Cholesky transformation (default)
  int int_exist ;       //!< Set to 0 calculate 2e integrals (default) set to 1 2e integrals exist
  int int_exist_no_sym ;//!< Set to 0 calculate 2e integrals (default) set to 1 2e integrals exist case of no symmetry
  int gw ;              //!< Set to 1 if GW correction for BSE is to be calculated
                        //!< Set to 2 if full self-energy matrix is to be calculated and diagonalised
  int self_plot ;       //!< Set to 1 if energy dependence diagonal matrix elements of GW self energy is to be calculated
  int bse_ham ;         //!< BSE calculation type: 0 = TDHF, 1 = BSE (default), 2 = RPA
  int bse_tda ;         //!< BSE calculation type: 0 = full BSE matrix (default), 1 = Tamm-Dancoff approximation
  int bse_spin ;        //!< BSE calculation spin: 0 = singlet (default), 1 = triplet
  int bse_spk ;         //!< BSE calculation type: 0 = LAPACK (default), 1 = SCALAPACK linear algebra routines in BSE and TDHF
  int bse_int ;         //!< BSE calculation type: 0 = calculate integrals, 1 = integrals exist on disc, 2 = calculate in core
  int bse_lim ;         //!< BSE calculation: number of BSE eigenvectors to write to disc and use in optical spectrum calculation
  int bse_cou ;         //!< BSE calculation: calculate SCF HF exchange to test q == 0 integrals 
  int bse_exc ;         //!< BSE calculation: calculate SCF HF exchange to test q != 0 integrals 
  int band_dim ;        //!< Number of Bloch function product arrays in density fitting integrals
  int rpa_lim ;         //!< RPA calculation: number of RPA eigenvectors to write to disc and use in screened interaction calculation
  int gw_spk ;          //!< GW calculation type: 0 = LAPACK (default), 1 = SCALAPACK linear algebra in GW Casida calculation
  int gw_int ;          //!< GW calculation type: 0 = calculate integrals, 1 = integrals exist on disc, 2 = calculate in core
  int bse_denfit ;      //!< Set to 1 if BSE integrals are calculated by density fitting (must be 1 for screened exchange)
  int bse_screxc ;      //!< Set to 1 if BSE screened exchange integrals are calculated by density fitting
  int vectors ;         //!< eigenvectors to be read from disk if set to 1 to be read or written if set to 2
  int values ;          //!< eigenvalues to be read from disk if set to 1 to be read or written if set to 2
  int density ;         //!< density matrix to be read from disk if set to 1 to be read or written if set to 2
  int dimp;             //!< dimension of reduced density matrix (multiply by 2 for spin-polarised density matrix)
  int dimf;             //!< dimension of full density matrix (multiply by 2 for spin-polarised density matrix)
  int kpoints ;         //!< K point net type: Monkhorst-Pack = 0; list of points = 1
  int pms ;             //!< Switch (= 0/1) for permutation symmetry off/on
  int trs ;             //!< Switch (= 0/1) for time reversal symmetry off/on
  int sgs ;             //!< Switch (= 0/1) for space group symmetry off/on
  int rss ;             //!< Switch (= 0/1) for real space symmetry off/on
  int kss ;             //!< Switch (= 0/1) for reciprocal space symmetry off/on
  int l_max ;           //!< Maximum l value for wave function or auxillary bases
  int lmax ;            //!< Maximum l value for multipole moment expansion
  int lmax_fac ;        //!< Maximum number of monomials for given l value plus 1 for spheropole moment
  int field_dirs ;      //!< Number of field directions for susceptibilities
  int npoints ;         //!< Number of energy points at which a spectrum is to be sampled
  int nspectra ;        //!< Number of spectra to be calculated
  int C09 ;             //!< Data read from Crystal09 if set to 1
  int mixing_order ;    //!< Number of generations of density matrix used in Pulay mixing
  int diis ;            //!< Pulay DIIS mixing switch: turned on if set to 1
  int mxr ;             //!< Maximum side of grid for R->vec_ai real space lattice vectors 
  int mxg ;             //!< Maximum side of grid for G->vec_bi reciprocal space lattice vectors 
  int fix_occ ;         //!< Occupancy for density matrix calculation fixed if set to 1, determined by Fermi level if set to 0
  int fermi_k ;         //!< K point of Fermi level
  int lumo_k ;          //!< K point of LUMO level
  int *memory ;         //!< Counts current memory allocated, node by node
  int *max_memory ;     //!< Counts maximum memory allocated, node by node
  double fermi_energy ; //!< Fermi energy
  double lumo_energy ;  //!< LUMO energy
  double total_energy ; //!< Total SCF energy
  double energy_change; //!< Total SCF energy change in an iteration
  double twoe_energy;   //!< Total Coulomb + exchange + exchange-correlation energy in SCF calculation
  double nuc_nuc ;      //!< Nuclear-nuclear repulsion energy
  double fock_mixing ;  //!< Mixing fraction in density matrix mixing between 0 (no damping) and 1
  double scf_tol ;      //!< Parameter for SCF convergence on energy
  double itol1 ;        //!< Parameter for cutoff of Coulomb integrals in Ewald summation
  double itol2 ;        //!< Parameter for cutoff of Coulomb integrals in Ewald summation
  double itol3 ;        //!< Parameter for cutoff of exchange integrals in real space
  double itol4 ;        //!< Parameter for cutoff of exchange integrals in real space
  double itol5 ;        //!< Parameter for cutoff of exchange integrals in real space
  double overlap_tol_1; //!< Overlap criterion used to decide whether to retain pair_p pairs
  double overlap_tol_2; //!< Overlap criterion used to decide whether to retain pair_c pairs
  double overlap_tol_3; //!< Overlap criterion used to decide whether to retain pair_e pairs
  double overlap_tol_4; //!< Overlap criterion used to decide whether charges penetrate in point multipole approximation
  double overlap_tol_5; //!< Overlap criterion used to decide whether to retain pair_t lattice triples
  double overlap_tol_6; //!< Overlap criterion used to decide whether to retain pair_t lattice quads
  double overlap_tol_7; //!< Overlap criterion used to decide whether to retain pair_t point triples
  double overlap_tol_8; //!< Overlap criterion used to decide whether to retain pair_t point quads
  double overlap_tol_9; //!< Overlap criterion used to decide whether to retain pair_e quads
  double electron_count;//!< Total number of electrons per primitive unit cell
  double ham_hyb_wt;    //!< Fock exchange weight in hybrid DFT Hamiltonian
  double energy_range[2]; //!< Energy range associated with job->field_dirs
  double total_time ;   //!< Accumulates total wall time expired
  double linewidth ;    //!< linewidth for spectra plotting in eV
  double scissor;       //!< Scissor shift in periodic TDHF
  double scalefac;      //!< Scaling factor for electron-hole attraction in TDHF
  VECTOR_DOUBLE e_field[3];//!< Electric field direction cosines associated with job->field_dirs
  clock_t start ;       //!< not used
  clock_t end ;         //!< not used

} JOB_PARAM ; //!< Various job control parameters

typedef struct {

FILE *out ;            //!< Main output file
FILE *in ;             //!< Input file for atomic structure and basis sets. Must be named INPUT
FILE *err ;            //!< Errors directed to this file
FILE *job ;            //!< Input file for tasks. First argument in command line input
MPI_File fh ;          //!< File pointer for SCF eigenvectors for periodic systems
MPI_File gh ;          //!< File pointer for SCF eigenvalues for periodic systems
char directory1[100];  //!< Name of directory where files such as datafile are to be written or found
char scf_eigvec[100];  //!< not used
char cas_eigvec[100];  //!< Directory for Casida equation eigenvector file
char bse_eigvec[100];  //!< Directory for BSE eigenvector file

} FILES ;   //!< File and directory names

typedef struct {

  double lattice_a ; //!< Lattice parameters read in in Angstroms
  double lattice_b ;
  double lattice_c ;
  double alpha ;     //!< Angles between lattice vectors read in in degrees and converted to Rad.
  double beta ;
  double gamma ;
  double primitive_cell_volume ; //!< Primitive cell volume/area/length printed in Angstroms and then converted to atomic units
  VECTOR_DOUBLE conventional_cell[3] ; //!< Conventional unit cell vectors printed in Angstroms and the converted to atomic units
  VECTOR_DOUBLE primitive_cell[3] ; //!< Primitive unit cell vectors printed in Angstroms and the converted to atomic units
  VECTOR_DOUBLE reciprocal_cell[3] ; //!< Reciprocal lattice cell vectors printed in Angstroms and the converted to atomic units
  char type[8] ;  //!< Set to C, S, P for 3-D, 2-D or 1-D systems and M for finite systems
  char system ;   //!< Set to A, M, O, T, R, H, C for triclinic, monoclinic, orthorhombic, tetragonal, rhombohedral, hexagonal, cubic crystal classes
  char centrosymm ; //!< Set to C, N for centrosymmetric, non-centrosymmetric space groups
  char centring ; //!< Set to P, F, A, B, C for primitive, face-centred, A, B or C face centred Bravais lattices
  int origin ;    //!< Set to 0 or 1 depending on choice of origin for a particular space group
  int cell_choice ; //!< Set to integer (0 to 5) depending on unit cell choice (monoclinic crystal classes)

} CRYSTAL ; //!< Parameters relating to the crystal structure but excluding atomic coordinates

typedef struct {

  int memory;               //!< Memory allocated to this array
  int number_of_operators ; //!< Number of symmetry operators in space/rod/layer/point group
  int number_of_permutations;//!< Number of permutation operators in space/rod/layer/point group 8 by default, set to 1 to turn off permutation symmetry
  int number_of_classes;    //!< Number of conjugacy classes in rotation group
  int grp_dim;              //!< Dimension of symmetry group
  int *inverse ;            //!< Array containing inverses of each symmetry operator. 
  int *ind_i ;              //!< Indices of non-zero elements of symmetry operators in packed arrays
  int *ind_j ;    
  int *num_ij ;             //!< Number of non-zero symmetry operators in packed arrays
  int *op_shift ;           //!< Array containing integer offsets to rotation operators for s, sp, p, d and f states
  int *inr ;                //!< Symmetry operators stored as integers in basis of primitive unit cell vectors
  int *grp_pm ;             //!< Group table for permutations
  int *grp_k ;              //!< Group table for rotations
  int *cls_num_pm ;         //!< Group class sizes for permutations
  int *cls_num_k ;          //!< Group class sizes for rotations
  int *cls_pos_pm ;         //!< Group class positions for permutations
  int *cls_pos_k ;          //!< Group class positions for rotations
  int *cls_pm ;             //!< Group class membership for permutations
  int *cls_k ;              //!< Group class membership for rotations
  int *irp_dim_k ;          //!< Dimension of irreps for rotations
  double *irr ;             //!< Symmetry operators stored as doubles in Cartesian basis
  double *rot ;             //!< Values of non-zero elements of symmetry operators in packed arrays
  double *character_table ; //!< Character table for symmetry group
  VECTOR_DOUBLE *taur ;

} SYMMETRY ;                //!< Structure containing information on symmetry operators

typedef struct {

  double *mag;              //!< Magnitude of each lattice vector
  int memory;               //!< Memory allocated to this array
  int *num;                 //!< Number of lattice vectors in each shell
  int *vec_ai_ord;          //!< Order of lattice vectors after grouping in shells
  int *vec_ai_inv;          //!< Inverts order of lattice vectors back to original order
  int number_of_shells;     //!< Number of shells of lattice vectors up to last_vector
  int number_of_ewald_shells;  //!< Number of shells of lattice vectors up to last_ewald_vector
  int last_vector;          //!< Number of lattice vectors
  int last_ewald_vector;    //!< Number of lattice vectors for ewald sums
  int margin_vector;        //!< Last_vector plus a margin guaranteed to allow for addition or subtraction of pairs of O_temp vectors to/from last_vector
  int max_vector;           //!< Number of lattice vectors defined by MXRX, MXRY and MXRZ
  double cutoff;            //!< Distance in Angstroms at which pair generation is cut off
  VECTOR_DOUBLE *vec_ai;    //!< Cartesian components of each lattice vector

} REAL_LATTICE;             //!< Structure containing real space lattice vector components, magnitudes and total number

typedef struct {

  int *sumvec;
  int *diffvec;
  int *lattvec;
  int max_vector;
  int last_vector;
  int margin_vector;        //!< Last_vector plus a margin guaranteed to allow for addition or subtraction of pairs of O_temp vectors to/from last_vector
  int last_ewald_vector;    //!< Number of lattice vectors for density matrix
  int memory;               //!< Memory allocated to this array
  //int last_overlap_vector;  //!< Number of lattice vectors for density matrix
  //double cutoff;
  //double ewald_cutoff;
  //double margin;
  //double aspect_ratio;

} REAL_LATTICE_TABLES;      //!< Tables for addition, subtraction and rotation of real lattice vectors

typedef struct {

  VECTOR_DOUBLE *vec;
  double *EXPFAC;
  double *sqr;
  double *invsqr;
  double *x;
  double *y;
  double *z;
  double cutoff;
  double gamma_0_inv ;
  int max_vector;
  int last_vector;
  int memory;               //!< Memory allocated to this array

} Q_LATTICE;                //!< Structure containing reciprocal space Q vector components and magnitudes

typedef struct {

  VECTOR_DOUBLE *vec_bi;
  VECTOR_DOUBLE *vec_b2;
  double *mag;
  double *EXPFAC;
  DoubleMatrix *A;
  DoubleMatrix *B;
  double *sqr;
  double *invsqr;
  double *shell_mag;
  double *x;
  double *y;
  double *z;
  double cutoff;
  double gamma_0_inv ;
  int *num;
  int number_of_shells;
  int max_vector;
  int last_vector;
  int memory;               //!< Memory allocated to this array

} RECIPROCAL_LATTICE;       //!< Structure containing reciprocal space lattice vector components, magnitudes and shell structure

typedef struct {

  double *sc ;              //!< Coefficients of s, p, d, f and g Gaussian primitives in s, sp, p, d, f, g shell structure
  double *pc ;
  double *dc ;
  double *fc ;
  double *gc ;
  double *expo ;            //!< Exponents of Gaussian primitives in s, sp, p, d, f shell structure
  double *c_sh ;            //!< Coefficients of s, p, d, f and g Gaussian primitives in sh shell structure (sp shells split into s and p)
  double *expo_sh ;         //!< Exponents of Gaussian primitives in sh shell structure
  int memory;               //!< Memory allocated to this array

} GAUSSIAN ;

typedef struct {

  int *imax ;               //!< 
  int *ng ;                 //!< 
  int *type ;               //!<
  int *type1 ;              //!<
  int *ord ;                //!<
  int *opp ;
  int *imax_sh ; 
  int *ng_sh ; // ng
  int *type_sh ; //shtype obsolete
  int *type1_sh ; //shtype obsolete
  int *cart ; //shtype
  int *shar ; //shtype
  int *ord_sh ;
  int *opp_sh ;
  int *ind_i;
  int *ind_j;
  int *num_ij;
  int tuv[5][15][3] ; 
  //CHANGES2015 int tuv[8][6][3] ; 
  int memory;               //!< Memory allocated to this array
  double *nele; // initial number of electrons
  double *pop_sh ;
  double *pop ;
  double *rot;
  double *min_expo_sh; // smallest gaussian exponent on shell
  double *min_coef_sh; // smallest gaussian exponent on shell

} SHELL ;

typedef struct {

  int num_salc ;                             //!< Number of symmetry adapted linear combinations in an atom star
  int num_atom ;                             //!< Number of atoms in an atom star
  int total_coef ;                           //!< Number of coefficients in all symmetry adapted linear combinations
  //int total_atom_coef ;                      //!< Number of atoms in star * number of coefficients in all symmetry adapted linear combinations
  int memory ;                               //!< Number of bytes allocated to SALC array
  int *atm ;                                 //!< Pointer to atom position for each salc coefficient
  //int *ind_i ;                               //!< Pointer to atom basis function position for each salc coefficient
  int *irp ;                                 //!< Irreducible representation of each symmetry adapted linear combination
  int *num_irp ;                             //!< Number of irreducible representation in each class for a given atom
  int *num_coef ;                            //!< Number of non-zero coefficients in a symmetry adapted linear combination
  //double *coef ;                             //!< Coefficients of symmetry adapted linear combinations of basis functions in a shell
  IntMatrix *bfn_posn ;                      //!< Pointer to atom basis function position for each salc/atom
  DoubleMatrix *coeff ;                      //!< Coefficients of symmetry adapted linear combinations of basis functions in a shell

} SALC ;

#define MXSH 40
#define MXGF 20
//#define MXPR 6000
//#define MXP  10000
#define MXP  120000
//#define MXP  70000
//#define MXP  130000
//#define MXT  10000
#define MXT  25000

typedef struct {

  int shell_type;
  int number_of_Gaussians;
  double number_of_electrons;
  double exponent[MXGF];
  double coefc[MXGF];
  double coefp[MXGF];

} SHELLTYPE;                                 //!< Temporary array used in reading basis set

typedef struct {

  int atomic_number;
  int number_of_shells;
  int number_of_sh_shells;
  VECTOR_DOUBLE posn;
  SHELLTYPE shell[MXSH];

} ATOMTYPE;                                  //!< Temporary array used in reading basis set

typedef struct {

  int memory;               //!< Memory allocated to this array
  int lmax;
  int *numb;
  int *l;
  int *n;
  DoubleMatrix *coef;
  DoubleMatrix *expo;

} PSEUDO;

typedef struct {

  int *nshel_sh ;                            //!< Number of shells on each atom with spherical harmonic d, f states
  int *shelposn_sh ;                         //!< First shell on each atom with spherical harmonic d, f states
  int *gausposn_sh ;                         //!< First gaussian primitive on each atom with spherical harmonic d, f states
  int *shlposn_sh ;                          //!< First shell on each atom with spherical harmonic d, f states
  int *bfnposn_sh ;                          //!< First basis function on each atom with spherical harmonic d, f states
  int *bfnnumb_sh ;                          //!< Number of basis functions on each atom with spherical harmonic d, f states
  int *nshel ;                               //!< Number of shells on each atom with Cartesian d, f states
  int *shelposn ;                            //!< First shell on each atom with Cartesian d, f state
  int *gausposn ;                            //!< First gaussian primitive on each atom with Cartesian d, f states
  int *shlposn ;                             //!< First shell on each atom with Cartesian d, f states
  int *bfnposn ;                             //!< First basis function on each atom with Cartesian d, f states
  int *bfnnumb ;                             //!< Number of basis functions on each atom with Cartesian d, f states
  int *atomic_number;                        //!< Atomic numbers of all atoms
  int *magnetic;                             //!< Array containing identities of magnetic ions +1/+2 non-magnetic/magnetic
  int *spin;                                 //!< Set to +1/-1 for spin up/spin down initial guess on particular atom
  int *basis_set;                            //!< Basis set attached to all atoms
  int *uniq;                                 //!< Unique atom from which this atom is derived by a symmetry operation
  int memory;                                //!< Memory allocated to this array
  int number_of_unique_atoms;                //!< Number of unique atoms in unit cell
  int number_of_basis_sets;                  //!< Number of basis sets
  int number_of_atoms_in_unit_cell;          //!< Total number of atoms in unit cell
  int number_of_shells_in_unit_cell;         //!< Number of shells in basis with Cartesian d, f states
  int number_of_sh_shells_in_unit_cell;      //!< Number of shells in basis with spherical harmonic d, f states
  int number_of_bfns_in_unit_cell;           //!< Number of basis functions with Cartesian d, f states
  int number_of_sh_bfns_in_unit_cell;        //!< Number of basis functions with spherical harmonic d, f states
  int number_of_gaussians_in_unit_cell;      //!< Number of gaussians in basis with Cartesian d, f states
  int number_of_sh_gaussians_in_unit_cell;   //!< Number of gaussians in basis with spherical harmonic d, f states
  double number_of_electrons_in_unit_cell;   //!< Number of electrons
  double *pop;                               //!< Electron population on each atom
  VECTOR_DOUBLE *cell_vector;                //!< All atomic positions in Cartesian coordinates

} ATOM;                                      //!< Structure containing basis set information on each atom

typedef struct {

  int *K;                                    //!< Symmetry operator which generates equivalent atom from unique atom
  int *O;                                    //!< Lattice vector which maps atom back to zeroth cell after symmetry operation
  int *P;                                    //!< Set to 1 if an atom transforms to itself modulo a lattice translation, zero otherwise
  int *numb;                                 //!< Number of atoms equivalent to a particular unique atom
  int *posn;                                 //!< Index for beginning of nth set of atoms which are equivalent by symmetry
  int memory;                                //!< Memory allocated to this array
  double *expo;                              //!< Smallest gaussian primitive exponent on each atom
  double *coef;                              //!< Expansion coefficient corresponding to *expo

} ATOM_TRAN;                                 //!< Information on unique atoms and how they transform under symmetry operations

typedef struct {

  int cell1[MXP];
  int cell2[MXP];
  int latt1[MXP];
  int latt2[MXP];
/*
  int ptr[MXP];
  int P[MXP];
  int numb;
  int orc;
  int ord;
  double rc1;
  double rc2;
  double rsq;
  double rs1;
  double rs2;
  double sab;
  double scd;
  double sac;
  double sbd;
  double ovp;
*/

} PAIR_TMP;                                  //!< Temporary array similar to PAIR_TRAN used to count pairs before allocation

typedef struct {

  int *cell1;              //!< Index for first atom in a pair
  int *cell2;              //!< Index for second atom in a pair
  int *latt1;              //!< Index for first lattice vector in a pair (zero in many cases)
  int *latt2;              //!< Index for second lattice vector in a pair
  int *lattd;              //!< Index for lattice vector corresponding to *latt2 - *latt1 for coulomb integrals
  int *k;                  //!< Symmetry operator which generates a pair from a unique pair
  int *p;                  //!< Permutation operator which generates a pair from a unique pair
  int *O;                  //!< Lattice vector which translates *cell1 for a pair to the zeroth cell after a symmetry operation
  int *P;                  //!< Set to 1 if a pair transforms to itself modulo a lattice vector (is fixed) under a symmetry operation
  int *posn;               //!< Index for beginning of nth set of pairs which are equivalent by symmetry
  //int *posi;               //!< Index for beginning of nth set of fixed pairs
  int *off;                //!< Offset in full density matrix at which nth pair begins
  int *Off;                //!< Offset in reduced density matrix at which n'th pair begins
  int *ptr;                //!< Index of a particular pair in a list of pairs i.e. (i,j,g) stored at *ptr
  int *Ptr;
  int *uniq;
  int *rot;                //!< Transformation of unique pairs under each permutation and space group operation
  int *numb;               //!< Number of pairs equivalent by symmetry to the nth unique pair
  int *K;                  //!< Symmetry operator which generates nth fixed pair
  //int *Q;                  //!< Permutation operator which generates nth fixed pair
  //int *I;                  //!< Lattice vector which translates *cell1 for nth **fixed** pair to the zeroth cell after a symmetry operation
  //int *numi;               //!< Number of pairs equivalent by symmetry to nth fixed pair
  //int *limit;              //!< Maximum lattice vector for each pair of unique atoms
  int nump;                //!< Number of unique pairs in this structure
  //int numq;                //!< Number of quads in this structure
  //int numt;                //!< Number of triples in this structure
  int tot;                 //!< Total number of pairs or triples in this structure
  int memory;              //!< Memory allocated to this array
  //int maxl;                //!< Maximum real space lattice vector in this structure
  double *dist;            //!< Distance between atoms in nth pair
  double cutoff;           //!< Cutoff distance between atoms for range selected pairs
  //double *ovp;             //!< Overlap between gaussians with smallest exponents in nth pair

} PAIR_TRAN;               //!< Structure containing atom/lattice indices for ion pairs and their transformations under symmetry operations

typedef struct {

  int *cell1;              //!< Index for first atom in a pair
  int *cell2;              //!< Index for second atom in a pair
  int *cell3;              //!< Index for first atom in a pair
  int *cell4;              //!< Index for second atom in a pair
  int *latt1;              //!< Index for first lattice vector in a pair (zero in many cases)
  int *latt2;              //!< Index for second lattice vector in a pair
  int *latt3;              //!< Index for first lattice vector in a pair
  int *latt4;              //!< Index for second lattice vector in a pair
  int *k;                  //!< Symmetry operator which generates a quad from a unique quad
  int *p;                  //!< Permutation operator which generates a quad from a unique quad
  int tot;                 //!< Total number of quads in this structure
  int memory;              //!< Memory allocated to this array

} QUAD_TRAN;               //!< Structure containing atom/lattice indices for ion quads and their transformations under symmetry operations

typedef struct {

  int *cell1;              //!< Index for first atom 
  int *cell2;              //!< Index for first atom in a pair
  int *cell3;              //!< Index for second atom in a pair
  int *latt1;              //!< Index for first lattice vector (zero in many cases)
  int *latt2;              //!< Index for first lattice vector in a pair
  int *latt3;              //!< Index for second lattice vector in a pair
  int *posn;               //!< Index for beginning of nth set of triples which are equivalent by symmetry
  int *k;                  //!< Symmetry operator which generates a triple from a unique triple
  int *p;                  //!< Permutation operator which generates a triple from a unique triple
  int *numb;               //!< Number of triples equivalent by symmetry to the nth unique triple
  int nump;                //!< Number of unique triples in this structure
  int tot;                 //!< Total number of pairs or triples in this structure
  int memory;              //!< Memory allocated to this array

} TRIPLE_TRAN;             //!< Structure containing atom/lattice indices for ion triples and their transformations under symmetry operations

typedef struct {

  int cell1[MXT];          //!< Index for first atom 
  int cell2[MXT];          //!< Index for first atom in a pair
  int cell3[MXT];          //!< Index for second atom in a pair
  int latt1[MXT];          //!< Index for first lattice vector (zero in many cases)
  int latt2[MXT];          //!< Index for first lattice vector in a pair
  int latt3[MXT];          //!< Index for second lattice vector in a pair

} TRIPLE_TMP;              //!< Structure containing atom/lattice indices for ion triples and their transformations under symmetry operations

/*
typedef struct {

  int ibz;
  int fbz;
  int opr;
  int num;
  int elem;
  int uniq;
  int fbzo;
  double wght;
  double kwgt;

} KPT_TRAN ;
*/

typedef struct {

  int unique;              //!< Number of nkunique k points
  int nktot;               //!< Total number of k points
  int *ibz;                //!< List of nkunique points in irreducible Brillouin zone
  int *fbz;                //!< Pointer from one of ksize points in full Brillouin zone to unique point
  int *bz;                 //!< List of ksize points in full Brillouin zone in order
  int *opr;                //!< Operator taking unique k point to point in full Brillouin zone
  int *num;                //!< Number of kpoints equivalent to each unique k point
  int *trs;                //!< trs = 1(0) point (not) generated by time reversal symmetry
  int memory;              //!< Memory allocated to this array
  double *weight;          //!< Weight of each unique k point for Brillouin zone integration
  VECTOR_INT *oblique;     //!< Coordinates of ksize k points in integer, oblique coordinates
  VECTOR_DOUBLE *cart;     //!< Coordinates of ksize k points in Cartesian coordinates

} KPOINT_TRAN ;            //<! Structure containing k point coordinates and transformations

typedef struct {

  int unique;              //!< Number of nkunique k,k+q point pairs
  int *ibz;                //!< q point bz value in corresponding knet Brillouin zone
  int *bz1;                //!< k point bz value in corresponding knet Brillouin zone
  int *bz2;                //!< k+q point bz value in corresponding knet Brillouin zone
  int *opr;                //!< Operator taking unique k,k+q point pair to equivalent pair
  int *num;                //!< Number of kpoints equivalent to each unique k,k+q point pair
  int *trs;                //!< trs = 1(0) point (not) generated by time reversal symmetry
  int memory;              //!< Memory allocated to this array

} KQPOINT_TRAN ;            //<! Structure containing k,k+q pair point coordinates and transformations

typedef struct {

  int unique;              //!< Number of nkunique k points
  int nktot;               //!< Total number of k points
  int *ibz1;               //!< List of nkunique points in irreducible Brillouin zone
  int *ibz2;               //!< List of nkunique points in irreducible Brillouin zone
  int *fbz1;               //!< Pointer from one of ksize points in full Brillouin zone to unique point
  int *fbz2;               //!< Pointer from one of ksize points in full Brillouin zone to unique point
  int *bz1;                //!< List of ksize points in full Brillouin zone in order
  int *bz2;                //!< List of ksize points in full Brillouin zone in order
  int *opr;                //!< Operator taking unique k point to point in full Brillouin zone
  int *num;                //!< Number of kpoints equivalent to each unique k point
  int *posn;               //!< Positions in bz lists at which sets of equivalent k point pairs begin
  //int *trs;                //!< trs = 1(0) point (not) generated by time reversal symmetry
  int memory;              //!< Memory allocated to this array
  //double *weight;          //!< Weight of each unique k point for Brillouin zone integration
  //VECTOR_INT *oblique1;    //!< Coordinates of ksize k points in integer, oblique coordinates
  //VECTOR_INT *oblique2;    //!< Coordinates of ksize k points in integer, oblique coordinates
  //VECTOR_DOUBLE *cart1;    //!< Coordinates of ksize k points in Cartesian coordinates
  //VECTOR_DOUBLE *cart2;    //!< Coordinates of ksize k points in Cartesian coordinates

} KPOINT_PAIR_TRAN ;       //<! Structure containing k point pair coordinates and transformations

/*
typedef struct {

  int pt1[48];
  int pt2[48];
  int numb;
  double wght;

} KPAIR_TMP;
*/

/*
typedef struct {

  int *pt1;
  int *pt2;
  int *k;
  int numb;
  int uniq;
  double wght;

} KPAIR_TRAN;
*/

typedef struct {

  int *i;
  int *j;
  int *k;
  int *l;
  int *gj;
  int *gk;
  int *gl;
  int *offset1;
  int *offset2;
  int num;
  int quad;
  int memory;              //!< Memory allocated to this array
  double *value;

} INTEGRAL_LIST;             //!< Structure containing basis function, lattice indices and values of 4-centre coulomb integrals

typedef struct {

  int *i;
  int *j;
  int *k;
  int *l;
  int *gj;
  int *gk;
  int *gl;
  int *offset1;
  int *offset2;
  int num;
  int quad;
  int memory;              //!< Memory allocated to this array
  Complex *value;

} INTEGRAL_LIST_COMPLEX;   //!< Structure containing basis function, lattice indices and values of 4-centre coulomb integrals

typedef struct {

  int is[3];
  int bands[4];
  int plot_bands[10];          //!< List of bands for plotting
  int homo[2];
  //int homo_up;
  //int homo_dn;
  int nkunique;
  int nktot;
  int ksize;
  int npoints;                      //!< Number of k points in list fermi->knet_list
  int target_states;                //!< Total number of electrons * number of k points
  int memory;                       //!< Memory allocated to this array
  int *occupied;
  double *occupation;
  KPOINT_TRAN *knet;
  //VECTOR_KNET *knet_unique;
  VECTOR_KNET *knet_list;

} FERMI;

typedef struct {

  int memory;                       //!< Memory allocated to this array
  double *Fock;
  double *Kinetic;
  double *ElecNuc;
  double *Coulomb;
  double *Momentum;
  double *Overlap;
  double *Grad_Grad;
  double *Dipole;

} INT_1E;                                    //!< Structure containing 1-body operators in real space

typedef struct {

  int npoints;
  int memory;                       //!< Memory allocated to this array
  int *posn;
  int *numb;
  double *r;
  double *y;
  double *weight;
  VECTOR_DOUBLE *x;

} DFT_GRID;                                    //!< Structure containing grid for DFT integration
