#include <cstring>
#include "mycomplex.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "ROTATIONS_MOLECULE.h"
#include "PRINT_MOLECULE.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "SCF_ATOM.h"
#include "DENSITY_MATRIX_MOLECULE.h"

using namespace std;

void initial_density_matrix(double *P0, double *F0, PAIR_TRAN *pair_p, FERMI *fermi, ATOM *atoms,  ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file){

  int i, j, k, l, m, p, q, r, s, u;
  int dim, dim1, dimk, dimp, dims;
  int bfposi;
  int nd1, nd2;
  int vector_size, offset, nkunique;
  int count;
  int num_sym = symmetry->number_of_operators;
  int taken[atoms->number_of_basis_sets];
  int occupied[job->spin_dim];
  double fac;
  double electron_count_up, electron_count_down, electron_occupation_up, electron_occupation_down;
  double total_electron_occupation_up, total_electron_occupation_down;
  double time1, time2;
  double *eigenvalues;
  double *occupation;
  DoubleMatrix *eigenvectors, *P;

  num_sym = 1;
  nkunique = fermi->nkunique;
  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  dimk = job->spin_dim * nkunique * atoms->number_of_sh_bfns_in_unit_cell;
  vector_size = dimk * dim1;
  time1 = MPI_Wtime();
  for (i = 0; i < atoms->number_of_basis_sets; i++) taken[i] = -1;
  total_electron_occupation_up = k_zero;
  total_electron_occupation_down = k_zero;
  job->electron_count = k_zero;

  // ******************************************************************************************
  // * Loop over all basis sets and unique atoms                                              *
  // ******************************************************************************************

  for (i = 0; i < atoms->number_of_basis_sets; i++) {
    for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
      if (atoms->basis_set[j] == i) {
        if (taken[i] == -1) {
          electron_occupation_up = k_zero;
          electron_occupation_down = k_zero;
          electron_count_up = k_zero;
          electron_count_down = k_zero;
          dim  = atoms->bfnnumb_sh[j];
          dimp = atoms->magnetic[j] * atoms->bfnnumb_sh[j];
          dims = job->spin_dim * atoms->bfnnumb_sh[j];
          AllocateDoubleMatrix(&P,&dimp,&dim,job);
          AllocateDoubleMatrix(&eigenvectors,&dimp,&dim,job);
          AllocateDoubleArray(&eigenvalues,&dims,job);
          AllocateDoubleArray(&occupation,&dims,job);
          ResetDoubleMatrix(P);
          ResetDoubleMatrix(eigenvectors);
          ResetDoubleArray(eigenvalues,&dims);
          ResetDoubleArray(occupation,&dims);
          atom_scf(atoms,j,eigenvectors,eigenvalues,occupation,occupied,shells,gaussians,job,file);
          if (atoms->magnetic[j] == 1)
            for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
              electron_occupation_up   += occupation[m];
             }
          if (atoms->magnetic[j] == 2 && atoms->spin[j] > 0)
          //if (atoms->magnetic[j] == 2 && atoms->spin[j] == 1)
            for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
              electron_occupation_up   += occupation[m];
              electron_occupation_down += occupation[atoms->bfnnumb_sh[j] + m];
             }
          if (atoms->magnetic[j] == 2 && atoms->spin[j] < 0)
          //if (atoms->magnetic[j] == 2 && atoms->spin[j] == -1)
            for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
              electron_occupation_up   += occupation[atoms->bfnnumb_sh[j] + m];
              electron_occupation_down += occupation[m];
             }
          for (s = 0; s < atoms->magnetic[j]; s++) {
            for (m = 0; m < occupied[s]; m++) {
              if (s == 0) electron_count_up   += k_one;
              if (s == 1) electron_count_down += k_one;
             }
            }
          //for (m = 0; m < dims; m++) fprintf(file.out,"%3d %10.4f\n",m,occupation[m]);
        fprintf(file.out,"%d %d Occup/dn %3d %3d occup/dn %10.4f %10.4f totup/dn  ctup/dn %10.4f %10.4f job->spin_dim %3d\n",\
        atoms->magnetic[j],atoms->spin[j],occupied[0],occupied[1],electron_occupation_up,electron_occupation_down,electron_count_up,\
        electron_count_down,job->spin_dim);

  // ******************************************************************************************
  // * Calculate density matrix for unique atom j                                             *
  // ******************************************************************************************

          for (s = 0; s < atoms->magnetic[j]; s++) {
            offset = s * atoms->bfnnumb_sh[j];
            fac = k_one;
            if (job->spin_dim == 1) fac = electron_occupation_up / electron_count_up;
            //if (s == 0 && occupied[s] > 0) fac = electron_occupation_up   / electron_count_up;
            //if (s == 1 && occupied[s] > 0) fac = electron_occupation_down / electron_count_down;
            //fprintf(file.out,"fac %f %f %f %f\n",electron_occupation_up,electron_count_up,electron_occupation_down,electron_count_down);
            for (k = 0; k < dim; k++) {
              for (l = 0; l < dim; l++) {
                for (m = 0; m < occupied[s]; m++) {
                  P->a[offset + k][l] += fac * eigenvectors->a[m][k] * eigenvectors->a[m][l];
                 }
                }
               }
              }
          if (job->taskid == 0 && job->verbosity > 1) {
            fprintf(file.out,"eigvec2 %3d %3d\n",dims,dim);
            ////if (job->spin_polarisation == 0)
            ////print_real_eigenvector_matrix(eigenvectors, eigenvalues, file);
            ////if (job->spin_polarisation == 1)
            ////PrintDoubleMatrix(eigenvectors, file);
            fprintf(file.out,"Atom density matrix\n");
            for (k = 0; k < P->iRows; k++) {
            for (l = 0; l < P->iCols; l++) {
            fprintf(file.out,"%6.3f ",P->a[k][l]); }
            fprintf(file.out,"\n"); }
            fflush(file.out);
           }
          taken[i] = 0;
         } // close if (taken[i]

  // ******************************************************************************************
  // * Loop over atoms symmetry equivalent to unique atom j                                   *
  // * Assemble reduced and full density matrices P0 and F0                                   *
  // ******************************************************************************************

  //for (s = 0; s < job->spin_dim; s++) {
  //for (s = 0; s < atoms->magnetic[j]; s++) {
    //if (atoms->spin[j] == -1) u = 1 - s;
    //if (atoms->spin[j] ==  0) u = 0;
    //if (atoms->spin[j] ==  1) u = s;
    //offset = s * atoms->bfnnumb_sh[j];
    //fprintf(file.out,"j %d s %d u %d\n",j,s,u);

    bfposi = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[q]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[q]];
      if (pair_p->cell1[q] == j && pair_p->cell2[q] == j && pair_p->latt1[q] == 0 && pair_p->latt2[q] == 0) {
        if (job->spin_dim == 1) {
        count = 0;
        for (k = 0; k < dim; k++) {
          for (l = 0; l < dim; l++) {
            //P0[s * job->dimp + bfposi + count] = P->a[offset + k][l];
            P0[bfposi + count] = P->a[k][l];
            F0[bfposi + count] = P->a[k][l];
            count++;
           }
          }
         }
        if (job->spin_dim == 2 && atoms->magnetic[j] == 1) {
        count = 0;
        for (k = 0; k < dim; k++) {
          for (l = 0; l < dim; l++) {
            P0[bfposi + count] = P->a[k][l];
            P0[bfposi + count + job->dimp] = P->a[k][l];
            F0[bfposi + count] = P->a[k][l];
            F0[bfposi + count + job->dimf] = P->a[k][l];
            count++;
           }
          }
         }
        if (job->spin_dim == 2 && atoms->magnetic[j] == 2) {
        count = 0;
        for (k = 0; k < dim; k++) {
          for (l = 0; l < dim; l++) {
            if (atoms->spin[j] > 0) {
            //if (atoms->spin[j] == 1) {
            P0[bfposi + count] = P->a[k][l];
            P0[bfposi + count + job->dimp] = P->a[atoms->bfnnumb_sh[j] + k][l];
            F0[bfposi + count] = P->a[k][l];
            F0[bfposi + count + job->dimf] = P->a[atoms->bfnnumb_sh[j] + k][l];
           }
            else if (atoms->spin[j] < 0) {
            //else if (atoms->spin[j] == -1) {
            P0[bfposi + count] = P->a[atoms->bfnnumb_sh[j] + k][l];
            P0[bfposi + count + job->dimp] = P->a[k][l];
            F0[bfposi + count] = P->a[atoms->bfnnumb_sh[j] + k][l];
            F0[bfposi + count + job->dimf] = P->a[k][l];
           }
            count++;
           }
          }
         }
        }
       bfposi += nd1 * nd2;
      }

  // ******************************************************************************************
  // * Count total number of spin up and spin down electrons                                  *
  // ******************************************************************************************

  if (atoms->magnetic[j] == 1)
  for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
    total_electron_occupation_up   += occupation[m] / two;
    total_electron_occupation_down += occupation[m] / two;
    job->electron_count += occupation[m];
   }
  if (atoms->magnetic[j] == 2 && atoms->spin[j] > 0)
  //if (atoms->magnetic[j] == 2 && atoms->spin[j] == 1)
  for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
    total_electron_occupation_up   += occupation[m];
    total_electron_occupation_down += occupation[atoms->bfnnumb_sh[j] + m];
    job->electron_count += occupation[m] + occupation[atoms->bfnnumb_sh[j] + m];
   }
  else if (atoms->magnetic[j] == 2 && atoms->spin[j] < 0)
  //else if (atoms->magnetic[j] == 2 && atoms->spin[j] == -1)
  for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
    total_electron_occupation_up   += occupation[atoms->bfnnumb_sh[j] + m];
    total_electron_occupation_down += occupation[m];
    job->electron_count += occupation[m] + occupation[atoms->bfnnumb_sh[j] + m];
   }
    //fprintf(file.out,"electron count %e occupationup/dn %10.4f %10.4f\n",\
    job->electron_count,total_electron_occupation_up,total_electron_occupation_down);
   } // close if (atoms
  } // close loop on j
   if (taken[i] == 0) {
     DestroyDoubleMatrix(&eigenvectors,job);
     DestroyDoubleArray(&eigenvalues,&dims,job);
     DestroyDoubleArray(&occupation,&dims,job);
     DestroyDoubleMatrix(&P,job);
     taken[i]++;
    }
   } // close loop on i

  // ******************************************************************************************
  // * Generate fermi->occupied and fermi->occupation arrays                                  *
  // ******************************************************************************************

  for (s = 0; s < job->spin_dim; s++) {
    for (k = 0; k < nkunique; k++) {
      fermi->occupied[s * nkunique + k] = 0;
      for (i = 0; i < fermi->occupied[s]; i++) {
      fermi->occupation[(s * nkunique  + k) * dim1 + i] = k_zero;
     }
    }
   }

  switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':

    for (s = 0; s < job->spin_dim; s++) {
      if      (s == 0) fermi->occupied[s] = (int)(total_electron_occupation_up   + 0.00001);
      else if (s == 1) fermi->occupied[s] = (int)(total_electron_occupation_down + 0.00001);
      for (k = 0; k < nkunique; k++) {
        for (i = 0; i < fermi->occupied[s]; i++) {
          fermi->occupation[(s * nkunique  + k) * dim1 + i] = k_one;
          //fprintf(file.out,"spin %3d occupied %10.4f\n",s,fermi->occupation[(s * nkunique  + k) * dim1 + i]);
         }
        }
       }

    for (s = 0; s < job->spin_dim; s++) {
      for (k = 0; k < nkunique; k++) {
        //if      (s == 0) fermi->occupied[k]            = 13;
        //else if (s == 1) fermi->occupied[nkunique + k] = 11;
        if      (s == 0) fermi->occupied[k]            = (int)(total_electron_occupation_up   + 0.00001);
        else if (s == 1) fermi->occupied[nkunique + k] = (int)(total_electron_occupation_down + 0.00001);
        //fprintf(file.out,"spin %3d k %3d occupied %3d\n",s,k,fermi->occupied[s * nkunique + k]);
         }
        }

   break;

    case 'M':

    if (job->spin_dim == 1) fermi->occupied[0] = (int)(total_electron_occupation_up + 0.00001);
    else if (job->spin_dim == 2) {
      fermi->occupied[0] = (int)(total_electron_occupation_up   + 0.00001);
      fermi->occupied[1] = (int)(total_electron_occupation_down + 0.00001);
      for (s = 0; s < job->spin_dim; s++) {
        for (i = 0; i < fermi->occupied[s]; i++) {
          fermi->occupation[s * dim1 + i] = k_one;
          //fprintf(file.out,"spin %3d occupied %3d occupation %lf\n",s,fermi->occupied[s],fermi->occupation[s * dim1 + i]);
         }
        }
       }

     break;
  
  } // close switch

  for (i = 0; i < atoms->number_of_basis_sets; i++) {
      if (taken[i] == -1) {
        if (job->taskid == 0 || job->verbosity > 1)
        fprintf(file.out,"Initial guess for atom %d failed in initial_density_matrix_crystal\n",j);
        printf("Initial guess for atom %d failed in initial_density_matrix_crystal\n",j);
        MPI_Finalize();
        exit(0);
      }
     }

  time2 = MPI_Wtime();

}

void read_density_matrix(FERMI *fermi, double **P0, int *dimp, int *dimf, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i;
  char buf2[110];
  char yy[20] = "/new_density_matrix";;
  MPI_File fh;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,yy);

  if (*dimp == 0) {
    int job_parameters[6];
    MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
    MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    MPI_File_read(fh, job_parameters, 6, MPI_INT, MPI_STATUS_IGNORE);
    fermi->is[0] = job_parameters[0];
    fermi->is[1] = job_parameters[1];
    fermi->is[2] = job_parameters[2];
    fermi->nkunique = job_parameters[3];
    job->dimp = job_parameters[4];
    job->dimf = job_parameters[5];
   *dimp = job_parameters[4];
   *dimf = job_parameters[5];
    MPI_File_close(&fh);
    return;
  }

  if (*dimp > 0) {
    int dim = job->spin_dim * fermi->nkunique * atoms->number_of_sh_bfns_in_unit_cell;
    int job_parameters[41 + job->spin_dim * fermi->nkunique];
    double job_parameters1[11 + dim];
    MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
    MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    MPI_File_read(fh, job_parameters,  41 + job->spin_dim * fermi->nkunique, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read(fh, job_parameters1, 11 + dim, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_read(fh, (*P0), job->spin_dim * job_parameters[4], MPI_DOUBLE, MPI_STATUS_IGNORE);
    read_JOB_PARAM_array(fermi, job_parameters, job_parameters1, job, file);
    for (i = 0; i < job->spin_dim * fermi->nkunique; i++) {
    fermi->occupied[i] = job_parameters[41 + i];
   }
    job->electron_count = k_zero;
    for (i = 0; i < dim; i++) {
    fermi->occupation[i] = job_parameters1[11 + i];
    job->electron_count += fermi->occupation[i] / fermi->nkunique; // will need weights here when k symmetry included
   }
    MPI_File_close(&fh);
  }

}

void read_JOB_PARAM_array(FERMI *fermi, int *job_parameters, double *job_parameters1, JOB_PARAM *job, FILES file)

{

int i;

  // ******************************************************************************************
  // * State of job is captured in job_parameters and job_parameters1 arrays                  *
  // * Parameters which are needed to recreate a particular state are written to disk         *
  // ******************************************************************************************

  //int monk_pack;        //!< not used
  //int verbosity ;       //!< Printing level during runs. Default value 1
  //int print_pairs ;     //!< Print pair_p, pair_c or pair_e when generated
  //int check_pairs ;     //!< Check pair_p, pair_c or pair_e when generated
  //int numtasks ;        //!< Number of tasks generated under MPI_WORLD_COMM
  //int taskid ;          //!< Number of current task
  //int *Fermi ;          //!< 
  //int max_cycle ;       //!< Maximum number of cycles in SCF calculation
  //int iter ;            //!< Current iteration in SCF calculation
  //int guess_type ;      //!< Type of guess used to start/restart SCF calculation
  //int mpi_io ;          //!< IO switch (= 0) use multiple local files (= 1) use single MPI files
  //int scf_direct ;      //!< Set to 1 if SCF integrals are recalculated each cycle, set to 0 if they are stored on disk
  //int vectors ;         //!< eigenvectors to be read from disk if set to 1 to be read or written if set to 2
  //int values ;          //!< eigenvalues to be read from disk if set to 1 to be read or written if set to 2
  //int density ;         //!< density matrix to be read from disk if set to 1 to be read or written if set to 2
  //int kpoints ;         //!< K point net type: Monkhorst-Pack = 0; list of points = 1
  //int field_dirs ;      //!< Number of field directions for susceptibilities
  //int npoints;          //!< Number of energy points at which a spectrum is to be sampled
  //int mixing_order ;    //!< Number of generations of density matrix used in Pulay mixing
  //int job->*memory ;    //!< Counts current memory allocated, node by node
  //int *max_memory ;     //!< Counts maximum memory allocated, node by node
  //int fix_occ ;         //!< Occupancy for density matrix calculation fixed if set to 1, determined by Fermi level if set to 0

  fermi->is[0]           = job_parameters[0];
  fermi->is[1]           = job_parameters[1];
  fermi->is[2]           = job_parameters[2];
  fermi->nkunique        = job_parameters[3];
  job->dimp              = job_parameters[4];  //!< dimension of reduced density matrix (multiply by 2 for spin-polarised density matrix)
  job->dimf              = job_parameters[5];  //!< dimension of full density matrix (multiply by 2 for spin-polarised density matrix)
  job->type              = job_parameters[6];  //!< Job type  0 = scf, 1 = bse
  job->xc_num            = job_parameters[7];  //!< Number of functionals for xclib
  job->xc_typ[0]         = job_parameters[8];  //!< XC_FAMILY parameter for xclib
  job->xc_typ[1]         = job_parameters[9];  //!< XC_FAMILY parameter for xclib
  job->xc_hfx            = job_parameters[10]; //!< Hamiltonian requires Hartree-Fock exchange
  job->xc_lmx            = job_parameters[11]; //!< Max L value in Lededev grid for angular integration in DFT
  job->xc_sph            = job_parameters[12]; //!< Number of angular integration points in DFT
  job->xc_rad            = job_parameters[13]; //!< Number of radial integration points in DFT
  job->xc_grd            = job_parameters[14]; //!< Grid type in DFT
  job->ham_type          = job_parameters[15]; //!< Hamiltonian type: Set to 0 for DFT, 1 for hybrid DFT, 2 for Hartree-Fock
  job->ham_dft           = job_parameters[16]; //!< DFT Hamiltonian type: Set to 0 for DFT, 1 for DFT-GGA
  job->ham_dft_num       = job_parameters[17]; //!< DFT Hamiltonian type: Number of density functionals needed - max value 2
  job->ham_dft_exc[0]    = job_parameters[18]; //!< DFT Hamiltonian type: Exchange functional list
  job->ham_dft_exc[1]    = job_parameters[19]; //!< DFT Hamiltonian type: Exchange functional list
  job->ham_dft_cor[0]    = job_parameters[20]; //!< DFT Hamiltonian type: Correlation functional list
  job->ham_dft_cor[1]    = job_parameters[21]; //!< DFT Hamiltonian type: Correlation functional list
  job->ham_hyb           = job_parameters[22]; //!< Hybrid DFT Hamiltonian type: Set to 1 for hybrid DFT
  job->coul_int          = job_parameters[23]; //!< Calculate Coulomb integrals without multipole approximation
  job->exch_int          = job_parameters[24]; //!< Calculate Exchange integrals without multipole approximation
  job->spin_dim          = job_parameters[25]; //!< Used for summation over spin in spin-polarised(= 2)/unpolarised(= 1) cases
  job->spin_fac          = job_parameters[26]; //!< Spin factor in spin-polarised ( =2)/unpolarised (= 1) cases
  job->spin_index        = job_parameters[27]; //!< Used to loop over spin up and down
  job->spin_pol          = job_parameters[28];
  job->spin_orb          = job_parameters[29]; //!< Spin orbit coupling switch - on iff set to 1  && spin_pol == 1 - off by default
  job->spin_polarisation = job_parameters[30]; //!< Switch (= 0) for spin-unpolarised cases (= 1) for spin-polarised cases
  job->mixing_type       = job_parameters[31]; //!< Parameter which determines type of mixing in SCF: 
  job->pms               = job_parameters[32]; //!< Switch (= 0/1) for permutation symmetry off/on
  job->trs               = job_parameters[33]; //!< Switch (= 0/1) for time reversal symmetry off/on
  job->sgs               = job_parameters[34]; //!< Switch (= 0/1) for space group symmetry off/on
  job->rss               = job_parameters[35]; //!< Switch (= 0/1) for real space symmetry off/on
  job->kss               = job_parameters[36]; //!< Switch (= 0/1) for reciprocal space symmetry off/on
  job->lmax              = job_parameters[37]; //!< Maximum l value for multipole moment expansion
  job->lmax_fac          = job_parameters[38]; //!< Maximum number of monomials for given l value plus 1 for spheropole moment
  job->C09               = job_parameters[39]; //!< Data read from Crystal09 if set to 1
  job->mxr               = job_parameters[40]; //!< Maximum side of grid for R->vec_ai real space lattice vectors 
  job->overlap_tol_1     = job_parameters1[0]; //!< Overlap criterion used to decide whether to retain pair_p pairs
  job->overlap_tol_2     = job_parameters1[1]; //!< Overlap criterion used to decide whether to retain pair_c pairs
  job->overlap_tol_3     = job_parameters1[2]; //!< Overlap criterion used to decide whether to retain pair_e pairs
  job->overlap_tol_4     = job_parameters1[3]; //!< Overlap criterion used to decide whether charges penetrate in point multipole approximation
  job->overlap_tol_5     = job_parameters1[4]; //!< Overlap criterion used to decide whether to retain pair_t lattice triples
  job->overlap_tol_6     = job_parameters1[5]; //!< Overlap criterion used to decide whether to retain pair_t lattice quads
  job->overlap_tol_7     = job_parameters1[6]; //!< Overlap criterion used to decide whether to retain pair_t point triples
  job->overlap_tol_8     = job_parameters1[7]; //!< Overlap criterion used to decide whether to retain pair_t point quads
  job->overlap_tol_9     = job_parameters1[8]; //!< Overlap criterion used to decide whether to retain pair_e quads
  job->electron_count    = job_parameters1[9]; //!< Total number of electrons per primitive unit cell
  job->ham_hyb_wt        = job_parameters1[10];//!< Fock exchange weight in hybrid DFT Hamiltonian
  //double fermi_energy ; //!< Fermi energy
  //double total_energy ; //!< Total SCF energy
  //double energy_change; //!< Total SCF energy change in an iteration
  //double nuc_nuc ;      //!< Nuclear-nuclear repulsion energy
  //double fock_mixing ;  //!< Mixing fraction in density matrix mixing between 0 (no damping) and 1
  //double scf_tol ;      //!< Parameter for SCF convergence on energy
  //double energy_range[2]; //!< Energy range associated with job->field_dirs
  //double total_time ;   //!< Accumulates total wall time expired
  //VECTOR_DOUBLE e_field[3];//!< Electric field direction cosines associated with job->field_dirs
  //clock_t start ;       //!< not used
  //clock_t end ;         //!< not used

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"read_JOB_PARAM_array\n");
  for (i = 0; i < 41; i++) {
  fprintf(file.out,"%3d %5d\n",i,job_parameters[i]);
 }
  for (i = 0; i < 11; i++) {
  fprintf(file.out,"%3d %10.4lf\n",i,job_parameters1[i]);
 }
 }

}

void write_JOB_PARAM_array(FERMI *fermi, int *job_parameters, double *job_parameters1, JOB_PARAM *job, FILES file)

{

int i;

  // ******************************************************************************************
  // * State of job is captured in job_parameters and job_parameters1 arrays                  *
  // * Parameters which are needed to recreate a particular state are written to disk         *
  // ******************************************************************************************

  //int monk_pack;        //!< not used
  //int verbosity ;       //!< Printing level during runs. Default value 1
  //int print_pairs ;     //!< Print pair_p, pair_c or pair_e when generated
  //int check_pairs ;     //!< Check pair_p, pair_c or pair_e when generated
  //int numtasks ;        //!< Number of tasks generated under MPI_WORLD_COMM
  //int taskid ;          //!< Number of current task
  //int *Fermi ;          //!< 
  //int max_cycle ;       //!< Maximum number of cycles in SCF calculation
  //int iter ;            //!< Current iteration in SCF calculation
  //int guess_type ;      //!< Type of guess used to start/restart SCF calculation
  //int mpi_io ;          //!< IO switch (= 0) use multiple local files (= 1) use single MPI files
  //int scf_direct ;      //!< Set to 1 if SCF integrals are recalculated each cycle, set to 0 if they are stored on disk
  //int vectors ;         //!< eigenvectors to be read from disk if set to 1 to be read or written if set to 2
  //int values ;          //!< eigenvalues to be read from disk if set to 1 to be read or written if set to 2
  //int density ;         //!< density matrix to be read from disk if set to 1 to be read or written if set to 2
  //int kpoints ;         //!< K point net type: Monkhorst-Pack = 0; list of points = 1
  //int field_dirs ;      //!< Number of field directions for susceptibilities
  //int npoints;          //!< Number of energy points at which a spectrum is to be sampled
  //int mixing_order ;    //!< Number of generations of density matrix used in Pulay mixing
  //int job->*memory ;    //!< Counts current memory allocated, node by node
  //int *max_memory ;     //!< Counts maximum memory allocated, node by node
  //int fix_occ ;         //!< Occupancy for density matrix calculation fixed if set to 1, determined by Fermi level if set to 0
  job_parameters[0]   = fermi->is[0];
  job_parameters[1]   = fermi->is[1];
  job_parameters[2]   = fermi->is[2];
  job_parameters[3]   = fermi->nkunique;
  job_parameters[4]   = job->dimp;             //!< dimension of reduced density matrix (multiply by 2 for spin-polarised density matrix)
  job_parameters[5]   = job->dimf;             //!< dimension of full density matrix (multiply by 2 for spin-polarised density matrix)
  job_parameters[6]   = job->type ;            //!< Job type  0 = scf, 1 = bse
  job_parameters[7]   = job->xc_num ;          //!< Number of functionals for xclib
  job_parameters[8]   = job->xc_typ[0] ;       //!< XC_FAMILY parameter for xclib
  job_parameters[9]   = job->xc_typ[1] ;       //!< XC_FAMILY parameter for xclib
  job_parameters[10]  = job->xc_hfx ;          //!< Hamiltonian requires Hartree-Fock exchange
  job_parameters[11]  = job->xc_lmx ;          //!< Max L value in Lededev grid for angular integration in DFT
  job_parameters[12]  = job->xc_sph ;          //!< Number of angular integration points in DFT
  job_parameters[13]  = job->xc_rad ;          //!< Number of radial integration points in DFT
  job_parameters[14]  = job->xc_grd ;          //!< Grid type in DFT
  job_parameters[15]  = job->ham_type ;        //!< Hamiltonian type: Set to 0 for DFT, 1 for hybrid DFT, 2 for Hartree-Fock
  job_parameters[16]  = job->ham_dft  ;        //!< DFT Hamiltonian type: Set to 0 for DFT, 1 for DFT-GGA
  job_parameters[17]  = job->ham_dft_num ;     //!< DFT Hamiltonian type: Number of density functionals needed - max value 2
  job_parameters[18]  = job->ham_dft_exc[0] ;  //!< DFT Hamiltonian type: Exchange functional list
  job_parameters[19]  = job->ham_dft_exc[1] ;  //!< DFT Hamiltonian type: Exchange functional list
  job_parameters[20]  = job->ham_dft_cor[0] ;  //!< DFT Hamiltonian type: Correlation functional list
  job_parameters[21]  = job->ham_dft_cor[1] ;  //!< DFT Hamiltonian type: Correlation functional list
  job_parameters[22]  = job->ham_hyb  ;        //!< Hybrid DFT Hamiltonian type: Set to 1 for hybrid DFT
  job_parameters[23]  = job->coul_int ;        //!< Calculate Coulomb integrals without multipole approximation
  job_parameters[24]  = job->exch_int ;        //!< Calculate Exchange integrals without multipole approximation
  job_parameters[25]  = job->spin_dim ;        //!< Used for summation over spin in spin-polarised(= 2)/unpolarised(= 1) cases
  job_parameters[26]  = job->spin_fac ;        //!< Spin factor in spin-polarised ( =2)/unpolarised (= 1) cases
  job_parameters[27]  = job->spin_index ;      //!< Used to loop over spin up and down
  job_parameters[28]  = job->spin_pol ;
  job_parameters[29]  = job->spin_orb ;        //!< Spin orbit coupling switch - on iff set to 1  && spin_pol == 1 - off by default
  job_parameters[30]  = job->spin_polarisation;//!< Switch (= 0) for spin-unpolarised cases (= 1) for spin-polarised cases
  job_parameters[31]  = job->mixing_type ;     //!< Parameter which determines type of mixing in SCF: 
  job_parameters[32]  = job->pms ;             //!< Switch (= 0/1) for permutation symmetry off/on
  job_parameters[33]  = job->trs ;             //!< Switch (= 0/1) for time reversal symmetry off/on
  job_parameters[34]  = job->sgs ;             //!< Switch (= 0/1) for space group symmetry off/on
  job_parameters[35]  = job->rss ;             //!< Switch (= 0/1) for real space symmetry off/on
  job_parameters[36]  = job->kss ;             //!< Switch (= 0/1) for reciprocal space symmetry off/on
  job_parameters[37]  = job->lmax ;            //!< Maximum l value for multipole moment expansion
  job_parameters[38]  = job->lmax_fac ;        //!< Maximum number of monomials for given l value plus 1 for spheropole moment
  job_parameters[39]  = job->C09;              //!< Data read from Crystal09 if set to 1
  job_parameters[40]  = job->mxr;              //!< Maximum side of grid for R->vec_ai real space lattice vectors 
  job_parameters1[0]  = job->overlap_tol_1; //!< Overlap criterion used to decide whether to retain pair_p pairs
  job_parameters1[1]  = job->overlap_tol_2; //!< Overlap criterion used to decide whether to retain pair_c pairs
  job_parameters1[2]  = job->overlap_tol_3; //!< Overlap criterion used to decide whether to retain pair_e pairs
  job_parameters1[3]  = job->overlap_tol_4; //!< Overlap criterion used to decide whether charges penetrate in point multipole approximation
  job_parameters1[4]  = job->overlap_tol_5; //!< Overlap criterion used to decide whether to retain pair_t lattice triples
  job_parameters1[5]  = job->overlap_tol_6; //!< Overlap criterion used to decide whether to retain pair_t lattice quads
  job_parameters1[6]  = job->overlap_tol_7; //!< Overlap criterion used to decide whether to retain pair_t point triples
  job_parameters1[7]  = job->overlap_tol_8; //!< Overlap criterion used to decide whether to retain pair_t point quads
  job_parameters1[8]  = job->overlap_tol_9; //!< Overlap criterion used to decide whether to retain pair_e quads
  job_parameters1[9]  = job->electron_count;//!< Total number of electrons per primitive unit cell
  job_parameters1[10] = job->ham_hyb_wt;    //!< Fock exchange weight in hybrid DFT Hamiltonian
  //double fermi_energy ; //!< Fermi energy
  //double total_energy ; //!< Total SCF energy
  //double energy_change; //!< Total SCF energy change in an iteration
  //double nuc_nuc ;      //!< Nuclear-nuclear repulsion energy
  //double fock_mixing ;  //!< Mixing fraction in density matrix mixing between 0 (no damping) and 1
  //double scf_tol ;      //!< Parameter for SCF convergence on energy
  //double energy_range[2]; //!< Energy range associated with job->field_dirs
  //double total_time ;   //!< Accumulates total wall time expired
  //VECTOR_DOUBLE e_field[3];//!< Electric field direction cosines associated with job->field_dirs
  //clock_t start ;       //!< not used
  //clock_t end ;         //!< not used

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"write_JOB_PARAM_array\n");
  for (i = 0; i < 41; i++) {
  fprintf(file.out,"%3d %5d\n",i,job_parameters[i]);
 }
  for (i = 0; i < 11; i++) {
  fprintf(file.out,"%3d %10.4lf\n",i,job_parameters1[i]);
 }
 }

}

void density_matrix_molecule2(FERMI *fermi, double *P0, double *F0, KPOINT_TRAN *knet, int *nkunique, REAL_LATTICE *R, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

  {

    int i, j, k, p, q, s;
    int dimp_spin, dimf_spin;
    int dim1, dim2, count;

    dimp_spin = job->dimp * job->spin_dim;
    dimf_spin = job->dimf * job->spin_dim;

    if (job->guess_type == 0 && job->density == 2)   // only zero P0 if beginning from atomic wave functions
      for (i=0;i<dimp_spin;i++) { P0[i] = k_zero;}
      for (i=0;i<dimf_spin;i++) { F0[i] = k_zero;}

    if (job->guess_type == 0 && job->density == 2) { 
    reduced_density_matrix_molecule2(fermi,P0,knet,nkunique,R,pair_p,atom_p,atoms,shells,symmetry,job,file);
   }

    expand_density_matrix(P0,F0,pair_p,atoms,shells,symmetry,job,file);

    if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"density matrix %d\n",pair_p->nump);
      count = 0;
      for (s = 0; s < job->spin_dim; s++) {
        for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
            fprintf(file.out,"pair %d [%3d %3d] spin %d \n",p,pair_p->cell1[q],pair_p->cell2[q],s);
              for (k = 0; k < pair_p->numb[p]; k++) {
          for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
            for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
              //fprintf(file.out,"%16.12f ",F0[count]);
              fprintf(file.out,"%6.2f ",F0[count]);
              count++;
              }
            fprintf(file.out,"\n");
           }
            fprintf(file.out,"\n");
           }
          fprintf(file.out,"\n\n");
         }
        fprintf(file.out,"\n\n");
       }
      }

}

void reduced_density_matrix_molecule2(FERMI *fermi, double *P, KPOINT_TRAN *knet, int *nkunique, REAL_LATTICE *R, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, l, n, p, q, s;
int dim, dim1, dimk;
int dimp_spin = job->dimp * job->spin_dim;
int bfposi, bfposj;
int vector_size;
int count;
size_t result;
double *P_buffer, *eigenvalues;
char zz[24] = "scf_evectors";
FILE *scf_evectors; 
ComplexMatrix *eigenvectors;
//char buf2[110];
//char yy[10] = "/scf_evec";;
//MPI_File fh;

    dim1 = atoms->number_of_sh_bfns_in_unit_cell;
    dimk = job->spin_dim * atoms->number_of_sh_bfns_in_unit_cell;
    //printf("dimp %3d %3d %3d\n",dimk,dim1,dimp_spin);
    vector_size = dimk * dim1;
    AllocateComplexMatrix(&eigenvectors,&dimk,&dim1,job);
    AllocateDoubleArray(&eigenvalues,&dimk,job);
    AllocateDoubleArray(&P_buffer,&dimp_spin,job); // check between these last two alternatives : second, seems correct but dumps
    ResetDoubleArray(P_buffer,&dimp_spin);

    //strcpy(buf2,file.scf_eigvec);
    //strcat(buf2,yy);
    //MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
    //MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    //MPI_File_read(fh, &eigenvectors->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    //MPI_File_close(&fh);
    //fprintf(file.out,"eigenvectors from disk\n");
    //PrintComplexMatrix(eigenvectors, file);
    //print_complex_eigenvector_matrix(eigenvectors, eigenvalues, file);

    if (job->taskid == 0) {

    scf_evectors = fopen(zz, "rb");
    result = fread(&eigenvectors->a[0][0], sizeof(double), 2 * vector_size, scf_evectors);
    fclose(scf_evectors);
    //print_complex_eigenvector_matrix(eigenvectors, eigenvalues, file);

    //printf("fermi->occupied density matrix2 %3d %3d\n",fermi->occupied[0],fermi->occupied[1]);

    dim = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (p = 0; p < pair_p->nump; p++) {
        q = pair_p->posn[p];
        bfposi = 0;
        bfposj = 0;
        for (l = 0; l < pair_p->cell1[q]; l++)
          bfposi += atoms->bfnnumb_sh[l];
          for (l = 0; l < pair_p->cell2[q]; l++)
            bfposj += atoms->bfnnumb_sh[l];
              for (n = 0; n < fermi->occupied[s]; n++) {
                count = 0;
                for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
                  for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
                    //fprintf(file.out,"s %3d p %3d %3d oc %3d i %3d j %3d dim %3d c %3d i %3d j %3d %16.11lf %16.11lf %16.11lf\n", \
                    s,p,l,n,i,j,dim,count,bfposi+i,bfposj+j,eigenvectors->a[s * dim1 + n][bfposi + i].real(),\
                    eigenvectors->a[s * dim1 + n][bfposj + j].real(),job->spin_fac*eigenvectors->a[s * dim1 + n][bfposi + i].real()*\
                    eigenvectors->a[s * dim1 + n][bfposj + j].real());
                    P_buffer[dim + count] += job->spin_fac * eigenvectors->a[s * dim1 + n][bfposi + i].real() * \
                    eigenvectors->a[s * dim1 + n][bfposj + j].real();
                    //P[dim + count] += job->spin_fac * eigenvectors->a[s * dim1 + n][bfposi + i].real() * \
                    eigenvectors->a[s * dim1 + n][bfposj + j].real();
                    count++;
                   }
                  }
                 }
               dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
              }
             //printf("%3d %3d\n",dim,count);
             }
            } // close if (job->taskid

            //MPI_Allreduce(P_buffer,P,vector_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(P_buffer,P,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            DestroyComplexMatrix(&eigenvectors,job);
            DestroyDoubleArray(&eigenvalues,&dimk,job);
            DestroyDoubleArray(&P_buffer,&dimp_spin,job);
            //free(eigenvalues);

   if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
            fprintf(file.out,"pair %d [%3d %3d] spin %d \n",p,pair_p->cell1[q],pair_p->cell2[q],s);
         //fprintf(file.out,"pair %d spin %d \n",p,s);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             //fprintf(file.out,"%16.12lf ",P[count]);
             fprintf(file.out,"%6.2f ",P[count]);
             //fprintf(file.out,"%9.2e ",P[count]);
             //fprintf(file.out,"%10.3e ",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
     }

}

void expand_density_matrix(double *P, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int count, dim1, dim2, i, j, p, q, r, s;

    for (s = 0; s < job->spin_dim; s++) {
      dim1 = s * job->dimp;
      dim2 = s * job->dimf;
      for (p = 0; p < pair_p->nump; p++) {
        q = pair_p->posn[p];
        rotate_permute_expand_pair(p, pair_p, &P[dim1], &F[dim2], atoms, shells, symmetry, job, file);
        dim1 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
        dim2 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]] * pair_p->numb[p];
       }
      }

   if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"full density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         for (r = 0; r < pair_p->numb[p]; r++) {
          fprintf(file.out,"pair %d [%3d %3d] gj %d \n",p,pair_p->cell1[q + r],pair_p->cell2[q + r],pair_p->latt2[q + r]);
          for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
            for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
              fprintf(file.out,"%16.12lf ",F[count]);
              //fprintf(file.out,"%9.2e ",P[count]);
              //fprintf(file.out,"%10.3e ",P[count]);
              count++;
             }
            fprintf(file.out,"\n");
           }
          fprintf(file.out,"\n");
         }
        }
       fprintf(file.out,"\n");
      }
     }

}

void atom_shell_populations2(INT_1E *one_ints, double *P, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, p, q, r, s, p1, count, count1;
int dim2, dim3, dim4;
int nd1, nd2;
int index_i, index_j, shelposi, shelposj, shlposi, sheli, shelj;
int i1, j1;
double *Overlap, *Density;

  for (i = 0; i < job->spin_dim * atoms->number_of_sh_shells_in_unit_cell; i++) {
    shells->pop_sh[i] = k_zero;
   }
  dim4 = 0;
  for (s = 0; s < job->spin_dim; s++) {
  dim2 = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q  = pair_p->posn[p];
      i1 = pair_p->cell1[q];
      j1 = pair_p->cell2[q];
      nd1 = atoms->bfnnumb_sh[i1];
      nd2 = atoms->bfnnumb_sh[j1];
      dim3 = atoms->bfnnumb_sh[i1] * atoms->bfnnumb_sh[j1] * pair_p->numb[p];
      AllocateDoubleArray(&Overlap,&dim3,job);
      AllocateDoubleArray(&Density,&dim3,job);
      ResetDoubleArray(Overlap,&dim3);
      ResetDoubleArray(Density,&dim3);
      rotate_permute_expand_pair(p, pair_p, &one_ints->Overlap[dim2], Overlap, atoms, shells, symmetry, job, file);
      rotate_permute_expand_pair(p, pair_p, &P[dim4], Density, atoms, shells, symmetry, job, file);

  if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"expanded pair density\n");
    count1 = 0;
    for (p1 = 0; p1 < pair_p->numb[p]; p1++) {
      fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p1,s,pair_p->cell1[q + p1],pair_p->cell2[q + p1],pair_p->latt2[q + p1]);
      for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
        for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
          fprintf(file.out,"%10.4lf",Density[count1]);
          count1++;
         }
        fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n");
     }
    }

      dim2 += nd1 * nd2;
      dim4 += nd1 * nd2;
      count = 0;
       for (r = 0; r < pair_p->numb[p]; r++) {
        i1 = pair_p->cell1[q + r];
        j1 = pair_p->cell2[q + r];
        shlposi = atoms->shlposn_sh[i1];
        shelposi = atoms->shelposn_sh[i1];
        shelposj = atoms->shelposn_sh[j1];
        //fprintf(file.out,"p %3d r %3d q %3d i1 %3d j1 %3d shlposi %3d shelposi %3d shelposj %3d\n",p,r,q,i1,j1,shlposi, \
        shelposi,shelposj);
        for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[i1]; index_i++) {
          sheli = shells->type_sh[index_i];
            for (i = 0; i < sheli; i++) {
              for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[j1]; index_j++) {
                shelj = shells->type_sh[index_j];
                  for (j = 0; j < shelj; j++) {
                    shells->pop_sh[s * atoms->number_of_sh_shells_in_unit_cell + shlposi] += Overlap[count] * Density[count] ;
                    //fprintf(file.out,"%3d %3d %3d %3d %3d %10.4lf %10.4lf\n",i,j,shlposi,shelposi,count,Overlap[count], \
                    Density[count]);
                      count++;
                     }
                    }
                   }
                  shlposi++;
                  }
                 }
                DestroyDoubleArray(&Overlap,&dim3,job);
                DestroyDoubleArray(&Density,&dim3,job);
               }
              }

}

double print_total_population2(ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

int i, j, s, count;
double total_population, total_population_up, total_population_down;

  count = 0;
  total_population      = k_zero;
  total_population_up   = k_zero;

  for (s = 0; s < job->spin_dim; s++) {
  total_population_down = k_zero;
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        for (j = 0; j < atoms->nshel_sh[i]; j++) {
          total_population      += shells->pop_sh[count];
          total_population_down += shells->pop_sh[count];
          //printf("%3d %3d %lf\n",i,j,shells->pop_sh[count]);
          count++;
         }
        }
       }

  total_population_up = total_population - total_population_down;

  //if (job->taskid == 0) 
  //printf("Total Population %15.9lf\n", total_population);
  //fprintf(file.out, "Total Population %15.9lf\n", total_population);
  //if (job->spin_polarisation == 1)
  //fprintf(file.out, "Total Population Up %15.9lf Total Population Down %15.9lf \n", total_population_up, total_population_down);
  //fflush(file.out);
 
  return total_population;

}

void print_atom_populations2(ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

int i, j, s, count, nrows, remdr;
int dim = atoms->number_of_atoms_in_unit_cell;
double total_pop = k_zero, total_spin = k_zero;

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        atoms->pop[s * atoms->number_of_atoms_in_unit_cell + i] = k_zero;
        for (j = 0; j < atoms->nshel_sh[i]; j++) {
          atoms->pop[s * atoms->number_of_atoms_in_unit_cell + i] += shells->pop_sh[count];
          count++;
         }
        }
       }

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        for (j = 0; j < atoms->nshel_sh[i]; j++) {
          total_pop += shells->pop_sh[count];
          count++;
         }
        }
       }

    count = 0;
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
      total_spin += atoms->pop[count] - atoms->pop[atoms->number_of_atoms_in_unit_cell + count];
      //fprintf(file.out,"%lf %lf %lf \n",total_spin, atoms->pop[count], atoms->pop[atoms->number_of_atoms_in_unit_cell + count]);
      count++;
     }

    nrows = atoms->number_of_atoms_in_unit_cell / 10;
    remdr = atoms->number_of_atoms_in_unit_cell - nrows * 10;

    if (job->taskid == 0) {

    if (job->spin_polarisation == 0) {
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|                                        ATOM POPULATION  %9.5lf                                       |\n",total_pop);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    for (i = 0; i < nrows; i++) {
      fprintf(file.out, "|   %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5f   |\n", \
      atoms->pop[i * 10 + 0],atoms->pop[i * 10 + 1], atoms->pop[i * 10 + 2],atoms->pop[i * 10 + 3], atoms->pop[i * 10 + 4], \
      atoms->pop[i * 10 + 5],atoms->pop[i * 10 + 6], atoms->pop[i * 10 + 7],atoms->pop[i * 10 + 8], atoms->pop[i * 10 + 9]);
     }
    if (remdr > 0) {
      fprintf(file.out,"|   ");
    for (i = 0; i < remdr; i++) 
      fprintf(file.out,"%9.5lf ",atoms->pop[nrows * 10 + i]);
    for (i = 0; i < 10 - remdr; i++) 
      fprintf(file.out,"          ");
      fprintf(file.out,"  |\n");
     }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"\n");
   }

  if (job->spin_polarisation == 1) {
     fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
     fprintf(file.out,"|                                        ATOM POPULATION  %9.5lf                                       |\n",total_pop);
     fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
     for (i = 0; i < nrows; i++) {
     fprintf(file.out, "|   %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5f   |\n", \
     atoms->pop[i * 10 + 0] + atoms->pop[dim + i * 10 + 0], atoms->pop[i * 10 + 1] + atoms->pop[dim + i * 10 + 1], \
     atoms->pop[i * 10 + 2] + atoms->pop[dim + i * 10 + 2], atoms->pop[i * 10 + 3] + atoms->pop[dim + i * 10 + 3], \
     atoms->pop[i * 10 + 4] + atoms->pop[dim + i * 10 + 4], atoms->pop[i * 10 + 5] + atoms->pop[dim + i * 10 + 5], \
     atoms->pop[i * 10 + 6] + atoms->pop[dim + i * 10 + 6], atoms->pop[i * 10 + 7] + atoms->pop[dim + i * 10 + 7], \
     atoms->pop[i * 10 + 8] + atoms->pop[dim + i * 10 + 8], atoms->pop[i * 10 + 9] + atoms->pop[dim + i * 10 + 9]);
     }
    if (remdr > 0) {
      fprintf(file.out,"|   ");
    for (i = 0; i < remdr; i++) 
      fprintf(file.out,"%9.5lf ",atoms->pop[nrows * 10 + i] + atoms->pop[dim + nrows * 10 + i]);
    for (i = 0; i < 10 - remdr; i++) 
      fprintf(file.out,"          ");
      fprintf(file.out,"  |\n");
     }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|                                        SPIN POPULATION  %9.5lf                                       |\n",total_spin);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    for (i = 0; i < nrows; i++) {
    fprintf(file.out, "|   %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5f   |\n", \
    atoms->pop[i * 10 + 0] - atoms->pop[dim + i * 10 + 0], atoms->pop[i * 10 + 1] - atoms->pop[dim + i * 10 + 1], \
    atoms->pop[i * 10 + 2] - atoms->pop[dim + i * 10 + 2], atoms->pop[i * 10 + 3] - atoms->pop[dim + i * 10 + 3], \
    atoms->pop[i * 10 + 4] - atoms->pop[dim + i * 10 + 4], atoms->pop[i * 10 + 5] - atoms->pop[dim + i * 10 + 5], \
    atoms->pop[i * 10 + 6] - atoms->pop[dim + i * 10 + 6], atoms->pop[i * 10 + 7] - atoms->pop[dim + i * 10 + 7], \
    atoms->pop[i * 10 + 8] - atoms->pop[dim + i * 10 + 8], atoms->pop[i * 10 + 9] - atoms->pop[dim + i * 10 + 9]);
     }
    if (remdr > 0) {
      fprintf(file.out,"|   ");
    for (i = 0; i < remdr; i++) 
      fprintf(file.out,"%9.5lf ",atoms->pop[nrows * 10 + i] - atoms->pop[dim + nrows * 10 + i]);
    for (i = 0; i < 10 - remdr; i++) 
      fprintf(file.out,"          ");
      fprintf(file.out,"  |\n");
     }
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"\n");
     }

      fflush(file.out);
    }
}

void print_shell_populations2(ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

int i, j, s, nrows, remdr, offset;
int count;

int dim = atoms->number_of_atoms_in_unit_cell;
double total_pop = k_zero, total_spin = k_zero;

  if (job->taskid == 0) {

    fprintf(file.out,"Shell Population\n");
    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        for (j = 0; j < atoms->nshel_sh[i]; j++) {
          fprintf(file.out,"%10.6lf",shells->pop_sh[count]);
          count++;
         }
        fprintf(file.out,"\n");
        }
       }

      }

}

void init_mix(DoubleMatrix *C_i, DoubleMatrix *C_o, DoubleMatrix *D, double *P0, double *P1)

{

  int i, j;

  for (i = C_i->iRows - 1; i > 0 ; i--) {
   for (j = 0; j < C_i->iCols; j++) {
     C_i->a[i][j] = C_i->a[i - 1][j];
     C_o->a[i][j] = C_o->a[i - 1][j];
     D->a[i][j]   =   D->a[i - 1][j];
    }
   }

   for (j = 0; j < C_i->iCols; j++) {
     C_i->a[0][j] = P0[j];
     C_o->a[0][j] = P1[j];
     D->a[0][j]   =     C_o->a[0][j] - C_i->a[0][j];
    }

}

void mix_density(DoubleMatrix *C_i, DoubleMatrix *C_o, DoubleMatrix *D, double *P, JOB_PARAM *job, FILES file)

{

  int i, j, k, nmix, nmix1, ione = 1;
  int *work;
  int info;
  double *beta, *d_i, *d_o;
  DoubleMatrix *C;

  nmix = (job->mixing_order < job->iter) ? job->mixing_order : job->iter;
  //nmix = (*s < job->iter) ? *s : job->iter;
  nmix--;
  nmix1 = nmix + 1;
  //printf("nmix %d job->iter %d mixing_order %d\n",nmix,job->iter,job->mixing_order);

  if (nmix == 0) {
    for (i = 0; i < C_i->iCols; i++) {
     P[i] = job->fock_mixing * C_i->a[0][i] + (k_one - job->fock_mixing) * C_o->a[0][i];
     }
    }

 else {
  AllocateDoubleMatrix(&C, &nmix1, &nmix1,job);
  AllocateDoubleArray(&beta, &nmix1, job);
  AllocateDoubleArray(&d_i, &(C_i->iCols), job);
  AllocateDoubleArray(&d_o, &(C_i->iCols), job);
  work = (int *) malloc(nmix1 * sizeof(int));
  ResetDoubleArray(d_i, &(C_i->iCols));
  ResetDoubleArray(d_o, &(C_i->iCols));
  ResetDoubleArray(beta, &nmix);

    for (i = 0; i < nmix; i++) {
     for (j = 0; j < nmix; j++) {
       DoubleDotProd(&C->a[i][j], &(C_i->iCols), &D->a[i], &ione, &D->a[j], &ione);
      }
     }

     beta[nmix] = -k_one;
      for (i = 0; i < nmix; i++) {
        C->a[i][nmix] = -k_one;
        C->a[nmix][i] = -k_one;
        work[i] = 0;
       }
        C->a[nmix][nmix] = k_zero;
        work[nmix] = 0;

    //fprintf(file.out,"C\n");
    //print_double_matrix(C, file);

    //for (i=0;i<=nmix;i++) {
    //fprintf(file.out,"beta[i] %e nmix %d\n",beta[i],nmix); }
    //fprintf(file.out,"\n");

    dgesv_(&nmix1,&ione,C->a[0],&nmix1,work,beta,&nmix1,&info);

    //for (i=0;i<=nmix;i++) {
    //fprintf(file.out,"beta[i] %e\n",beta[i]); }
    //fprintf(file.out,"\n");

    if (info == 0) {
      for (i = 0;i < nmix; i++) {
        for (j = 0;j < C_i->iCols; j++) {
          d_i[j] += beta[i] * C_i->a[i][j];
          d_o[j] += beta[i] * C_o->a[i][j];
         }
        }
         for (i = 0; i < C_i->iCols; i++) {
           P[i] = job->fock_mixing * d_i[i] + (k_one - job->fock_mixing) * d_o[i];
         }
        } // close if (info == 0)

    //for (i=0;i<C_i->iCols;i++) {
    //fprintf(file.out,"d_o[i]d_i[j] %8.3lf %8.3lf C_o %8.3lf %8.3lf C_i %8.3lf %8.3lf\n",d_o[i],d_i[i],C_o->a[1][i],\
    C_o->a[0][i],C_i->a[1][i],C_i->a[0][i]);}
    //fprintf(file.out,"\n");


        DestroyDoubleMatrix(&C, job);
        DestroyDoubleArray(&beta, &nmix1, job);
        DestroyDoubleArray(&d_i, &(C_i->iCols), job);
        DestroyDoubleArray(&d_o, &(C_i->iCols), job);

       } // close else

}

/* 
void density_matrix_crystal2(FERMI *fermi, double *P0, double *F0, KPOINT_TRAN *knet, char *filename, int *nkunique, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

  {

    int i, j, k, p, q, s;
    int dimp_spin, dimf_spin;
    int dim1, dim2, count;

    dimp_spin = job->dimp * job->spin_dim;
    dimf_spin = job->dimf * job->spin_dim;

    if (job->guess_type == 0 && job->density == 2)   // only zero P0 if beginning from atomic wave functions
      for (i=0;i<dimp_spin;i++) { P0[i] = k_zero;}
      for (i=0;i<dimf_spin;i++) { F0[i] = k_zero;}


    if (job->guess_type == 0 && job->density == 2) { 
    reduced_density_matrix_crystal(fermi,P0,knet,filename,R,R_tables,pair_p,atom_p,atoms,shells,symmetry,job,file);
   }

    expand_density_matrix(P0,F0,pair_p,atoms,shells,symmetry,job,file);

    if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"density matrix %d\n\n",pair_p->nump);
      count = 0;
      for (s = 0; s < job->spin_dim; s++) {
        for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
            fprintf(file.out,"unique pair %d spin %d\n",p,s);
              for (k = 0; k < pair_p->numb[p]; k++) {
              fprintf(file.out,"%3d   [%d] [%d]   [%d]\n",q+k,pair_p->cell1[q+k],pair_p->cell2[q+k],pair_p->latt2[q+k]);
          for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
            for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
              //fprintf(file.out,"%16.12lf ",F0[count]);
              //fprintf(file.out,"%10.3e ",F0[count]);
              fprintf(file.out,"%6.3f ",F0[count]);
              count++;
              }
            fprintf(file.out,"\n");
           }
            fprintf(file.out,"\n");
           }
          fprintf(file.out,"\n\n");
         }
        fprintf(file.out,"\n\n");
       }
      }

}

void reduced_density_matrix_crystal2(FERMI *fermi, double *P, KPOINT_TRAN *knet, int *nkunique, REAL_LATTICE *R, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, n, p, q, r, s;
int dim, dim1, dimk;
int bfposi, bfposj;
int vector_size;
int count, countk, dimp_spin = job->dimp * job->spin_dim;
double *eigenvalues;
char buf2[110];
char yy[10] = "/scf_evec";;
//char yy[10] = "/datafile";;
Complex *Pi;
ComplexMatrix *eigenvectors;
MPI_File fh;

    strcpy(buf2,file.scf_eigvec);
    strcat(buf2,yy);

    dim1 = atoms->number_of_sh_bfns_in_unit_cell;
    dimk = job->spin_dim * *nkunique * atoms->number_of_sh_bfns_in_unit_cell;
    vector_size = dimk * dim1;
    AllocateComplexMatrix(&eigenvectors,&dimk,&dim1,job);
    AllocateDoubleArray(&eigenvalues,&dimk,job);

    MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
    MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    MPI_File_read(fh, &eigenvectors->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    MPI_File_close(&fh);

    //fprintf(file.out,"eigenvectors from disk\n");
    //PrintComplexMatrix(eigenvectors, file);
    //print_complex_eigenvector_matrix(eigenvectors, eigenvalues, file);

    Pi = (Complex *) calloc(dimp_spin, sizeof(Complex));
    if (Pi == NULL) {  fprintf(stderr, "ERROR: not enough memory for Complex Pi! \n"); exit(1); }

    dim = 0;
    //2019for (s = 1; s < 2; s++) {
    for (s = 0; s < job->spin_dim; s++) {
      for (p = 0; p < pair_p->nump; p++) {
        q = pair_p->posn[p];
        r = pair_p->latt2[q];
        bfposi = 0;
        bfposj = 0;
        for (l = 0; l < pair_p->cell1[q]; l++)
          bfposi += atoms->bfnnumb_sh[l];
          for (l = 0; l < pair_p->cell2[q]; l++)
            bfposj += atoms->bfnnumb_sh[l];
            countk = 0;
            for (k = 0; k < *nkunique; k++) {
              for (m = 0; m < dimp_spin; m++) 
                Pi[m] = Complex(k_zero,k_zero);
                //temp = -double_vec_dot(&knet->cart[knet->ibz[k]],&R->vec_ai[r]) ;
                //temp1 = (double)job->spin_fac * Complex(cos(temp),sin(temp)) ;
                //fprintf(file.out,"temp temp1 %d %d %d %lf %lf %lf\n",p,q,k,temp,temp1.real(),temp1.imag());
                for (n = 0; n < fermi->occupied[s * *nkunique + k]; n++) {
                  count = 0;
                  for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
                    for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
  //fprintf(file.out,"pair %3d spin %3d kpt %3d %3d occn %3d i %3d j %3d dim %3d count %3d i %3d j %3d %3d %3d %10.4lf %10.4lf %10.4lf\n",\
  p,s,knet->ibz[k],k,n,i,j,\
  dim,count,bfposi+i,bfposj+j,fermi->nktot,fermi->knet->num[k],eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposi + i].real(),\
  eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposj + j].real(),(conj(eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposi + i]) * \
  eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposj + j]).real() / (double)fermi->nktot);
                            //P[dim + count] += (conj(eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposi + i]) * \
                           eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposj + j] * temp1).real() * fermi->knet->num[k] / fermi->nktot;
                      Pi[count] += conj(eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposi + i]) * \
                      eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposj + j] / (double)fermi->nktot;
                      //eigenvectors->a[(s * *nkunique + k) * dim1 + n][bfposj + j] * knet->weight[k];
                      count++;
                     }
                    }
                   }
                  rotate_reduced_density2(Pi, &P[dim], pair_p, &p, knet, &k, &countk, R, atom_p, atoms, shell, symmetry, job, file);
                 } // close loop on k
                dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
               }
              }

  DestroyComplexMatrix(&eigenvectors,job);
  DestroyDoubleArray(&eigenvalues,&dimk,job);

   if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         fprintf(file.out,"pair %d spin %d \n",p,s);
         q = pair_p->posn[p];
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%16.12lf ",P[count]);
             //fprintf(file.out,"%9.2e ",P[count]);
             //fprintf(file.out,"%10.3e ",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
     }

}

void reduced_density_matrix_crystal(FERMI *fermi, double *P, KPOINT_TRAN *knet, char *filename, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, gj, k, kp, l, m, n, p, q, s;
int dim, dim1, dim2, dimk;
int bfposi, bfposj;
int vector_size;
int count, count1, countk, kunique, dimp_spin = job->dimp * job->spin_dim;
int nbands, bands[4];
double temp;
Complex temp1;
char buf2[110];
ComplexMatrix *eigenvectors, *eigenvectors1;
MPI_File fh;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,filename);

  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  if (job->fix_occ == 0) {
  dim2 = dim1;
 }
  else if (job->fix_occ == 1) {
  s = 0; // need to fix this
  bands[0] = fermi->bands[2 * s] - 1; 
  bands[1] = fermi->bands[2 * s + 1] - 1;
  nbands = bands[1] - bands[0] + 1;
  dim2 = nbands;
 }
  dimk = job->spin_dim * fermi->nkunique * dim2;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&eigenvectors,&dimk,&dim1,job);
  AllocateComplexMatrix(&eigenvectors1,&dimk,&dim1,job);

  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_read(fh, &eigenvectors->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  MPI_File_close(&fh);

  count1 = 0;
  for (s = 0; s < job->spin_dim; s++) {
    countk = 0;
    for (k = 0; k < fermi->nkunique; k++) {
      if (job->kss == 0) kunique = countk;
      else if (job->kss == 1) kunique = k;
      for (kp = 0; kp < knet->num[k]; kp++) {
        if (job->fix_occ == 0) {
        bands[0] = 0; 
        bands[1] = fermi->occupied[s * fermi->nkunique + k];
        nbands = bands[1] - bands[0];
        dim2 = dim1;
       }
        else if (job->fix_occ == 1) {
        bands[0] = fermi->bands[2 * s] - 1; 
        bands[1] = fermi->bands[2 * s + 1] - 1;
        nbands = bands[1] - bands[0] + 1;
        dim2 = nbands;
       }
        rotate_psi(&eigenvectors->a[(s * fermi->nkunique + kunique) * dim2][0],&eigenvectors1->a[0][0],nbands,knet->bz[countk],knet,\
        atom_p,atoms,R,shell,symmetry,job,file);
        dim = s * job->dimp;
        for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
          gj = pair_p->latt2[q]; 
          temp = double_vec_dot(&knet->cart[knet->bz[countk]],&R->vec_ai[gj]) ;
          temp1 = (double)job->spin_fac * Complex(cos(temp),sin(temp)) ;
          bfposi = 0;
          bfposj = 0;
          for (l = 0; l < pair_p->cell1[q]; l++)
          bfposi += atoms->bfnnumb_sh[l];
          for (l = 0; l < pair_p->cell2[q]; l++)
          bfposj += atoms->bfnnumb_sh[l];
           for (n = 0; n < nbands; n++) {
            count = 0;
            for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
              for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
                P[dim + count] += (conj(eigenvectors1->a[n][bfposi + i]) * eigenvectors1->a[n][bfposj + j] * \
                fermi->occupation[count1 + n] / (double)fermi->nktot * temp1).real();
                count++;
               }
              }
             }
            dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
           } // close loop on p
          countk++;
         } // close loop on kp
        count1 += dim1;
       } // close loop on k
      } // close loop on s

  DestroyComplexMatrix(&eigenvectors,job);
  DestroyComplexMatrix(&eigenvectors1,job);

  if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p,s,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%7.3lf",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
       fflush(file.out);
     }

}

void reduced_density_matrix_crystal3(FERMI *fermi, double *P, KPOINT_TRAN *knet, int s, int k, ComplexMatrix *eigenvectors, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, gj, kp, l, n, p, q;
int dim, dim1;
int bfposi, bfposj;
int count, count1, countk, dimp_spin = job->dimp * job->spin_dim;
int nbands;
double temp;
Complex temp1;
ComplexMatrix *eigenvectors1;

  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  AllocateComplexMatrix(&eigenvectors1,&dim1,&dim1,job);

  countk = 0;
  for (kp = 0; kp < k; kp++) 
  countk += knet->num[kp];
  count1 = (s * fermi->nkunique + k ) * dim1;
  nbands = fermi->bands[1] - fermi->bands[0] + 1;
  for (kp = 0; kp < knet->num[k]; kp++) {
  //fprintf(file.out,"s %d k %d kp %d countk %d count1 %d nbands %d \n",s,k,kp,countk,count1,nbands);
      rotate_psi(&eigenvectors->a[0][0],&eigenvectors1->a[0][0],nbands,knet->bz[countk],knet,atom_p,atoms,R,shell,symmetry,job,file);
      dim = s * job->dimp;
      for (p = 0; p < pair_p->nump; p++) {
        q = pair_p->posn[p];
        gj = pair_p->latt2[q]; 
        temp = double_vec_dot(&knet->cart[knet->bz[countk]],&R->vec_ai[gj]) ;
        temp1 = (double)job->spin_fac * Complex(cos(temp),sin(temp)) ;
        bfposi = 0;
        bfposj = 0;
        for (l = 0; l < pair_p->cell1[q]; l++)
        bfposi += atoms->bfnnumb_sh[l];
        for (l = 0; l < pair_p->cell2[q]; l++)
        bfposj += atoms->bfnnumb_sh[l];
         for (n = 0; n < nbands; n++) {
          count = 0;
          for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
            for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
              P[dim + count] += (conj(eigenvectors1->a[n][bfposi + i]) * eigenvectors1->a[n][bfposj + j] * fermi->occupation[count1 + n] \
              / (double)fermi->nktot * temp1).real();
              count++;
             }
            }
           }
          dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
         } // close loop on p
       countk++;
       } // close loop on kp

  DestroyComplexMatrix(&eigenvectors1,job);

  if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"MPP reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p,s,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%7.3lf",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
       fflush(file.out);
     }

}

void reduced_density_matrix_crystal1(FERMI *fermi, double *P, KPOINT_TRAN *knet, char *filename, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, gj, k, kp, l, m, n, p, q, s;
int dim, dim1, dim2, dimk;
int bfposi, bfposj;
int vector_size;
int count, countk, kunique, dimp_spin = job->dimp * job->spin_dim;
int nbands, bands[4];
double sign, temp;
Complex temp1;
double *eigenvalues;
char buf2[110];
//char yy[10] = "/datafile";;
//char yy[10] = "/scf_evec";;
ComplexMatrix *eigenvectors, *eigenvectors1;
MPI_File fh;

  //strcat(buf2,yy);
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,filename);
  //printf("%s\n",buf2);

  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  if (job->fix_occ == 0) 
  dim2 = dim1;
  else if (job->fix_occ == 1) {
  s = 0; // need to fix this
  bands[0] = fermi->bands[2 * s] - 1; 
  bands[1] = fermi->bands[2 * s + 1] - 1;
  nbands = bands[1] - bands[0] + 1;
  dim2 = nbands;
 }
  dimk = job->spin_dim * fermi->nkunique * dim2;
  ////dimk = job->spin_dim * fermi->nkunique * atoms->number_of_sh_bfns_in_unit_cell;
 

  vector_size = dimk * dim1;
//printf("%3d %3d %3d\n",dimk,dim1,vector_size);
  AllocateComplexMatrix(&eigenvectors,&dimk,&dim1,job);
  AllocateComplexMatrix(&eigenvectors1,&dimk,&dim1,job);
  AllocateDoubleArray(&eigenvalues,&dimk,job);

  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_read(fh, &eigenvectors->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  MPI_File_close(&fh);

//for(j=0;j<dimk;j++) {
//for(i=0;i<dim1;i++) {
//fprintf(file.out,"%3d %3d %lf %lf\n",j,i,(eigenvectors->a[j][i]).real(),(eigenvectors->a[j][i]).imag());
//}}

  //fprintf(file.out,"fermi %s %4d %4d %4d %4d\n",buf2,fermi->nkunique,job->fix_occ,bands[0],bands[1]);
  //fprintf(file.out,"eigenvectors from disk  %5d %5d %5d \n",dimk, dim1, 2 * vector_size);
  //PrintComplexMatrix(eigenvectors, file);
  //if (job->taskid == 0)
  //print_complex_eigenvector_matrix(eigenvectors, eigenvalues, file);

  countk = 0;
  for (k = 0; k < fermi->nkunique; k++) {
         if (job->kss == 0) kunique = countk;
    else if (job->kss == 1) kunique = k;
    for (kp = 0; kp < knet->num[k]; kp++) {
      for (s = 0; s < job->spin_dim; s++) {
      //fprintf(file.out,"nkunique %3d k %3d kp %3d num %3d bz %3d\n",fermi->nkunique,k,kp,knet->num[k],knet->bz[countk]);
      // only doing spin down here
      ////for (s = 1; s < 2; s++) {
        //printf("k %3d kp %3d kuniq %3d bz %3d %3d %5d\n",k,kp,kunique,knet->bz[countk],fermi->nkunique,(s * fermi->nkunique + kunique) * dim1);
        //fprintf(file.out,"k %3d kp %3d kuniq %3d bz %3d %3d %5d\n",k,kp,kunique,knet->bz[countk],fermi->nkunique,(s * fermi->nkunique + kunique) * dim1);
        if (job->fix_occ == 0) {
        bands[0] = 0; 
        bands[1] = fermi->occupied[s * fermi->nkunique + k];
        nbands = bands[1] - bands[0];
dim2 = dim1;
       }
        else if (job->fix_occ == 1) {
        bands[0] = fermi->bands[2 * s] - 1; 
        bands[1] = fermi->bands[2 * s + 1] - 1;
        nbands = bands[1] - bands[0] + 1;
dim2 = nbands;
       }
        //fprintf(file.out,"fix_occ %3d bands[0] %3d bands[1] %3d nbands %3d\n",job->fix_occ,bands[0],bands[1],nbands);
        ////rotate_psi(&eigenvectors->a[(s * fermi->nkunique + kunique) * dim1 + bands[0]][0],&eigenvectors1->a[0][0],bands,knet->bz[countk],knet,atom_p, \
        atoms,R,shell,symmetry,job,file);
        //CHANGE2014rotate_psi(&eigenvectors->a[(s*fermi->nkunique+kunique)*dim2 + bands[0]][0],&eigenvectors1->a[0][0],bands,knet->bz[countk],knet,atom_p, \
        atoms,R,shell,symmetry,job,file);
        rotate_psi(&eigenvectors->a[(s * fermi->nkunique + kunique) * dim2 + bands[0]][0],&eigenvectors1->a[0][0],nbands,knet->bz[countk],knet,atom_p, \
        atoms,R,shell,symmetry,job,file);
        //printf("bands %3d    %3d %3d  offset %7d    %3d %3d %3d\n", \
        nbands,bands[0] + 1,bands[1] + 1,(s * fermi->nkunique + kunique) * dim1 + bands[0],k,kp,s);
        //fprintf(file.out,"bands %3d    %3d %3d  offset %7d    %3d %3d %3d\n", \
        nbands,bands[0] + 1,bands[1] + 1,(s * fermi->nkunique + kunique) * dim1 + bands[0],k,kp,s);
        //rotate_psi(&eigenvectors->a[(s * fermi->nkunique + kunique) * dim1][0],&eigenvectors1->a[0][0],bands,knet->bz[countk],knet,atom_p,atoms,R, \
        shell,symmetry,job,file);
        //if      (knet->trs[knet->bz[countk]] == 1) sign = -k_one;           // trs from time reversal symmetry
        //else if (knet->trs[knet->bz[countk]] == 0) sign =  k_one;
sign = k_one;
        dim = s * job->dimp;
        //dim = 0;
        //*
        for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
          gj = pair_p->latt2[q]; 
          temp = sign * double_vec_dot(&knet->cart[knet->bz[countk]],&R->vec_ai[gj]) ;
          // change 3 temp = -sign * double_vec_dot(&knet->cart[knet->bz[countk]],&R->vec_ai[gj]) ;
          temp1 = (double)job->spin_fac * Complex(cos(temp),sin(temp)) ;
          bfposi = 0;
          bfposj = 0;
          for (l = 0; l < pair_p->cell1[q]; l++)
          bfposi += atoms->bfnnumb_sh[l];
          for (l = 0; l < pair_p->cell2[q]; l++)
          bfposj += atoms->bfnnumb_sh[l];
          //for (n = 0; n < 1; n++) {
          // only doing 1st band here
           for (n = 0; n < nbands; n++) {
          //for (n = 0; n < fermi->occupied[s * fermi->nkunique + k]; n++) {
            count = 0;
            for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
              for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
                //P[dim + count] += (conj(eigenvectors->a[(s * fermi->nkunique + countk) * dim1 + n][bfposi + i]) * \
                eigenvectors->a[(s * fermi->nkunique + countk) * dim1 + n][bfposj + j] / (double)fermi->nktot * temp1).real();
                P[dim + count] += (conj(eigenvectors1->a[n][bfposi + i]) * \
                eigenvectors1->a[n][bfposj + j] / (double)fermi->nktot * temp1).real();
                count++;
               }
              }
             }
            dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
           } // close loop on p
          } // close loop on s
         countk++;
        } // close loop on kp
       } // close loop on k

  DestroyComplexMatrix(&eigenvectors,job);
  DestroyComplexMatrix(&eigenvectors1,job);
  DestroyDoubleArray(&eigenvalues,&dimk,job);

  //if (job->verbosity > 1) {
  if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p,s,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%7.3lf",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
       fflush(file.out);
     }

}
*/

/*
void expand_density_matrix_complex(Complex *P, Complex *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int count, dim1, dim2, i, j, p, q, r, s;

    for (s = 0; s < job->spin_dim; s++) {
      dim1 = s * job->dimp;
      dim2 = s * job->dimf;
      for (p = 0; p < pair_p->nump; p++) {
        q = pair_p->posn[p];
        rotate_permute_expand_pair_complex(p, pair_p, &P[dim1], &F[dim2], atoms, shells, symmetry, job, file);
        dim1 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
        dim2 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]] * pair_p->numb[p];
       }
      }

   if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"full matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         for (r = 0; r < pair_p->numb[p]; r++) {
           fprintf(file.out,"pair %d [%3d %3d] gj %d \n",p,pair_p->cell1[q + r],pair_p->cell2[q + r],pair_p->latt2[q + r]);
           for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
             for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
               fprintf(file.out,"%9.2e %9.2e ",(F[count]).real(), (F[count]).imag());
               count++;
              }
             fprintf(file.out,"\n");
            }
           fprintf(file.out,"\n");
          }
         }
        fprintf(file.out,"\n");
       }
      }
}
*/
/*
void expand_screening_integral_matrix_complex(Complex *P, Complex *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int dim1, dim2, count, i, j, p, q, r, s;
  
  dim1 = 0;
  dim2 = 0;
  for (p = 0; p < pair_p->nump; p++) {
    q = pair_p->posn[p];
    rotate_permute_expand_pair_complex(p, pair_p, &P[dim1], &F[dim2], atoms, shells, symmetry, job, file);
    dim1 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
    dim2 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]] * pair_p->numb[p];
   }

  if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"full screening integral matrix\n");
    count = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      for (r = 0; r < pair_p->numb[p]; r++) {
        fprintf(file.out,"pair %d [%3d %3d] gj %d \n",p,pair_p->cell1[q + r],pair_p->cell2[q + r],pair_p->latt2[q + r]);
        for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
          for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
          fprintf(file.out,"%5.2e %5.2e",(F[count]).real(),(F[count]).imag());
          count++;
         }
        fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n");
     }
    }
   fprintf(file.out,"\n");
  }

}

double normalise_density_matrix2(double *P, double *F, double total_population, JOB_PARAM *job, FILES file)

{

int i;
double fac;
 
    fac = total_population / job->electron_count;

    for (i = 0; i < job->spin_dim * job->dimp; i++) {
      P[i] /= fac;
   }

    for (i = 0; i < job->spin_dim * job->dimf; i++) {
      F[i] /= fac;
   }

   return fac;

}
*/
/*

void count_density_matrix_shells2(int *nshells, PAIR_TRAN *pair_c, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int p, q;
int ip, jp;

  *nshells = 0;
    for (p = 0; p < pair_c->nump; p++) {
      q  = pair_c->posn[p];
        ip = pair_c->cell1[q];
        jp = pair_c->cell2[q];
        *nshells += atoms->nshel_sh[ip] * atoms->nshel_sh[jp];
       } 
      *nshells *= job->spin_dim;

}

void density_matrix_largest_shell_elements2(double *lse, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_c, ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

int c, i, j, q;
int ip, jp, gj;
int index_i, index_j, shelposi, shelposj;
int count, count1;
int F_offset;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
double lse_max, lse_tmp;

  count1 = 0;
  for (c = 0; c < pair_c->nump; c++) {
    q  = pair_c->posn[c];
    ip = pair_c->cell1[q];
    jp = pair_c->cell2[q];
    gj = pair_c->latt2[q];
    F_offset = pair_p->off[pair_p->ptr[dim2 * gj + dim1 * ip + jp]];
    count = 0;
    shelposi = atoms->shelposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      for (i = 0; i < shells->type_sh[index_i]; i++) {
        shelposj = atoms->shelposn_sh[jp];
        for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
          lse_max = k_zero;
          for (j = 0; j < shells->type_sh[index_j]; j++) {
            lse_tmp = fabs(F[F_offset + count]);
            lse_max = (lse_tmp > lse_max ? lse_tmp : lse_max); 
            count++;
           }
          lse[count1] = lse_max; 
          count1++;
         }
        }
       }
      }

}

void orbital_product_populations(INT_1E *one_ints, double *P, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, p, q, r, s, count, count1;
int dim2, dim3, dim4;
int nd1, nd2;
int index_i, index_j, shelposi, shelposj, shlposi, sheli, shelj;
int i1, j1;
double *Overlap, *Density;

int numq = 3;
int g1, offset;
int nbfn = atoms->number_of_sh_shells_in_unit_cell;
int dimij = nbfn * nbfn * numq;
Complex *N_ij;
double q_dot_g;
AllocateComplexArray(&N_ij,&dimij,job);
VECTOR_DOUBLE qvec[3];
qvec[0].comp1 = 0.0;
qvec[0].comp2 = 0.0;
qvec[0].comp3 = 0.0;

  dim4 = 0;
  for (s = 0; s < job->spin_dim; s++) {
  dim2 = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q  = pair_p->posn[p];
      i1 = pair_p->cell1[q];
      j1 = pair_p->cell2[q];
      nd1 = atoms->bfnnumb_sh[i1];
      nd2 = atoms->bfnnumb_sh[j1];
      dim3 = atoms->bfnnumb_sh[i1] * atoms->bfnnumb_sh[j1] * pair_p->numb[p];
      AllocateDoubleArray(&Overlap,&dim3,job);
      AllocateDoubleArray(&Density,&dim3,job);
      ResetDoubleArray(Overlap,&dim3);
      ResetDoubleArray(Density,&dim3);
      rotate_permute_expand_pair(p, pair_p, &one_ints->Overlap[dim2], Overlap, atoms, shells, symmetry, job, file);
      rotate_permute_expand_pair(p, pair_p, &P[dim4], Density, atoms, shells, symmetry, job, file);

  if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"expanded density matrix\n");
     count1 = 0;
     //for (s = 0; s < job->spin_dim; s++) {
       //for (p = 0; p < pair_p->nump; p++) {
         //q = pair_p->posn[p];
       for (r = 0; r < pair_p->numb[p]; r++) {
         fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p,s,pair_p->cell1[q + r],pair_p->cell2[q + r],pair_p->latt2[q + r]);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%10.4lf",Density[count1]);
             count1++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
        //}
       //fprintf(file.out,"\n");
      //}
       //fflush(file.out);
     }

  dim2 += nd1 * nd2;
  dim4 += nd1 * nd2;
  offset = 0;
  for (q = 0; q < numq; q++) {
   count = 0;
   for (r = 0; r < pair_p->numb[p]; r++) {
    i1 = pair_p->cell1[q + r];
    j1 = pair_p->cell2[q + r];
    g1 = pair_p->latt2[q + r];
    //q_dot_g = double_vec_dot(&qvec[q],&R->vec_ai[r]) ;
    //q_dot_g = double_vec_dot(&knet->cart[knet->ibz[k]],&R->vec_ai[r]) ;
    shlposi = atoms->shlposn_sh[i1];
    shelposi = atoms->shelposn_sh[i1];
    shelposj = atoms->shelposn_sh[j1];
    //fprintf(file.out,"p %3d r %3d q %3d i1 %3d j1 %3d shlposi %3d shelposi %3d shelposj %3d\n",p,r,q,i1,j1,shlposi, \
    shelposi,shelposj);
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[i1]; index_i++) {
      sheli = shells->type_sh[index_i];
        for (i = 0; i < sheli; i++) {
          for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[j1]; index_j++) {
            shelj = shells->type_sh[index_j];
              for (j = 0; j < shelj; j++) {
                N_ij[offset + count] += Overlap[count] * Density[count] * Complex(cos(q_dot_g), sin(q_dot_g));
                //fprintf(file.out,"%3d %3d %3d %3d %3d %10.4lf %10.4lf\n",i,j,shlposi,shelposi,count,Overlap[count], \
                Density[count]);
                count++;
               }
              }
             }
            shlposi++;
           }
          }
         offset += nbfn * nbfn;
        }
       DestroyDoubleArray(&Overlap,&dim3,job);
       DestroyDoubleArray(&Density,&dim3,job);
      }
     }

}

void fetch_scf_vectors(ComplexMatrix *scf_eigenvectors, FERMI* fermi, ATOM *atoms, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int dima, dimk;
int vector_size;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
char buf2[110], buf3[110];
char xx1[4];
char yy[10] = "/scf_evec";
char zz1[24] = "scf_evectors.";
char zz2[24] = "scf_evec";
MPI_File fh;
FILE *scf_evectors;

  strcpy(buf2,file.directory1);
  strcat(buf2,yy);
  sprintf(xx1, "%d", job->taskid);
  strcat(zz1,xx1);

  dima = atoms->number_of_sh_bfns_in_unit_cell;
  dimk = job->spin_dim * nbands;
  vector_size = dimk * dima;

  if (crystal->type[0] != 'M') {
    scf_evectors = fopen(zz1, "rb");
    fclose(scf_evectors);
   }

    scf_evectors = fopen(zz2, "rb");
    fseek(scf_evectors, dima * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
    fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
    fclose(scf_evectors);

  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, dima * (fermi->bands[0] - 1) * sizeof(Complex), MPI_SEEK_SET) ;
  ////MPI_File_read(fh, &scf_eigenvectors->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  MPI_File_close(&fh);

}

void fetch_scf_eigenvalues(double *scf_eigenvalues, FERMI* fermi, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int nbands = fermi->bands[1] - fermi->bands[0] + 1;
char buf2[110], buf3[110];
char xx1[4];
char yy[10] = "/scf_evec", zz[10] = "/scf_eval";
char zz3[24] = "scf_evalues.";
char zz4[24] = "scf_eval";
double *evalues_temp;
MPI_File gh;
FILE *scf_evalues;

  AllocateDoubleArray(&evalues_temp,&nbands,job);
  ResetDoubleArray(evalues_temp,&nbands);

  if (job->mpi_io == 0) {
  sprintf(xx1, "%d", job->taskid);
  strcat(zz3,xx1);
  if (crystal->type[0] != 'M') {
    scf_evalues = fopen(zz3, "rb");
   }
  if (crystal->type[0] == 'M' && job->taskid == 0) {
    scf_evalues = fopen(zz4, "rb");
    fseek(scf_evalues,(fermi->bands[0] - 1) * sizeof(double),SEEK_SET);
    fread(evalues_temp,sizeof(double),job->spin_dim * nbands,scf_evalues);
    //fread(scf_eigenvalues,sizeof(double),job->spin_dim * nbands,scf_evalues);
    fclose(scf_evalues);
   }  
  MPI_Allreduce(evalues_temp,scf_eigenvalues,job->spin_dim * nbands,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }
  else if (job->mpi_io == 1) {
  strcpy(buf3,file.directory1);
  strcat(buf3,zz);
  MPI_File_open(MPI_COMM_WORLD,buf3,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, (fermi->bands[0] - 1) * sizeof(double), MPI_SEEK_SET) ;
  MPI_File_read(gh, scf_eigenvalues, job->spin_dim * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  MPI_File_close(&gh) ;
 }  
  else {
    if (job->taskid == 0)
    fprintf(file.out,"ERROR: MPI_IO has incorrect value %d in bethe_salpeter\n",job->mpi_io);
    MPI_Finalize(); 
    exit(1);
   }

    DestroyDoubleArray(&evalues_temp,&nbands,job);

}
*/
