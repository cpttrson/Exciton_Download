#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


#include "mycomplex.h"
#include "mylogical.h"
#include "conversion_factors.h"
#include "myconstants.h"

#include "USER_DATA.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "PRINT_UTIL.h"
#include "KPOINTS.h"
#include "CRYSTAL09.h"
#include "PARALLEL.h"
#include "ALLOCATE_MEMORY.h"
#include "ROTATION_OPERATORS.h"
#include "PAIRS_QUADS.h"
#include "INTEGRALS1.h"
#include "INTEGRALS_TWO_CENTRE.h"
#include "INTEGRALS_THREE_CENTRE.h"
#include "DENSITY_MATRIX.h"
#include "DIELECTRIC_FUNCTION.h"
#include "ANALYSIS.h"

using namespace std;

void fermi_surface(int is_fermi[3], int bands[2], CRYSTAL *crystal, SYMMETRY *symmetry, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, m, n, l1, m1, n1, s;
  int rem, num, dim;
  int nbands = bands[1] - bands[0] + 1;
  int k_fbz[(is_fermi[0] + 1) * (is_fermi[1] + 1) * (is_fermi[2] + 1)];
  int nkpts = (is_fermi[0] + 1) * (is_fermi[1] + 1) * (is_fermi[2] + 1);
  double *p_eigval, *eigval;
  char buf3[110];
  char zz1[10] = "/evalfile";
  const char *f;
  KPOINT_TRAN knet;
  FILE *file_fs;
  MPI_File gh ;

  // ******************************************************************************************
  // * Routine calculates eigenvalues on grid for Fermi surface construction by XCrysDen      *
  // ******************************************************************************************

  f = "fermi_surface";
 
    if (job->taskid == 0) {

    if (crystal->type[0] != 'C') {
    fprintf(file.out,"Fermi Surface programmed only for 3-D Crystal\n");
    MPI_Finalize();
    exit(1);
   }

    file_fs = fopen("fermi.bxsf", "w");
    if (file_fs == NULL) { fprintf(file.err, "ERROR: Cannot open file fermi.bxsf\n"); exit(1); }

  // ******************************************************************************************
  // * Generate k points                                                                      *
  // ******************************************************************************************

    count_k_points(&knet,is_fermi,crystal,symmetry,job,file);
    allocate_k_points(&knet,crystal,job,file);
    if (job->C09 == 1)
    read_XCBD_crystal_09(&knet,crystal,job,file);
    generate_k_points(&knet,is_fermi,crystal,symmetry,job,file);

  // ******************************************************************************************
  // * open MPI files and load eigenvalues                                                    *
  // ******************************************************************************************

    dim = job->spin_dim * knet.unique * nbands ;
    AllocateDoubleArray(&eigval,&dim,job);

    strcpy(buf3,file.directory1);
    strcat(buf3,zz1);
    MPI_File_open(MPI_COMM_WORLD,buf3,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;
    MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
    MPI_File_read(gh, eigval, job->spin_dim * knet.unique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"Eigenvalues\n");
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < knet.unique; i++) {
        for (j = 0; j < nbands; j++) {
          fprintf(file.out,"%5d %5d %10.4lf\n",i,j,eigval[s * knet.unique * nbands + i * nbands + j]);
         }
        }
       }
      }

  // ******************************************************************************************
  // * Eigenvalues printed on grid in full Brillouin zone                                     *
  // ******************************************************************************************

      i = 0;
      for (s = 0; s < job->spin_dim; s++) {
        for (l = 0; l < is_fermi[0] + 1; l++) {
          for (m = 0; m < is_fermi[1] + 1; m++) {
            for (n = 0; n < is_fermi[2] + 1; n++) {
              l1 = l;
              m1 = m;
              n1 = n;
              if (l == is_fermi[0]) l1 = 0;
              if (m == is_fermi[1]) m1 = 0;
              if (n == is_fermi[2]) n1 = 0;
              k = l1 * is_fermi[1] * is_fermi[2] + m1 * is_fermi[2] + n1 ;
              k_fbz[i] = knet.fbz[k];
              //fprintf(file.out,"%4d    %4d %4d %4d    %4d %4d\n",i,l1,m1,n1,k,k_fbz[i]);
              i++;
             }
            }
           }
          }

    fprintf(file_fs, " BEGIN_INFO\n   Fermi Energy:%13.6lf\n END_INFO\n BEGIN_BLOCK_BANDGRID_3D\n band_energies\n",
    job->fermi_energy);
    fprintf(file_fs, " BANDGRID_3D_BANDS\n");
    fprintf(file_fs, "%2d\n%2d %2d %2d\n", job->spin_dim * nbands, is_fermi[0] + 1, is_fermi[1] + 1, is_fermi[2] + 1);
    fprintf(file_fs, "      0.00000000      0.00000000      0.00000000\n");
    for (i = 0; i < 3; i++) {
    fprintf(file_fs, "%16.8lf%16.8lf%16.8lf\n", crystal->reciprocal_cell[i].comp1, crystal->reciprocal_cell[i].comp2, \
    crystal->reciprocal_cell[i].comp3);
    }

    for (s = 0; s < job->spin_dim ; s++) {
      for (j = 0; j < nbands ; j++) {
        fprintf(file_fs, "  BAND:%6d\n",j + 1);
        rem = nkpts / 6;
        num = nkpts - rem * 6;
        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_fs, "%14.6e",eigval[s * nbands * knet.unique + nbands * k_fbz[i + 6 * k] + j]) ;
           }
          fprintf(file_fs, "\n");
         }
          for (i = 0; i < num; i++) {
            fprintf(file_fs, "%14.6e",eigval[s * nbands * knet.unique + nbands * k_fbz[i + 6 * k] + j]) ;
           }
          if (num > 0)
          fprintf(file_fs, "\n");
         }
        }
       fprintf(file_fs, " END_BANDGRID_3D\n END_BLOCK_BANDGRID_3D\n");

    DestroyDoubleArray(&eigval,&dim,job);
    free_k_points(&knet,job);

    MPI_File_close(&gh) ;

  } // end if(job->taskid == 0)

}

void state_density(int is[3], int bands[4], int npoints, int nprojections, double atom_proj[][5], double energy_range[2], char spectrum_type, ATOM *atoms,  ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int begin_k1[job->numtasks], end_k1[job->numtasks];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int i, j, k, l, m, p, q, s, v, n;
  int dim, dim1, dimf, dimp, dimf_spin, dimp_spin;
  int dim4 = atoms->number_of_sh_bfns_in_unit_cell;
  int num_sym = symmetry->number_of_operators;
  int nbands = bands[1] - bands[0] + 1;
  int start, finish, vertices, polygons;
  int k_vertex[6][4];
  int Function[8] ;
  int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
  int block_size = vector_size * sizeof(Complex);
  double *eigval, *p_eigenvalues;
  double energy, range = energy_range[1] - energy_range[0];
  double increment = range / (double) (npoints + 1);
  double E, E_tmp, E_vertex[6][4][job->spin_dim * nbands], vol, fac;
  double dos[job->spin_dim * (nprojections + 1)][npoints + 1];
  double dos1[2 * job->spin_dim * (nprojections + 1)][npoints + 1];
  //CHANGES2014double dos1[job->spin_dim * (nprojections + 1)][npoints + 1];
  double projection_weight[nbands][nprojections + 1];
  //double projection_weight[nprojections + 1][nkpoints][nbands];
  double time1, time2;
  double *P, *F;
  char buf2[110], buf3[110];
  char zz1[10] = "/evalfile", yy1[10] = "/datafile";;
  FILE *spect;
  ComplexMatrix *eigvec0, *eigvec1, *eigvec2, *S_k, *S_x, *xtmp ;
  FERMI fermi;
  PAIR_TRAN pair_p;
  KPOINT_TRAN knet;
  MPI_File fh ;
  MPI_File gh ;
  INT_1E one_ints, one_ints_buffer;

      int info = 0;
      char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
      char uplo = 'U';
      char jobz = 'V';
      Complex alpha, beta;
      alpha.real() = k_one;
      alpha.imag() = k_zero;
      beta.real() = k_zero;
      beta.imag() = k_zero;

  // ******************************************************************************************
  // * Routine calculates electronic density of states and atom and shell populations         *
  // ******************************************************************************************

   //CHANGES2014for (i = 0; i < job->spin_dim * (nprojections + 1); i++) {
   for (i = 0; i < job->spin_dim * (nprojections + 1); i++) {
    for (j = 0; j < npoints + 1; j++) {
     dos[i][j] = k_zero;
    }
   }
   for (i = 0; i < 2 * job->spin_dim * (nprojections + 1); i++) {
    for (j = 0; j < npoints + 1; j++) {
     dos1[i][j] = k_zero;
    }
   }

  // ******************************************************************************************
  // * Generate k points                                                                      *
  // ******************************************************************************************

  count_k_points(&knet,is,crystal,symmetry,job,file);
  allocate_k_points(&knet,crystal,job,file);
  if (job->C09 == 1)
  read_XCBD_crystal_09(&knet,crystal,job,file);
  generate_k_points(&knet,is,crystal,symmetry,job,file);
  fermi.is[0] = is[0];
  fermi.is[1] = is[1];
  fermi.is[2] = is[2];
  fermi.bands[0] = bands[0];
  fermi.bands[1] = bands[1];
  if (job->spin_polarisation == 1) {
  bands[2] = bands[0];
  bands[3] = bands[1];
  fermi.bands[2] = bands[2];
  fermi.bands[3] = bands[3];
 }
  fermi.nkunique = knet.unique;
  fermi.nktot = knet.nktot;
  //print_knet(&knet, fermi.is, crystal, job, file);
 
  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->print_pairs == 1)
  print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Allocate Memory for population analysis                                                *
  // ******************************************************************************************

  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file);
  dimp_spin = job->spin_dim * job->dimp;
  dimf_spin = job->spin_dim * job->dimf;
  AllocateDoubleArray(&P,&dimp_spin,job);
  AllocateDoubleArray(&F,&dimf_spin,job);

  // ******************************************************************************************
  // * Allocate Memory for density of states                                                  *
  // ******************************************************************************************

  dim = job->spin_dim * knet.unique * nbands ;
  AllocateDoubleArray(&eigval,&dim,job);
  AllocateComplexMatrix(&eigvec0,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec1,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec2,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&S_k,&atoms->number_of_sh_bfns_in_unit_cell,&atoms->number_of_sh_bfns_in_unit_cell,job);

  int dim1a = atoms->number_of_sh_bfns_in_unit_cell;
  AllocateComplexMatrix(&S_x,&nbands,&nbands,job);
  AllocateComplexMatrix(&xtmp,&nbands,&dim1a,job);


  // ******************************************************************************************
  // * Initialise parameters for parallel run                                                 *
  // ******************************************************************************************

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);
  mpi_begin_end(begin_k1,end_k1,knet.nktot,job->numtasks,job,file);

  // ******************************************************************************************
  // * Generate overlap matrix elements                                                       *
  // ******************************************************************************************

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);

  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  fock_element_1e(&one_ints_buffer, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allgatherv(&one_ints_buffer.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);
  //fock_element_1e(&one_ints, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //MPI_Allgatherv(&one_ints.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  if (job->verbosity > 1) {
  dim = 0 ;
  for (p=0;p<pair_p.nump;p++) {
  fprintf(file.out,"%d %d %d %d\n",p,pair_p.cell1[pair_p.posn[p]],pair_p.cell2[pair_p.posn[p]],pair_p.latt2[pair_p.posn[p]]);
  for(i=0;i<atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]];i++) {
  for(j=0;j<atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];j++) {
  fprintf(file.out,"%5.2lf",one_ints.Overlap[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  fprintf(file.out,"\n") ; }
  fflush(file.out);
 }

  // ******************************************************************************************
  // * open MPI files and load eigenvalues                                                    *
  // ******************************************************************************************

  //strcpy(buf2,file.directory1);
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,yy1);
  //strcpy(buf3,file.directory1);
  strcpy(buf3,file.scf_eigvec);
  strcat(buf3,zz1);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  MPI_File_open(MPI_COMM_WORLD,buf3,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
  MPI_File_read(gh, eigval, job->spin_dim * knet.unique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Eigenvalues\n");
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < knet.unique; i++) {
      for (j = 0; j < nbands; j++) {
        fprintf(file.out,"%5d %5d %10.4lf\n",i,j,*(eigval + s * knet.unique * nbands + i * nbands + j));
       }
      }
     }
    }

  // ******************************************************************************************
  // * Calculate Fermi level                                                                  *
  // ******************************************************************************************
 
  allocate_fermi(&fermi,atoms,job,file);
  calculate_fermi_level_metal_old(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);
  //int working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2);
  //int target_states = working_value * knet.nktot;

  int count = 0;
  int count1 = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < knet.unique; j++) {
      count += job->spin_fac * fermi.occupied[count1] * knet.num[j];
      count1++;
     }
    }
  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Homo for each k point\n");
  count1 = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < knet.unique; j++) {
      fprintf(file.out,"fermi %d %d %d\n",i,j, fermi.occupied[i * knet.unique + j]);
      //for (k = 0; k <= fermi.occupied[i * knet.unique + j]; k++)
        //fprintf(file.out,"%12.6lf\n",fermi.occupation[i * knet.unique * nbands + j * nbands + k]);
       }
      }
     }
               
  // ******************************************************************************************
  // * Calculate atomic populations                                                           *
  // ******************************************************************************************

printf("begin density_matrix\n");
fflush(stdout);

  job->guess_type = 0;
  job->density = 2;
  job->fix_occ = 1;
  fflush(file.out);
  density_matrix_crystal2(&fermi,P,F,&knet,yy1,&fermi.nkunique,R,R_tables,&pair_p,atom_p,atoms,shells,symmetry,job,file);
  //atom_shell_populations2(&one_ints, P, &pair_p, atoms, shells, symmetry, job, file);

printf("finished density_matrix\n");
fflush(stdout);

        if (job->taskid == 0) {
        //fprintf(file.out,"===========================================================================================================\n");
        //fprintf(file.out,"|                          DENSITY OF STATES AND ATOM POPULATION CALCULATION                              |\n");
        //fprintf(file.out,"===========================================================================================================\n");
        //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| E  %9.2e to %9.2e | BAND ENERGY RANGE     %10.3e to %10.3e eV | FERMI ENERGY %10.4lf |\n", \
        energy_range[0] * au_to_eV,energy_range[1] * au_to_eV,*eigval * au_to_eV,*(eigval + nbands - 1) * au_to_eV,job->fermi_energy);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| UNIQUE PAIRS       %6d | TOTAL PAIRS      %6d | INTERACTION RANGE %5.2lf |                         |\n", \
        is[0],is[1],is[2],pair_p.nump,pair_p.tot,R->cutoff * bohr_to_AA);
        print_atom_populations2(atoms, shells, job, file);
        print_shell_populations2(atoms, shells, job, file);
        fflush(file.out);
        fprintf(file.out,"===========================================================================================================\n");
        fflush(file.out);
       }

  // ******************************************************************************************
  // * Tetrahedron method parameters                                                          *
  // ******************************************************************************************

  polygons_vertices_state_density(is,&polygons,&vertices,&vol,crystal,job,file);
  //changes2014polygons_vertices_vol(is,&polygons,&vertices,&vol,crystal,job,file);

  // ******************************************************************************************
  // * Calculate density of states                                                            *
  // ******************************************************************************************

   for (i = 0; i < job->spin_dim; i++) {
    for (j = begin_k1[job->taskid]; j < end_k1[job->taskid]; j++) {

     for (n = 0; n < nbands; n++) {
       projection_weight[n][0] = k_one;
      for (m = 1; m <= nprojections; m++) {
       projection_weight[n][m] = k_zero;
     }
    }

    tetrahedron_vertices2(is, j, k_vertex, &knet, crystal, job, file);
    if (job->verbosity >= 1) {
    fprintf(file.out,"k_vertex %d %d\n",j,k_vertex[0][0]);
    fflush(file.out);
   }

    MPI_File_seek(fh, (i * knet.unique + k_vertex[0][0]) * block_size, MPI_SEEK_SET) ;
    MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],nbands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);
    //for (l=0;l<nbands;l++) {
    //for (m=0;m<atoms->number_of_sh_bfns_in_unit_cell;m++) {
    //fprintf(file.out,"spin %d band %d eigvec %d %10.4lf %10.4lf\n",i,l,m,eigvec1->a[l][m].real(),eigvec1->a[l][m].imag());
    //}}

    int nk[2];
    nk[0] = j;
    nk[1] = j;

    fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

  // ******************************************************************************************
  // * Check that overlap matrix is correctly Fourier transformed                             *
  // ******************************************************************************************

      if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"S_k2 orthogonalised %d\n",nk[0]);
      int ii, i1, j1;
      int info = 0;
      char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
      char uplo = 'U';
      char jobz = 'V';
      Complex alpha, beta;
      alpha.real() = k_one;
      alpha.imag() = k_zero;
      beta.real() = k_zero;
      beta.imag() = k_zero;
      ComplexMatrix *xtrn1, *S_k1, *S_k2, *xtmp1, *eigenvectors;
      double *eigenvalues;
      AllocateComplexMatrix(&S_k1,&dim4,&dim4,job);
      AllocateComplexMatrix(&S_k2,&dim4,&dim4,job);
      fourier_transform(&one_ints.Overlap[0], &S_k1->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      fourier_transform(&one_ints.Overlap[0], &S_k2->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      AllocateComplexMatrix(&eigenvectors,&dim4,&dim4,job);
      AllocateDoubleArray(&eigenvalues,&dim4,job);
      ResetDoubleArray(eigenvalues,&dim4);
      DiagonaliseHermitian(&S_k1, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
      //fprintf(file.out,"Eigenvectors of Fourier transformed overlap matrix at k point %3d\n",k);
      //print_complex_matrix(S_k1, file);
      n = 0;
      for (ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
      if (eigenvalues[ii] < 1.0e-04 && job->taskid == 0)  fprintf(file.out,"small eigenvalue %d %d %e\n",nk[0],ii,eigenvalues[ii]);
      if (eigenvalues[ii] < 1.0e-04)  n++;
     }
      int dim11 = dim4 - n;
      AllocateComplexMatrix(&xtrn1,&dim11,&dim4,job);
      AllocateComplexMatrix(&xtmp1,&dim11,&dim4,job);
      for (i1 = n; i1 < dim4; i1++) {
      //fprintf(file.out,"%3d %10.4lf \n",k,*(eigenvalues + k));
      for (j1 = 0; j1 < dim4; j1++) {
      xtrn1->a[i1 - n][j1] = eigenvectors->a[i1][j1] / sqrt(*(eigenvalues + i1));
      //fprintf(file.out,"%10.4lf ", (eigenvectors->a[k][j]).real());
     }
      //fprintf(file.out,"\n");
     }
      //print_complex_matrix(xtrn1, file);
      ResetComplexMatrix(xtmp1);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &S_k2, &beta, &xtmp1);
      ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &xtrn1, &beta, &S_k2);
      fprintf(file.out,"S_k2 orthogonalised %d\n",nk[0]);
      print_complex_matrix(S_k2, file);
      DestroyComplexMatrix(&S_k1,job);
      DestroyComplexMatrix(&S_k2,job);
      DestroyComplexMatrix(&xtrn1,job);
      DestroyComplexMatrix(&xtmp1,job);
      DestroyComplexMatrix(&eigenvectors,job);
      DestroyDoubleArray(&eigenvalues,&dim4,job);
    }

  // ******************************************************************************************
  // * Calculate projection weights                                                           *
  // ******************************************************************************************


      for (p = 1; p <= nprojections; p++) {
        for (n = 0; n < nbands; n++)     {
          for (l = 0; l < atoms->number_of_sh_bfns_in_unit_cell; l++) {
            eigvec2->a[n][l] = eigvec1->a[n][l] * atom_proj[l][p];
            //fprintf(file.out,"%3d %3d %3d %3d   %10.4lf   %10.4lf  %10.4lf\n",i,p,n,l,eigvec2->a[n][l].real(),eigvec2->a[n][l].imag(),atom_proj[l][p]);
           }
          }
           ResetComplexMatrix(S_x);
           ResetComplexMatrix(xtmp);
           ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigvec1, &S_k, &beta, &xtmp);
           ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp, &eigvec2, &beta, &S_x);
           //fprintf(file.out,"S_x %d\n",k);
           //print_complex_matrix(S_x, file);
           for (n = 0; n < nbands; n++) {
             projection_weight[n][p] += S_x->a[n][n].real();
             //fprintf(file.out,"s %3d k %3d proj %3d band %3d %16.10lf %16.10lf %16.4e\n",\
             i,j,p,n,S_x->a[b][b].real(),S_x->a[b][b].imag(),fabs(k_one - S_x->a[b][b].real()));
            } // close loop over b
           } // close loop over p

        //for (n = 0; n < nbands; n++)     {
          //for (p = 0; p <= nprojections; p++) {
            //fprintf(file.out," %14.8lf",projection_weight[n][p]);
           //}
            //fprintf(file.out,"   sum %14.8lf \n", projection_weight[n][1] + projection_weight[n][2] + projection_weight[n][3]  );
           //}

          dim1 = i * knet.unique * nbands;
           for (p = 0; p < polygons; p++) {
            for (v = 0; v < vertices ; v++) {
             p_eigenvalues = eigval;
             for (n = 0; n < nbands; n++) {
              E_vertex[p][v][n] = *(dim1 + k_vertex[p][v] * nbands + p_eigenvalues) ;
              //fprintf(file.out,"E_vertex %d %d %d %lf\n",p,v,n,dim1 + k_vertex[p][v] * nbands + p_eigenvalues,\
              E_vertex[p][v][n]);
              p_eigenvalues++ ;
              }
             }
            }

           for (p = 0; p < polygons; p++) {
            for (v = 0; v < vertices - 1; v++) {
             for (k = 0; k < vertices - 1; k++) {
              for (n = 0; n < nbands; n++) {
                if (E_vertex[p][k][n] < E_vertex[p][k + 1][n]) {
                  E_tmp = E_vertex[p][k][n];
                  E_vertex[p][k][n] = E_vertex[p][k + 1][n];
                  E_vertex[p][k+1][n] = E_tmp;
                }
               }
              }
             }
            }

   for (p = 0; p < polygons; p++) {
    //changes2014for (v = 0; v < vertices; v++) {
     for (n = 0; n < nbands; n++) {

    start = (int) ((E_vertex[p][vertices - 1][n] - energy_range[0]) / increment);
    finish = (int) ((E_vertex[p][0][n] - energy_range[0]) / increment) + 1;
    if (start < 0)
      start = 0;
    if (start > npoints)
      continue;
    if (finish < 0)
      continue;
    if (finish > npoints)
      finish = npoints;

    for (k = start; k <= finish; k++) {
      E = increment * (double) k + energy_range[0];
      fac = k_zero;

      switch (crystal->type[0]) {

        case 'C':

      if (E > E_vertex[p][3][n] && E < E_vertex[p][2][n]) {
      fac = (E - E_vertex[p][3][n]) * (E - E_vertex[p][3][n]) / (E_vertex[p][2][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
                - E_vertex[p][3][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]);
      //fprintf(file.out,"fac k %d %lf %lf\n",k,E,fac);
      //fflush(file.out);
          }

      if (E > E_vertex[p][2][n] && E < E_vertex[p][1][n]) {
      fac = -(E - E_vertex[p][1][n]) * (E - E_vertex[p][3][n]) / (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][1][n]
                - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]) -\
                 (E - E_vertex[p][0][n]) * (E - E_vertex[p][2][n])
               / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][2][n]) / \
                 (E_vertex[p][0][n] - E_vertex[p][3][n]);
      //fprintf(file.out,"fac k %d %lf %lf\n",k,E,fac);
      //fflush(file.out);
          }

     if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {
     fac = (E - E_vertex[p][0][n]) * (E - E_vertex[p][0][n]) / (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
              - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]);
      //fprintf(file.out,"fac k %d %lf %lf\n",k,E,fac);
      //fflush(file.out);
          }

          fac *= vol;

          break;

        case 'S':

     if (E > E_vertex[p][2][n] && E < E_vertex[p][1][n]) {
     fac = (E - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][2][n]);
          }

     if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {
     fac = -(E - E_vertex[p][0][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]);
          }

          fac *= vol;

          break;

        case 'P':

          fprintf(file.out, "1-D DOS code not complete\n");
          exit(1);

          break;

      } // close switch (crystal->type

         for (l = 0; l <= nprojections; l++) {
          dos[i * (nprojections + 1) + l][k] += fac * projection_weight[n][l];
         }

              } // close loop over k
             } // close loop over n
            //} // close loop over v
           } // close loop over p
          } // close loop over j
         } // close loop on i

printf("finished calc\n");
fflush(stdout);

      for (j = 0; j < job->spin_dim * (nprojections + 1); j++) {
      MPI_Reduce(&dos[j][0],&dos1[j][0],npoints + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     }

printf("finished MPI_Reduce\n");
fflush(stdout);

     if (job->taskid == 0) {
      for (j = 1; j < npoints + 1; j++) {
        energy = energy_range[0] + range * (double) j / double(npoints);
        for (s = 0; s < job->spin_dim; s++) {
          for (l = 0; l < (nprojections + 1); l++) {
            dos1[(job->spin_dim + s) * (nprojections + 1) + l][j] = dos1[(job->spin_dim + s) * (nprojections + 1) + l][j - 1] + \
           (dos1[s * (nprojections + 1) + l][j - 1] + dos1[s * (nprojections + 1) + l][j]) / two * increment * au_to_eV;
           }
          }
         }
      spect = fopen("state_density.dat", "w");
       if (spect == NULL) { fprintf(file.out, "cannot open file state_density.dat\n"); exit(1); }
        for (j = 0; j < npoints; j++) {
          energy = energy_range[0] + range * (double) j / double(npoints);
          fprintf(spect, "%10.4e   ", energy * au_to_eV);
          //CHANGES2014for (l = 0; l < job->spin_dim * (nprojections + 1); l++) {
          fac = k_one;
          for (s = 0; s < job->spin_dim; s++) {
            for (l = 0; l < 2; l++) {
              for (m = 0; m < nprojections + 1; m++) {
            //CHANGES2014for (l = 0; l < (nprojections + 1); l++) {
              //CHANGES2014fprintf(spect, "%12.4e", fac * dos1[s * (nprojections + 1) + l][j]);
              fprintf(spect, "%12.4e", fac * dos1[s * 2 * (nprojections + 1) + l * (nprojections + 1) + m][j]);
              //fprintf(spect, "%12.4e", dos[l][j]);
             }
            fac *= -k_one;
            }
           }
          fprintf(spect,"\n");
         }
        fflush(spect);;
        fclose(spect);
       } // end if job->taskid

  DestroyDoubleArray(&eigval,&dim,job);
  DestroyDoubleArray(&P,&dimp_spin,job);
  DestroyDoubleArray(&F,&dimf_spin,job);
  DestroyComplexMatrix(&S_k,job);
  DestroyComplexMatrix(&S_x,job);
  DestroyComplexMatrix(&xtmp,job);
  DestroyComplexMatrix(&eigvec0,job);
  DestroyComplexMatrix(&eigvec1,job);
  DestroyComplexMatrix(&eigvec2,job);
  free_INT_1E(&one_ints, Function, job, file);
  free_k_points(&knet,job);
  free_PAIR_TRAN(&pair_p,job);

  MPI_File_close(&fh) ;
  MPI_File_close(&gh) ;

}

void susceptibility(int is[3], int bands[2], int npoints, double energy_range[2], VECTOR_DOUBLE *field, ATOM *atoms,  ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int begin_k1[job->numtasks], end_k1[job->numtasks];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int receive_s[job->numtasks], offset_s[job->numtasks];
  int nbands = bands[1] - bands[0] + 1;
  int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
  int block_size = vector_size * sizeof(Complex);
  int i, j, k, l, m, p, q, s, v, n, nq;
  int i1, j1, j3;
  int dim, dim1, dim3, dimf;
  int dim5;
  int start, finish, vertices, polygons;
  int k_vertex[6][4], kq_vertex[6][4];
  int nk[2];
  int Function[8] ;
  int count, count1;
  double energy, range = energy_range[1] - energy_range[0];
  double increment = range / (double) (npoints + 1);
  double E, E_tmp, vol, fac,  fac_real, fac_imag;
  double E_i, E_j;
  double spectrum[job->spin_dim][2 * job->field_dirs][npoints + 1];
  double spectrum1[job->spin_dim][3 * job->field_dirs][npoints + 1];
  double tempr, tempi;
  double sum_rule_fac;
  double time1, time2;
  double *eigval, *eigvec;
  char buf2[110], buf3[110];
  char zz1[10] = "/evalfile", yy1[10] = "/datafile";
  char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
  FILE *spect;
  ComplexMatrix *eigvec0, *eigvec1, *eigvec2, *M_k, *tmp, *M_x, *S_k ;
  VECTOR_INT k_point, q_point, kq_point;
  PAIR_TRAN pair_p;
  KPOINT_TRAN knet;
  INT_1E one_ints, one_ints_buffer;
  FERMI fermi;
  MPI_File fh ;
  MPI_File gh ;
  Complex alpha, beta;
  alpha.real() = k_one;
  alpha.imag() = k_zero;
  beta.real() = k_zero;
  beta.imag() = k_zero;

  // ******************************************************************************************
  // * Routine calculates long wavelength dielectric susceptibility as function of frequency  *
  // ******************************************************************************************

   for (i = 0; i < job->spin_dim; i++) {
    for (j3 = 0; j3 < 2 * job->field_dirs; j3++) { // 3 is for real, imaginary parts, sum rule
     for (j = 0; j < npoints + 1; j++) {
     spectrum[i][j3][j] = k_zero;
     spectrum1[i][j3][j] = k_zero;
    }
   }
  }

  // ******************************************************************************************
  // * Generate k points                                                                      *
  // ******************************************************************************************

  count_k_points(&knet,is,crystal,symmetry,job,file);
  allocate_k_points(&knet,crystal,job,file);
  if (job->C09 == 1)
  read_XCBD_crystal_09(&knet,crystal,job,file);
  generate_k_points(&knet,is,crystal,symmetry,job,file);
  fermi.nkunique = knet.unique;

  // ******************************************************************************************
  // * Allocate Memory                                                                        *
  // ******************************************************************************************

  dim = job->spin_dim * knet.unique * nbands ;
  dim3 = atoms->number_of_sh_bfns_in_unit_cell * 3;
  dim5 = nbands * 3;

  AllocateDoubleArray(&eigval,&dim,job);
  AllocateComplexMatrix(&eigvec0,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec1,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec2,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&S_k,&atoms->number_of_sh_bfns_in_unit_cell,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&M_k,&dim3,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&tmp,&atoms->number_of_sh_bfns_in_unit_cell,&nbands,job);
  AllocateComplexMatrix(&M_x,&dim5,&nbands,job);
  ResetComplexMatrix(M_x);
  ResetComplexMatrix(tmp);

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->print_pairs == 1)
  print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Initialise parameters for parallel run                                                 *
  // ******************************************************************************************

  //mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  //mpi_receive_offset_momentum(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);
  mpi_begin_end(begin_k1,end_k1,knet.nktot,job->numtasks,job,file);

  // ******************************************************************************************
  // * open MPI files and load eigenvalues                                                    *
  // ******************************************************************************************

  strcpy(buf2,file.directory1);
  strcat(buf2,yy1);
  strcpy(buf3,file.directory1);
  strcat(buf3,zz1);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  MPI_File_open(MPI_COMM_WORLD,buf3,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
  MPI_File_read(gh, eigval, job->spin_dim * knet.unique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Eigenvalues\n");
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < knet.unique; i++) {
      for (j = 0; j < nbands; j++) {
        fprintf(file.out,"%5d %5d %10.4lf\n",i,j,*(eigval + s * knet.unique * nbands + i * nbands + j)*au_to_eV);
       }
        fprintf(file.out,"\n");
      }
     }
    }

  // ******************************************************************************************
  // * Generate momentum matrix elements                                                      *
  // ******************************************************************************************

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 1 ;
  Function[4] = 0 ;
  //Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_momentum(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);
  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  fock_element_1e(&one_ints_buffer, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allgatherv(&one_ints_buffer.Momentum[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Momentum[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);

/*
  INT_1E one_ints_buffer;
  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  time1 = MPI_Wtime();
  fock_element_1e1(&one_ints_buffer, dim, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  time2 = MPI_Wtime() - time1;
  printf("Time to generate 1e matrix %10.4e\n",(double)(time2));
  MPI_Allreduce(&one_ints_buffer.Momentum[0],&one_ints.Momentum[0],3*dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);
*/

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Momentum Matrix Elements in Real Space\n");
  dim = 0 ;
  for (p = 0; p < pair_p.nump; p++) {
  fprintf(file.out,"Pair[%4d]  %4d %4d   %4d\n",p,pair_p.cell1[pair_p.posn[p]],pair_p.cell2[pair_p.posn[p]],pair_p.latt2[pair_p.posn[p]]);
  for(i = 0; i < atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]]; i++) {
  for(j = 0; j < atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]]; j++) {
  fprintf(file.out,"%5.2lf",one_ints.Momentum[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  fprintf(file.out,"\n") ;
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  for(i = 0; i < atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]]; i++) {
  for(j = 0; j < atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]]; j++) {
  fprintf(file.out,"%5.2lf",one_ints.Momentum[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  fprintf(file.out,"\n") ;
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  for(i = 0; i < atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]]; i++) {
  for(j = 0; j < atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]]; j++) {
  fprintf(file.out,"%5.2lf",one_ints.Momentum[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  fprintf(file.out,"\n") ;
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  fprintf(file.out,"\n") ; }
 }

  // ******************************************************************************************
  // * Calculate Fermi level                                                                  *
  // ******************************************************************************************

  allocate_fermi(&fermi,atoms,job,file);
  calculate_fermi_level_metal_old(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);

  ////calculate_fermi_level_metal_old(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);
  //FIX calculate_fermi_level_metal(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);
  ////allocate_fermi(&fermi,atoms,job,file);
  ////occupation_fermi(&fermi, &knet, eigval, knet.unique, nbands, job, file);
  //CHANGES2014occupation_fermi(&fermi, eigval, knet.unique, nbands, job, file);
  ////int working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2);
  ////int target_states = working_value * knet.nktot;

  count = 0;
  count1 = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < knet.unique; j++) {
      count += fermi.occupied[count1] * knet.num[j];
      count1++;
      //count += fermi.occupied[i] * knet.num[i] * job->spin_fac;
     }
    }
  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Homo for each k point\n");
  for (i = 0; i < job->spin_dim * knet.unique; i++)
  fprintf(file.out,"fermi %d %d\n",i,fermi.occupied[i]);
 }
  //fprintf(file.out,"occupied states %d target states %d total electrons %d\n",count, target_states, count/knet.nktot + 2 * bands[0] - 2);

  // ******************************************************************************************
  // * Tetrahedron method parameters                                                          *
  // ******************************************************************************************

  polygons_vertices_susceptibility(is,&polygons,&vertices,&vol,crystal,job,file);
  //changes2014polygons_vertices_vol(is,&polygons,&vertices,&vol,crystal,job,file);

  // ******************************************************************************************
  // * Calculate dielectric susceptibilty                                                     *
  // ******************************************************************************************

        if (job->taskid == 0) {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                     DIELECTRIC FUNCTION CALCULATION                                     |\n");
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"| E  %9.3e to %9.3e | BAND ENERGY RANGE      %9.3e to %9.3e eV | FERMI ENERGY   %8.5lf |\n", \
        energy_range[0] * au_to_eV,energy_range[1] * au_to_eV,*eigval * au_to_eV,*(eigval + nbands - 1) * au_to_eV,job->fermi_energy);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| BAND RANGE   %4d to %4d | OCCUPIED STATES         | TARGET STATES           | TOTAL ELECTRONS         |\n", \
        bands[0],bands[1]);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| MONKHORST-PACK   %2d %2d %2d | UNIQUE PAIRS     %6d | TOTAL PAIRS      %6d | INTERACTION RANGE %5.2lf |\n", \
        is[0],is[1],is[2],pair_p.nump,pair_p.tot,R->cutoff * bohr_to_AA);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                                     APPLIED FIELD DIRECTION COSINES                                     |\n");
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        for (i = 0; i < job->field_dirs; i++) {
        fprintf(file.out,"| %3d                       |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        i + 1,field[i].comp1,field[i].comp2,field[i].comp3);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
       }
      fflush(file.out);
     }

   nq = 1;

  for (q = 0; q < nq; q++) {
    for (j = begin_k1[job->taskid]; j < end_k1[job->taskid]; j++) {

    //changes2014a int ntransitions = (nbands - fermi.occupied[knet.fbz[j]]) * fermi.occupied[knet.fbz[j]];
    //if (ntransitions <= 0) {
    //fprintf(file.out,"Error in susceptibility. Fermi energy not in range of bands. ntransitions = %d\n",ntransitions);
   //}

    //double E_vertex[polygons][vertices][job->spin_dim * ntransitions];
    //double matrix_element_fac[job->field_dirs][ntransitions];
    //changes2014a double matrix_element[job->field_dirs][ntransitions];

    k_point = decompose_k_point(is,j,crystal,job,file);
    q_point = decompose_k_point(is,0,crystal,job,file);
    p = compose_k_point(is,k_point,q_point,crystal,job,file);

    tetrahedron_vertices2(is, j,  k_vertex, &knet, crystal, job, file);
    tetrahedron_vertices2(is, j, kq_vertex, &knet, crystal, job, file);

    if (job->taskid == 0 && job->verbosity >= 1) {
    fprintf(file.out,"j %d k_point %3d %3d %3d k_vertex %3d kq_vertex %d\n",j, \
    k_point.comp1,k_point.comp2,k_point.comp3,k_vertex[0][0], kq_vertex[0][0]);
    //fprintf(file.out,"j %d k_point %3d %3d %3d k_vertex %3d kq_vertex %d\n",j,k_point.comp1,k_point.comp2,k_point.comp3,k_vertex[0][0], kq_vertex[0][0]);
    fflush(file.out);
   }

    nk[0] = j;
    nk[1] = j;

  // ******************************************************************************************
  // * Fourier transform momentum matrix elements                                             *
  // ******************************************************************************************

    fourier_transform_3(&one_ints.Momentum[0], &M_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

    //fprintf(file.out,"M_k\n");
    //print_complex_matrix(M_k,file);

    for (i = 0; i < job->spin_dim; i++) {

    int ntransitions = (nbands - fermi.occupied[i * knet.unique + knet.fbz[j]]) * fermi.occupied[i * knet.unique + knet.fbz[j]];
    if (ntransitions <= 0) {
    fprintf(file.out,"Error in susceptibility. Fermi energy not in range of bands. ntransitions = %d\n",ntransitions);
   }

    double E_vertex[polygons][vertices][job->spin_dim * ntransitions];
    double matrix_element_fac[job->field_dirs][ntransitions];
    double matrix_element[job->field_dirs][ntransitions];

  // ******************************************************************************************
  // * Read and rotate wave functions                                                         *
  // ******************************************************************************************

        //bands[0] = fermi->bands[2 * i] - 1;
        //bands[1] = fermi->bands[2 * i + 1] - 1;
        //CHANGE2014 fix for spin pol case
        //nbands = bands[1] - bands[0] + 1;

    MPI_File_seek(fh, (i * knet.unique +  k_vertex[0][0]) * block_size, MPI_SEEK_SET) ;
    MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],nbands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);
    MPI_File_seek(fh, (i * knet.unique + kq_vertex[0][0]) * block_size, MPI_SEEK_SET) ;
    MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    rotate_psi(&eigvec0->a[0][0],&eigvec2->a[0][0],nbands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);

    for (int i2 = 0; i2 < nbands; i2++) {
      for (int j2 = 0; j2 < nbfn; j2++) {
        eigvec1->a[i2][j2].imag() *= -k_one;
       }
      }

  // ******************************************************************************************
  // * Contract momentum matrix elements with wave functions                                  *
  // ******************************************************************************************

    for (j3 = 0; j3 < 3; j3++) {
    dim = j3 * atoms->number_of_sh_bfns_in_unit_cell;
    dim5 = j3 * nbands;
    ComplexGEMM3(&NoTrans,&Trans,&nbfn,&nbands,&nbfn,&alpha,&(M_k->a[dim]),&nbfn,&(eigvec2->a[0]),&nbfn,&beta,&(tmp->a[0]),&nbands);
    ComplexGEMM3(&NoTrans,&NoTrans,&nbands,&nbands,&nbfn,&alpha,&(eigvec1->a[0]),&nbfn,&(tmp->a[0]),&nbands,&beta,&(M_x->a[dim5]),&nbands);
   }

    if (job->verbosity > 1 && j < 6) {
    fprintf(file.out,"M_k\n");
    print_complex_matrix(M_k,file);
    fprintf(file.out,"\neigvec1\n");
    print_complex_matrix(eigvec1,file);
    fprintf(file.out,"\neigvec2\n");
    print_complex_matrix(eigvec2,file);
    fprintf(file.out,"tmp\n");
    print_complex_matrix(tmp,file);
    fprintf(file.out,"\nM_x\n");
    print_complex_matrix(M_x,file);
    for (int ii = 0; ii < nbands; ii++) {
      for (int jj = 0; jj < atoms->number_of_sh_bfns_in_unit_cell; jj++) {
        fprintf(file.out,"k %2d %2d %16.10lf %16.10lf   %16.10lf %16.10lf\n",\
        k_vertex[0][0],j, \
        (eigvec0->a[ii][jj]).real(),(eigvec0->a[ii][jj]).imag(), \
        (eigvec2->a[ii][jj]).real(),(eigvec2->a[ii][jj]).imag());
       }
      fprintf(file.out,"\n");
     }
   }

    fprintf(file.out,"M_k %3d\n",j);
    print_complex_matrix(M_k,file);
    fprintf(file.out,"\nM_x\n");
    print_complex_matrix(M_x,file);

    for (j3 = 0; j3 < job->field_dirs; j3++) {
      count = 0;
      matrix_element[j3][0] = k_zero;
      for (i1 = 0; i1 < fermi.occupied[i * knet.unique + knet.fbz[j]]; i1++) {
        for (j1 = fermi.occupied[i * knet.unique + knet.fbz[j]]; j1 < nbands; j1++) {
          matrix_element_fac[j3][count] = (field[j3].comp1 * M_x->a[i1][j1].real() + \
                                           field[j3].comp2 * M_x->a[nbands + i1][j1].real() + \
                                           field[j3].comp3 * M_x->a[2 * nbands + i1][j1].real()) * \
                                          (field[j3].comp1 * M_x->a[i1][j1].real() + \
                                           field[j3].comp2 * M_x->a[nbands + i1][j1].real() + \
                                           field[j3].comp3 * M_x->a[2 * nbands + i1][j1].real()) + \
                                          (field[j3].comp1 * M_x->a[i1][j1].imag() + \
                                           field[j3].comp2 * M_x->a[nbands + i1][j1].imag() + \
                                           field[j3].comp3 * M_x->a[2 * nbands + i1][j1].imag()) * \
                                          (field[j3].comp1 * M_x->a[i1][j1].imag() + \
                                           field[j3].comp2 * M_x->a[nbands + i1][j1].imag() + \
                                           field[j3].comp3 * M_x->a[2 * nbands + i1][j1].imag());
           count++;
          }
         }
        }

//fprintf(file.out,"mat ele %3d %10.4lf\n",j,matrix_element_fac[0][7]);

    dim1 = i * knet.unique * nbands;
    for (p = 0; p < polygons; p++) {
      for (v = 0; v < vertices ; v++) {
        count = 0;
        //changes2014a for (m = 0; m < fermi.occupied[knet.fbz[j]]; m++) {
          //for (n = fermi.occupied[knet.fbz[j]]; n < nbands; n++) {
        for (m = 0; m < fermi.occupied[i * knet.unique + knet.fbz[j]]; m++) {
          for (n = fermi.occupied[i * knet.unique + knet.fbz[j]]; n < nbands; n++) {
            //E_vertex[p][v][count] = *(eigval + dim1 + kq_vertex[p][v] * nbands + n) - *(eigval + dim1 + k_vertex[p][v] * nbands + m) + 0.0001 * v;
            E_vertex[p][v][count] = *(eigval + dim1 + kq_vertex[p][v] * nbands + n) - *(eigval + dim1 + k_vertex[p][v] * nbands + m);

//if (p == 0 && v == 0 && E_vertex[p][v][count] > 4.4/27.211 && E_vertex[p][v][count] < 4.6/27.211 && (matrix_element_fac[0][count] > 1e-02 || matrix_element_fac[1][count] > 1e-02 || matrix_element_fac[2][count] > 1e-00))  fprintf(file.out,"point1 %3d %3d %3d %3d %e %e %e %e\n",j,i,m,n,E_vertex[p][v][count],matrix_element_fac[0][count],matrix_element_fac[1][count],matrix_element_fac[2][count]);
            count++;
           }
          }
         }
        }

    for (p = 0; p < polygons; p++) {
      for (v = 0; v < vertices - 1; v++) {
        for (k = 0; k < vertices - 1; k++) {
          for (n = 0; n < ntransitions; n++) {
            if (E_vertex[p][k][n] < E_vertex[p][k + 1][n]) {
              E_tmp = E_vertex[p][k][n];
              E_vertex[p][k][n] = E_vertex[p][k + 1][n];
              E_vertex[p][k+1][n] = E_tmp;
             }
            }
           }
          }
         }

    for (p = 0; p < polygons; p++) {
      for (v = 0; v < vertices; v++) {
        for (n = 0; n < ntransitions; n++) {
          start = (int) ((E_vertex[p][vertices - 1][n] - energy_range[0]) / increment);
          finish = (int) ((E_vertex[p][0][n] - energy_range[0]) / increment) + 1;
          //fprintf(file.out,"start %d finish %d p v n %d %d %d E_vertex %lf %lf %lf %lf\n",start,finish,p,v,n,E_vertex[p][0][n],E_vertex[p][1][n],\
          E_vertex[p][2][n], increment * (double) finish + energy_range[0]);
          //start = 0; finish = npoints;
          if (start < 0)
          start = 0;
          if (start > npoints)
          continue;
          if (finish < 0)
          continue;
          if (finish > npoints)
          finish = npoints;
          for (k = start; k <= finish; k++) {
            E = increment * (double) k + energy_range[0];
            fac = k_zero;
            fac_real = k_zero;
            fac_imag = k_zero;

      switch (crystal->type[0]) {

        case 'C':

    if (E > E_vertex[p][3][n] && E < E_vertex[p][2][n]) {
 fac = (E - E_vertex[p][3][n]) * (E - E_vertex[p][3][n]) / (E_vertex[p][2][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
          - E_vertex[p][3][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]);
          }

    if (E > E_vertex[p][2][n] && E < E_vertex[p][1][n]) {
fac = -(E - E_vertex[p][1][n]) * (E - E_vertex[p][3][n]) / (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][1][n]
          - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]) -\
       (E - E_vertex[p][0][n]) * (E - E_vertex[p][2][n])
         / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][2][n]) / \
           (E_vertex[p][0][n] - E_vertex[p][3][n]);
          }

    if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {
 fac = (E - E_vertex[p][0][n]) * (E - E_vertex[p][0][n]) / (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
          - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]);
          }

          break;

        case 'S':

    if (E > E_vertex[p][2][n] && E < E_vertex[p][1][n]) {
 fac = (E - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][2][n]);
          }

    if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {
fac = -(E - E_vertex[p][0][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]);
          }

          break;

        case 'P':

          fprintf(file.out, "1-D DOS code not complete\n");
          exit(1);

          break;

      } // close switch (crystal->type

         for (j3 = 0; j3 < job->field_dirs; j3++) {
           //fprintf(file.out,"matrix elem %d %d %d %d %e %e  %e %e %e\n",p,n,k,j3,fac,matrix_element_fac[j3][n], \
           E_vertex[p][0][n],E_vertex[p][1][n],E_vertex[p][2][n]);
           spectrum[i][2 * j3][k] += fac * vol * matrix_element_fac[j3][n] / E / E;
          }

              } // close loop over k
             } // close loop over n
            }  // close loop over v
           } // close loop over p
          } // close loop on i

          } // close loop over j
         } // close loop on q

    for (j = 0; j < job->spin_dim; j++) {
      for (j3 = 0; j3 < 2 * job->field_dirs; j3++) { // factor of 2 is for real, imaginary parts
        MPI_Reduce(&spectrum[j][j3][0],&spectrum1[j][j3][0],npoints + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
       }
      }

    for (j3 = 0; j3 < job->field_dirs; j3++) {
      for (i = 0; i < npoints; i++) {
        E_i = increment * (double) i + energy_range[0];
        if (E_i < 0.0001) continue;
        if ((i/2)*2 == i) {
        for (j = 1; j < npoints; j+=2) {
          E_j = increment * (double) j + energy_range[0];
          spectrum1[0][2 * j3 + 1][i] += two * increment / pi * spectrum1[0][2 * j3][j] * (k_one / (E_j - E_i) + k_one / (E_j + E_i));
         }
        }
      else {
        for (j = 0; j < npoints; j+=2) {
          E_j = increment * (double) j + energy_range[0];
          spectrum1[0][2 * j3 + 1][i] += two * increment / pi * spectrum1[0][2 * j3][j] * (k_one / (E_j - E_i) + k_one / (E_j + E_i));
         }
       }
      }
     }

  // ******************************************************************************************
  // * Trapezoid rule for spectrum sum rule                                                   *
  // ******************************************************************************************
    
    if (job->taskid == 0) {
      sum_rule_fac = two * epsilon_0 * m0 * crystal->primitive_cell_volume * a0 * a0 * a0 * au_to_eV * au_to_eV / hbar / hbar / pi;
      for (j = 0; j < npoints; j++) {
        energy = energy_range[0] + range * (double) j / double(npoints);
        for (l = 0; l < job->spin_dim; l++) {
          for (j3 = 0; j3 < job->field_dirs; j3++) {
            spectrum1[l][j3 + 2 * job->field_dirs][j] = spectrum1[l][j3 + 2 * job->field_dirs][j - 1] + \
            energy * spectrum1[l][2 * j3][j] * increment * sum_rule_fac;
           }
          }
         }
        }

    if (job->taskid == 0) {
      spect = fopen("optical_spectrum.dat", "w");
      if (spect == NULL) { fprintf(file.out, "cannot open file optical_spectrum.dat\n"); exit(1); }
      for (j = 0; j < npoints; j++) {
        energy = energy_range[0] + range * (double) j / double(npoints);
        fprintf(spect, "%10.4e   ", energy * au_to_eV);
        for (l = 0; l < job->spin_dim; l++) {
          for (j3 = 0; j3 < 3 * job->field_dirs; j3++) {
            fprintf(spect, "%12.4e", spectrum1[l][j3][j]);
           }
          }
         fprintf(spect,"\n");
        }
       fflush(spect);;
       fclose(spect);
      } // end if job->taskid

  DestroyComplexMatrix(&S_k,job);
  DestroyComplexMatrix(&M_k,job);
  DestroyComplexMatrix(&M_x,job);
  DestroyComplexMatrix(&tmp,job);
  DestroyComplexMatrix(&eigvec0,job);
  DestroyComplexMatrix(&eigvec1,job);
  DestroyComplexMatrix(&eigvec2,job);
  DestroyDoubleArray(&eigval,&dim,job);
  free_INT_1E(&one_ints, Function, job, file);
  free_k_points(&knet,job);
  free_PAIR_TRAN(&pair_p,job);
  free_fermi(&fermi,job);

  MPI_File_close(&fh) ;
  MPI_File_close(&gh) ;

}

void susceptibility1(int is[3], int bands[2], int npoints, double energy_range[2], VECTOR_DOUBLE *field, ATOM *atoms,  ATOM_TRAN *atom_p, ATOM_TRAN *atom_i, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int begin_k1[job->numtasks], end_k1[job->numtasks];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int num_sym = symmetry->number_of_operators;
  int nbands = bands[1] - bands[0] + 1;
  int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
  int block_size = vector_size * sizeof(Complex);
  int i, j, k, l, m, p, q, v, n, nq;
  int j3;
  int dim, dim1, dim3, dimf;
  int dim5;
  int start, finish, vertices, polygons;
  int k_vertex[6][4], kq_vertex[6][4];
  int nktot, nkunique, ksize;
  int nk[2];
  int Function[8] ;
  int num_pairs[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int pair_limits[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int num_p, max_pairs;
  int count;
  double energy, range = energy_range[1] - energy_range[0];
  double increment = range / (double) (npoints + 1);
  double E, E_tmp, vol, fac,  fac_real, fac_imag;
  double E_i, E_j;
  double spectrum[job->spin_dim][2 * job->field_dirs][npoints + 1];
  double spectrum1[job->spin_dim][3 * job->field_dirs][npoints + 1];
  //double spectrum1[job->spin_dim][2 * job->field_dirs][npoints + 1];
  double tempr, tempi;
  double sum_rule_fac;
  double time1, time2;
  double *eigval, *eigvec;
  char buf2[110], buf3[110];
  char zz1[10] = "/evalfile", yy1[10] = "/datafile";
  char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
  FILE *spect;
  ComplexMatrix *eigvec0, *eigvec1, *eigvec2, *M_k, *tmp, *M_x ;
  VECTOR_INT k_point, q_point, kq_point;
  PAIR_TRAN pair_p;
  KPOINT_TRAN knet;
  INT_1E one_ints;
  FERMI fermi;
  MPI_File fh ;
  MPI_File gh ;
  Complex alpha, beta;
  alpha.real() = k_one;
  alpha.imag() = k_zero;
  beta.real() = k_zero;
  beta.imag() = k_zero;

        //FILE *mat, *mat1, *mat2;
        //mat = fopen("data", "w");
        //mat1 = fopen("data1","w");
        //mat2 = fopen("data2","w");

   for (i = 0; i < job->spin_dim; i++) {
    for (j3 = 0; j3 < 2 * job->field_dirs; j3++) { // 3 is for real, imaginary parts, sum rule
     for (j = 0; j < npoints + 1; j++) {
     spectrum[i][j3][j] = k_zero;
     spectrum1[i][j3][j] = k_zero;
    }
   }
  }

  // Generate set of k points for density of states

  count_k_points(&knet,is,crystal,symmetry,job,file);
  allocate_k_points(&knet,crystal,job,file);
  if (job->C09 == 1)
  read_XCBD_crystal_09(&knet,crystal,job,file);
  generate_k_points(&knet,is,crystal,symmetry,job,file);
  fermi.nkunique = knet.unique;
nkunique = knet.unique;
nktot = fermi.nktot;

/*
  time1 = MPI_Wtime();
  knet_size(&nktot,is,crystal);
  count_k_points(num_sym,is,&nkunique,crystal,symmetry,job,file);
  allocate_k_points(&knet,nkunique,nktot,crystal,job,file);
  if (job->C09 == 1)
  read_XCBD_crystal_09(nkunique,&knet,crystal,job,file);
  generate_k_points(&knet,num_sym,is,nkunique,crystal,symmetry,job,file);
  time2 = MPI_Wtime();
*/
  //if (job->verbosity >= 1)
  //fprintf(file.out, "Time to generate (total %d unique %d) kpoints %10.4e\n",nktot,nkunique,time2 - time1);
  //VECTOR_KNET *knet_unique, *p_knet_unique, *p_knet;
  //KPT_TRAN knet_rotated[(is[0] + 3) * (is[1] + 3) * (is[2] + 3)];
  //knet_size(&ksize,is,crystal);
  //knet_unique = (VECTOR_KNET *) malloc(ksize * sizeof(VECTOR_KNET));
  //if (knet_unique == NULL) { fprintf(stderr, "ERROR: not enough memory for VECTOR_KNET knet_unique! \n"); exit(1); }
  //unique_kpoints1(&knet,knet_unique, &nkunique, &nktot, knet_rotated, num_sym, is, crystal, symmetry, job, file);

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->print_pairs == 1)
  print_pairs(&pair_p,atoms,R,job,file);

/*
  time1 = MPI_Wtime();
  count_pairs(atoms,atom_p,atom_i,symmetry,R,R_tables,&num_p,num_pairs,&max_pairs,pair_limits,job,file);
  allocate_density_pairs(&pair_p,num_p,max_pairs,symmetry,atoms,R_tables,file);
  generate_pairs(atoms,atom_p,atom_i,symmetry,R,R_tables,&pair_p,&num_p,&max_pairs,pair_limits,job,file);
  time2 = MPI_Wtime();
  if (job->verbosity > 1) {
  fprintf(file.out, "Time to generate %d pairs %10.4e\n",max_pairs,time2 - time1);
  print_pairs(&num_p, &pair_p, job, file);
 }
*/

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 1 ;
  Function[4] = 0 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_momentum(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);

  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  //array_dimensions(&dim, &dimf, &num_p, &pair_p, atoms, job, file); // don't change
  fflush(file.out);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e(&one_ints, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //fock_element_1e(&one_ints, dim, &pair_p, num_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(file.out);

  MPI_Allgatherv(&one_ints.Momentum[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Momentum[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);

  if (job->verbosity >= 1) {
  fprintf(file.out,"Momentum\n");
  dim = 0 ;
  for (p=0;p<pair_p.nump;p++) {
  fprintf(file.out,"%d %d %d %d\n",p,pair_p.cell1[pair_p.posn[p]],pair_p.cell2[pair_p.posn[p]],pair_p.latt2[pair_p.posn[p]]);
  for(i=0;i<atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]];i++) {
  for(j=0;j<atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];j++) {
  fprintf(file.out,"%5.2lf",one_ints.Momentum[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  fprintf(file.out,"\n") ;
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  for(i=0;i<atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]];i++) {
  for(j=0;j<atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];j++) {
  fprintf(file.out,"%5.2lf",one_ints.Momentum[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  fprintf(file.out,"\n") ;
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  for(i=0;i<atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]];i++) {
  for(j=0;j<atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];j++) {
  fprintf(file.out,"%5.2lf",one_ints.Momentum[dim + i * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]] + j]) ; }
  fprintf(file.out,"\n") ; }
  fprintf(file.out,"\n") ;
  dim += atoms->bfnnumb_sh[pair_p.cell1[pair_p.posn[p]]] * atoms->bfnnumb_sh[pair_p.cell2[pair_p.posn[p]]];
  fprintf(file.out,"\n") ; }
  fflush(file.out);
 }

  mpi_begin_end(begin_k1,end_k1,knet.nktot,job->numtasks,job,file);
  //mpi_begin_end(begin_k1,end_k1,nktot,job->numtasks,job,file);

  dim = job->spin_dim * knet.unique * nbands ;
  //dim = job->spin_dim * nkunique * nbands ;
  //dim = job->spin_dim * nkunique * nbands ;
  dim3 = atoms->number_of_sh_bfns_in_unit_cell * 3;
  dim5 = nbands * 3;

  AllocateDoubleArray(&eigval,&dim,job);
  AllocateComplexMatrix(&eigvec0,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec1,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec2,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&M_k,&dim3,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&tmp,&atoms->number_of_sh_bfns_in_unit_cell,&nbands,job);
  AllocateComplexMatrix(&M_x,&dim5,&nbands,job);
  ResetComplexMatrix(M_x);
  ResetComplexMatrix(tmp);

  polygons_vertices_susceptibility(is,&polygons,&vertices,&vol,crystal,job,file);
  //changes2014polygons_vertices_vol(is,&polygons,&vertices,&vol,crystal,job,file);

  // open MPI file, load eigenvectors and eigenvalues

  strcpy(buf2,file.directory1);
  strcat(buf2,yy1);
  strcpy(buf3,file.directory1);
  strcat(buf3,zz1);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  MPI_File_open(MPI_COMM_WORLD,buf3,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;

  MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
  MPI_File_read(gh, eigval, job->spin_dim * knet.unique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;

  if (job->verbosity >= 1) {
  for(i=0;i<knet.unique;i++) {
  for(j=0;j<nbands;j++) {
  fprintf(file.out,"eigval1 %d %d %lf %lf\n",i,j,*(eigval + i * nbands + j),*(eigval + 0*nkunique * nbands + i * nbands + j)); }}
  fflush(file.out);
 }

  int nk1 = 1;

  //calculate_fermi_level(&fermi,eigval,bands,&knet,nkunique,nktot,atoms,job,file);
  calculate_fermi_level_metal_old(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);
  //calculate_fermi_level_metal(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);
  allocate_fermi(&fermi,atoms,job,file);
  //allocate_fermi(&fermi,&knet.unique,&nk1,atoms,job,file);
  occupation_fermi(&fermi, &knet, eigval, knet.unique, nbands, job, file);
  //CHANGES2014occupation_fermi(&fermi, eigval, knet.unique, nbands, job, file);

  if (job->verbosity >= 1) {
  count = 0;
  //int working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) / 2 ;
  //int target_states = 2 * working_value * nktot;
  int working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2);
  int target_states = working_value * nktot;

  for(i=0;i<nkunique * job->spin_dim;i++) {
//fermi.occupied[i] = 179;
  //for(j=0;j<nbands;j++) {
  count += fermi.occupied[i] * knet.num[i] * job->spin_fac;
  //fprintf(file.out,"fermi %d %d %d %lf\n",i,j,fermi.occupied[i],fermi.occupation[i*nbands+j]);}
  fprintf(file.out,"fermi %d %d\n",i,fermi.occupied[i]);
  }
  fprintf(file.out,"occupied states %d target states %d total electrons %d\n",count, target_states, count/nktot + 2 * bands[0] - 2);
  }

   nq = 1;

   for (q = 0; q < nq; q++) {
    for (j = begin_k1[job->taskid]; j < end_k1[job->taskid]; j++) {

    //int ntransitions = (nbands - fermi.occupied[knet.fbz[j]]) * fermi.occupied[knet.fbz[j]];
    //if (ntransitions <= 0) {
    //fprintf(file.out,"Error in susceptibility. Fermi energy not in range of bands. ntransitions = %d\n",ntransitions);
   //}

    //double E_vertex[polygons][vertices][job->spin_dim * ntransitions];
    //double matrix_element_fac[job->field_dirs][ntransitions];

    k_point = decompose_k_point(is,j,crystal,job,file);
    q_point = decompose_k_point(is,0,crystal,job,file);
    p = compose_k_point(is,k_point,q_point,crystal,job,file);

    tetrahedron_vertices2(is, j,  k_vertex, &knet, crystal, job, file);
    tetrahedron_vertices2(is, j, kq_vertex, &knet, crystal, job, file);
    if (job->verbosity >= 1) {
    fprintf(file.out,"j %d k_point %3d %3d %3d k_vertex %3d kq_vertex %d\n",j,k_point.comp1,k_point.comp2,k_point.comp3,k_vertex[0][0], kq_vertex[0][0]);
    fflush(file.out);
   }

    //nk[0] = k_vertex[0][0];
    //nk[1] = k_vertex[0][0];
    nk[0] = j;
    nk[1] = j;

    fourier_transform_3(&one_ints.Momentum[0], &M_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
    //fourier_transform_3(&one_ints.Momentum[0], &M_k->a[0][0], &knet, nk, &pair_p, num_p, R, atoms, shells, symmetry, job, file);

    //fprintf(file.out,"M_k\n");
    //print_complex_matrix(M_k,file);

    for (i = 0; i < job->spin_dim; i++) {

    int ntransitions = (nbands - fermi.occupied[i * nkunique + knet.fbz[j]]) * fermi.occupied[i * nkunique + knet.fbz[j]];
    if (ntransitions <= 0) {
    fprintf(file.out,"Error in susceptibility. Fermi energy not in range of bands. ntransitions = %d\n",ntransitions);
   }

    double E_vertex[polygons][vertices][job->spin_dim * ntransitions];
    double matrix_element_fac[job->field_dirs][ntransitions];

    MPI_File_seek(fh, (i * nkunique +  k_vertex[0][0]) * block_size, MPI_SEEK_SET) ;
    MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    // CHECK new MPI_File_seek(fh, (i * nkunique +  j) * block_size, MPI_SEEK_SET) ;
    // CHECK MPI_File_read(fh, &eigvec1->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;

        //bands[0] = fermi->bands[2 * i] - 1;
        //bands[1] = fermi->bands[2 * i + 1] - 1;
        //CHANGE2014 fix for spin pol case
        nbands = bands[1] - bands[0] + 1;

    rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],nbands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);
    //CHANGE2014rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],bands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);

    //if (knet.trs[j] == 1) {
    //fprintf(file.out,"\neigvec0\n");
    //print_complex_matrix(eigvec0,file);
    //print_complex_matrix1(eigvec0,1,file);
    //}
    //fprintf(file.out,"\neigvec1\n");
    //print_complex_matrix(eigvec1,file);
    //print_complex_matrix1(eigvec1,1,file);

    MPI_File_seek(fh, (i * nkunique + kq_vertex[0][0]) * block_size, MPI_SEEK_SET) ;
    MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
    // CHECK new MPI_File_seek(fh, (i * nkunique +  j) * block_size, MPI_SEEK_SET) ;
    // CHECK MPI_File_read(fh, &eigvec2->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;

    rotate_psi(&eigvec0->a[0][0],&eigvec2->a[0][0],nbands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);
    //CHANGE2014rotate_psi(&eigvec0->a[0][0],&eigvec2->a[0][0],bands,j,&knet,atom_p,atoms,R,shells,symmetry,job,file);

    for (int i2 = 0; i2 < nbands; i2++) {
      for (int j2 = 0; j2 < nbfn; j2++) {
        eigvec1->a[i2][j2].imag() *= -k_one;
      }
     }

    //fprintf(file.out,"\neigvec2\n");
    //print_complex_matrix(eigvec1,file);
    //fprintf(file.out,"\neigvec1\n");
    //print_complex_matrix1(eigvec1,1,file);
    //fprintf(file.out,"\neigvec2\n");
    //print_complex_matrix1(eigvec2,1,file);

    //m1 = atoms->number_of_sh_bfns_in_unit_cell;
    //k1 = m1;
    //n1 = nbands;
    //lda1 = m1;
    //ldb1 = m1;
    //ldc1 = n1;
    //ComplexGEMM2(&NoTrans,&Trans,&m1,&n1,&k1,&alpha,&M_k,&lda1,&eigvec2,&ldb1,&beta,&tmp,&ldc1);
    //m1 = nbands;
    //n1 = nbands;
    //k1 = atoms->number_of_sh_bfns_in_unit_cell;
    //lda1 = k1;
    //ldb1 = nbands;
    //ldc1 = n1;
    //ComplexGEMM2(&NoTrans,&Trans,&m1,&n1,&k1,&alpha,&M_k,&lda1,&eigvec2,&ldb1,&beta,&tmp,&ldc1);
    //cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,m1,n1,k1,&alpha,&M_k->a[0][0],lda1,&eigvec2->a[0][0],ldb1,&beta,&tmp->a[0][0],ldc1);
    //m1 = nbands; ldb1 = nbands;
    //ComplexGEMM2(&NoTrans,&NoTrans,&m1,&n1,&k1,&alpha,&eigvec1,&lda1,&tmp,&ldb1,&beta,&M_x,&ldc1);
  //cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m1,n1,k1,&alpha,&eigvec1->a[0][0],lda1,&tmp->a[0][0],ldb1,&beta,&M_x->a[0][0],ldc1);

    for (j3 = 0; j3 < 3; j3++) {
    dim = j3 * atoms->number_of_sh_bfns_in_unit_cell;
    dim5 = j3 * nbands;
    ComplexGEMM3(&NoTrans,&Trans,&nbfn,&nbands,&nbfn,&alpha,&(M_k->a[dim]),&nbfn,&(eigvec2->a[0]),&nbfn,&beta,&(tmp->a[0]),&nbands);
    ComplexGEMM3(&NoTrans,&NoTrans,&nbands,&nbands,&nbfn,&alpha,&(eigvec1->a[0]),&nbfn,&(tmp->a[0]),&nbands,&beta,&(M_x->a[dim5]),&nbands);
   }

/*
        fprintf(file.out,"EIGENVEC2\n");

      for (int i2 = 0; i2 < 4; i2++) {
      for (int j2 = 0; j2 < nbfn; j2++) {
        fprintf(file.out,"%3d %3d %10.4lf %10.4lf\n",i2,j2,eigvec2->a[i2][j2].real(),eigvec2->a[i2][j2].imag()) ;
      }}

        fprintf(file.out,"M_k\n");

      for (int i2 = 0; i2 < 4; i2++) {
      for (int j2 = 0; j2 < nbfn; j2++) {
        fprintf(file.out,"%3d %3d %10.4lf %10.4lf\n",i2,j2,M_k->a[i2][j2].real(),M_k->a[i2][j2].imag()) ;
      }}

        fprintf(file.out,"tmp\n");

      for (int i2 = 0; i2 < 4; i2++) {
      for (int j2 = 0; j2 < nbfn; j2++) {
        fprintf(file.out,"%3d %3d %10.4lf %10.4lf\n",i2,j2,tmp->a[i2][j2].real(),tmp->a[i2][j2].imag()) ;
      }}

        fprintf(file.out,"M_x\n");

      for (int i2 = 0; i2 < 4; i2++) {
      for (int j2 = 0; j2 < nbfn; j2++) {
        fprintf(file.out,"%3d %3d %10.4lf %10.4lf\n",i2,j2,M_x->a[i2][j2].real(),M_x->a[i2][j2].imag()) ;
      }}
*/

    ////ComplexGEMM2(&NoTrans,&Trans,&nbfn,&nbands,&nbfn,&alpha,&M_k,&nbfn,&eigvec2,&nbfn,&beta,&tmp,&nbands);
    ////ComplexGEMM2(&NoTrans,&NoTrans,&nbands,&nbands,&nbfn,&alpha,&eigvec1,&nbfn,&tmp,&nbands,&beta,&M_x,&nbands);

    //fprintf(file.out,"M_k\n");
    //print_complex_matrix(M_k,file);
    //fprintf(file.out,"\nM_x\n");
    //print_complex_matrix(M_x,file);

    //if (job->verbosity > 1) {
    if (job->verbosity > 1) {
    fprintf(file.out,"M_k\n");
    print_complex_matrix(M_k,file);

    fprintf(file.out,"\neigvec1\n");
    print_complex_matrix(eigvec1,file);

    fprintf(file.out,"\neigvec2\n");
    print_complex_matrix(eigvec2,file);

    fprintf(file.out,"tmp\n");
    print_complex_matrix(tmp,file);

    fprintf(file.out,"\nM_x\n");
    print_complex_matrix(M_x,file);
   }
 
        for (j3 = 0; j3 < job->field_dirs; j3++) {
        count = 0;
         for (int i1 = 0; i1 < fermi.occupied[i * nkunique + knet.fbz[j]]; i1++) {
          for (int j1 = fermi.occupied[i * nkunique + knet.fbz[j]]; j1 < nbands; j1++) {
            matrix_element_fac[j3][count] = (field[j3].comp1 * M_x->a[i1][j1].real() + \
                                             field[j3].comp2 * M_x->a[nbands + i1][j1].real() + \
                                             field[j3].comp3 * M_x->a[2 * nbands + i1][j1].real()) * \
                                            (field[j3].comp1 * M_x->a[i1][j1].real() + \
                                             field[j3].comp2 * M_x->a[nbands + i1][j1].real() + \
                                             field[j3].comp3 * M_x->a[2 * nbands + i1][j1].real()) + \
                                            (field[j3].comp1 * M_x->a[i1][j1].imag() + \
                                             field[j3].comp2 * M_x->a[nbands + i1][j1].imag() + \
                                             field[j3].comp3 * M_x->a[2 * nbands + i1][j1].imag()) * \
                                            (field[j3].comp1 * M_x->a[i1][j1].imag() + \
                                             field[j3].comp2 * M_x->a[nbands + i1][j1].imag() + \
                                             field[j3].comp3 * M_x->a[2 * nbands + i1][j1].imag());
            //matrix_element_fac[j3][count] = k_one;
           //fprintf(file.out,"%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n\
           %3d %3d %3d %e %e\n%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n\
           %3d %3d %3d %e %e\n%3d %3d %3d %e %e\n%3d %3d %3d %e %e\n", \
           j3,i1,j1,field[j3].comp1 , M_x->a[i1][j1].real(), \
           j3,i1,j1,field[j3].comp2 , M_x->a[nbands + i1][j1].real(), \
           j3,i1,j1,field[j3].comp3 , M_x->a[2 * nbands + i1][j1].real(), \
           j3,i1,j1,field[j3].comp1 , M_x->a[i1][j1].real() , \
           j3,i1,j1,field[j3].comp2 , M_x->a[nbands + i1][j1].real(), \
           j3,i1,j1,field[j3].comp3 , M_x->a[2 * nbands + i1][j1].real(), \
           j3,i1,j1,field[j3].comp1 , M_x->a[i1][j1].imag() , \
           j3,i1,j1,field[j3].comp2 , M_x->a[nbands + i1][j1].imag(), \
           j3,i1,j1,field[j3].comp3 , M_x->a[2 * nbands + i1][j1].imag(), \
           j3,i1,j1,field[j3].comp1 , M_x->a[i1][j1].imag(), \
           j3,i1,j1,field[j3].comp2 , M_x->a[nbands + i1][j1].imag(), \
           j3,i1,j1,field[j3].comp3 , M_x->a[2 * nbands + i1][j1].imag());
           //if (j3 == 0) fprintf(mat, "%d %d %e %e %e\n",k_point.comp1,k_point.comp2,knet.cart[j].comp1,knet.cart[j].comp2,matrix_element_fac[0][count]);
           //if (j3 == 1) fprintf(mat1,"%d %d %e %e %e\n",k_point.comp1,k_point.comp2,knet.cart[j].comp1,knet.cart[j].comp2,matrix_element_fac[1][count]);
           //if (j3 == 0) fprintf(mat2,"%d %d %e %e %e\n",k_point.comp1,k_point.comp2, knet.cart[j].comp1,knet.cart[j].comp2, \
           *(eigval + kq_vertex[0][0] * nbands + 1) - *(eigval + k_vertex[0][0] * nbands));
           count++;
          }
         }
        }

          //double shift = 0.018375; // 0.5 eV
          //double shift = 0.011025; // 0.3 eV

          dim1 = i * nkunique * nbands;
           for (p = 0; p < polygons; p++) {
            for (v = 0; v < vertices ; v++) {
             count = 0;
             for (m = 0; m < fermi.occupied[i * nkunique + knet.fbz[j]]; m++) {
               for (n = fermi.occupied[i * nkunique + knet.fbz[j]]; n < nbands; n++) {
                 //E_vertex[p][v][count] = shift + *(eigval + dim1 + kq_vertex[p][v] * nbands + n) - \
                 *(eigval + dim1 + k_vertex[p][v] * nbands + m);
                   E_vertex[p][v][count] = *(eigval + dim1 + kq_vertex[p][v] * nbands + n) - *(eigval + dim1 + k_vertex[p][v] * nbands + m);
                   //fprintf(file.out,"E_vertex %d %d %d %d %d %lf %lf %lf %lf %lf\n",p,v,m,n,kq_vertex[p][v], \
                   *(eigval + dim1 + kq_vertex[p][v] * nbands + n),\
                   *(eigval + dim1 + k_vertex[p][v] * nbands + m),E_vertex[p][v][count] ,matrix_element_fac[0][count], \
                   matrix_element_fac[1][count]);
if (E_vertex[p][v][count] > 1.45/27.211 && E_vertex[p][v][count] < 1.55/27.211 && (matrix_element_fac[0][count] > 6e-02 || matrix_element_fac[1][count] > 6e-02))  fprintf(file.out,"point1 %3d %3d %3d %3d %e %e %e\n",j, m, n, count,E_vertex[p][v][count],\
 matrix_element_fac[0][count],matrix_element_fac[1][count]);
if (E_vertex[p][v][count] > 1.70/27.211 && E_vertex[p][v][count] < 1.80/27.211 && (matrix_element_fac[0][count] > 3e-02 || matrix_element_fac[1][count] > 3e-02))  fprintf(file.out,"point2 %3d %3d %3d %3d %e %e %e\n",j, m, n, count,E_vertex[p][v][count], \
 matrix_element_fac[0][count],matrix_element_fac[1][count]);
                   count++;
                  }
                 }
                }
               }

           for (p = 0; p < polygons; p++) {
            for (v = 0; v < vertices - 1; v++) {
             for (k = 0; k < vertices - 1; k++) {
              for (n = 0; n < ntransitions; n++) {
                if (E_vertex[p][k][n] < E_vertex[p][k + 1][n]) {
                  E_tmp = E_vertex[p][k][n];
                  E_vertex[p][k][n] = E_vertex[p][k + 1][n];
                  E_vertex[p][k+1][n] = E_tmp;
                }
               }
              }
             }
            }

             //for (n = 0; n < ntransitions; n++) {
             //for (p = 0; p < polygons; p++) {
             //for (v = 0; v < vertices; v++) {
             //fprintf(file.out,"check %d %d %d %lf\n",p,v,n,E_vertex[p][v][n]);
             //}}
             //fprintf(file.out,"\n");
             //}

   for (p = 0; p < polygons; p++) {
    for (v = 0; v < vertices; v++) {
     for (n = 0; n < ntransitions; n++) {

    start = (int) ((E_vertex[p][vertices - 1][n] - energy_range[0]) / increment);
    finish = (int) ((E_vertex[p][0][n] - energy_range[0]) / increment) + 1;
  //fprintf(file.out,"start %d finish %d p v n %d %d %d E_vertex %lf %lf %lf %lf\n",start,finish,p,v,n,E_vertex[p][0][n],E_vertex[p][1][n],\
    E_vertex[p][2][n], increment * (double) finish + energy_range[0]);
    //start = 0; finish = npoints;
    if (start < 0)
      start = 0;
    if (start > npoints)
      continue;
    if (finish < 0)
      continue;
    if (finish > npoints)
      finish = npoints;

    for (k = start; k <= finish; k++) {
      E = increment * (double) k + energy_range[0];
      fac = k_zero;
      fac_real = k_zero;
      fac_imag = k_zero;

      switch (crystal->type[0]) {

        case 'C':

     //if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {

/*
      //if (E_vertex[p][3][n] < E_vertex[p][2][n] && E_vertex[p][2][n] < E_vertex[p][1][n] && E_vertex[p][1][n] < E_vertex[p][0][n]) {
      if (E_vertex[p][3][n] < E_vertex[p][2][n] && E_vertex[p][2][n] < E_vertex[p][1][n] && E_vertex[p][1][n] < E_vertex[p][0][n]) {

     //fac = (E + E_vertex[p][0][n]) * (E + E_vertex[p][0][n]) * log(fabs((E_vertex[p][0][n] + E) / (E_vertex[p][3][n] + E))) /
     fac = (E + E_vertex[p][0][n]) * (E + E_vertex[p][0][n]) * log(fabs((E + E_vertex[p][0][n]) /(E + E_vertex[p][3][n]))) /
             (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]) ;
     //fprintf(file.out,"%d %lf %lf %lf %lf\n",k,E,E_vertex[p][0][n],E_vertex[p][1][n],log(fabs((E + E_vertex[p][0][n]) /(E + E_vertex[p][1][n]))));
       //(E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
              //- E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]);

// *
     fac += (E - E_vertex[p][1][n]) * (E - E_vertex[p][1][n]) * log(fabs((E_vertex[p][1][n] - E) / (E_vertex[p][3][n] - E))) /
       (E_vertex[p][1][n] - E_vertex[p][3][n]) / (E_vertex[p][1][n]
              - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][0][n]);

     fac += (E - E_vertex[p][2][n]) * (E - E_vertex[p][2][n]) * log(fabs((E_vertex[p][2][n] - E) / (E_vertex[p][3][n] - E))) /
       (E_vertex[p][2][n] - E_vertex[p][3][n]) / (E_vertex[p][2][n]
              - E_vertex[p][1][n]) / (E_vertex[p][2][n] - E_vertex[p][0][n]);
 * /

             }

     //fac += (E - E_vertex[p][1][n]) * (E - E_vertex[p][1][n]) * log(fabs((E_vertex[p][1][n] - E) / (E_vertex[p][3][n] - E))) /
       //(E_vertex[p][1][n] - E_vertex[p][0][n]) / (E_vertex[p][1][n]
              //- E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]);
//}

// *
     //fac = (E - E_vertex[p][0][n]) * (E - E_vertex[p][0][n]) / \
*/
 

      if (E > E_vertex[p][3][n] && E < E_vertex[p][2][n]) {
      fac = (E - E_vertex[p][3][n]) * (E - E_vertex[p][3][n]) / (E_vertex[p][2][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
                - E_vertex[p][3][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]);
      //fprintf(file.out,"fac k %d %lf %lf %lf\n",k,E,fac,matrix_element_fac[n]);
      //fprintf(file.out,"fac k %d %lf %lf\n",k,E,fac);
      //fflush(file.out);
          }

      if (E > E_vertex[p][2][n] && E < E_vertex[p][1][n]) {
      fac = -(E - E_vertex[p][1][n]) * (E - E_vertex[p][3][n]) / (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][1][n]
                - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][3][n]) -\
                 (E - E_vertex[p][0][n]) * (E - E_vertex[p][2][n])
               / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][2][n]) / \
                 (E_vertex[p][0][n] - E_vertex[p][3][n]);
      //fprintf(file.out,"fac k %d %lf %lf\n",k,E,fac);
      //fprintf(file.out,"fac k %d %lf %lf %lf\n",k,E,fac,matrix_element_fac[n]);
      //fflush(file.out);
          }

     if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {
     fac = (E - E_vertex[p][0][n]) * (E - E_vertex[p][0][n]) / (E_vertex[p][0][n] - E_vertex[p][3][n]) / (E_vertex[p][0][n]
              - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]);
      //fprintf(file.out,"fac k %d %lf %lf\n",k,E,fac);
      //fprintf(file.out,"fac k %d %lf %lf %lf\n",k,E,fac,matrix_element_fac[n]);
      //fflush(file.out);
          }
//  * /

          break;

        case 'S':

     if (E > E_vertex[p][2][n] && E < E_vertex[p][1][n]) {
     fac = (E - E_vertex[p][2][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]) / (E_vertex[p][1][n] - E_vertex[p][2][n]);
          }

     if (E > E_vertex[p][1][n] && E < E_vertex[p][0][n]) {
     fac = -(E - E_vertex[p][0][n]) / (E_vertex[p][0][n] - E_vertex[p][1][n]) / (E_vertex[p][0][n] - E_vertex[p][2][n]);
          }
      //fprintf(file.out,"fac k %d %lf %lf   %lf %lf %lf   %lf %lf\n",k,E,fac,E_vertex[p][0][n],E_vertex[p][1][n], \
      //E_vertex[p][2][n],matrix_element_fac[0][n],matrix_element_fac[1][n]);

          break;

        case 'P':

          fprintf(file.out, "1-D DOS code not complete\n");
          exit(1);

          break;

      } // close switch (crystal->type

         for (j3 = 0; j3 < job->field_dirs; j3++) {
           //fprintf(file.out,"matrix elem %d %d %d %d %e %e  %e %e %e\n",p,n,k,j3,fac,matrix_element_fac[j3][n], \
           E_vertex[p][0][n],E_vertex[p][1][n],E_vertex[p][2][n]);
           spectrum[i][2 * j3][k] += fac * vol * matrix_element_fac[j3][n] / E / E;
          }

              } // close loop over k
             } // close loop over n
            }  // close loop over v
           } // close loop over p
          } // close loop on i

          } // close loop over j
         } // close loop on q

      for (j = 0; j < job->spin_dim; j++) {
      for (j3 = 0; j3 < 2 * job->field_dirs; j3++) { // factor of 2 is for real, imaginary parts
      MPI_Reduce(&spectrum[j][j3][0],&spectrum1[j][j3][0],npoints + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     }
    }

 for (j3 = 0; j3 < job->field_dirs; j3++) {
  for (i = 0; i < npoints; i++) {
    E_i = increment * (double) i + energy_range[0];
      if (E_i < 0.0001) continue;
      if ((i/2)*2 == i) {
        for (j = 1; j < npoints; j+=2) {
          E_j = increment * (double) j + energy_range[0];
          spectrum1[0][2 * j3 + 1][i] += two * increment / pi * spectrum1[0][2 * j3][j] * (k_one / (E_j - E_i) + k_one / (E_j + E_i));
         }
        }
      else {
        for (j = 0; j < npoints; j+=2) {
          E_j = increment * (double) j + energy_range[0];
          spectrum1[0][2 * j3 + 1][i] += two * increment / pi * spectrum1[0][2 * j3][j] * (k_one / (E_j - E_i) + k_one / (E_j + E_i));
         }
       }
      }
     }

       // trapezoid rule for spectrum sum rule
    
       if (job->taskid == 0) {
        sum_rule_fac = two * epsilon_0 * m0 * crystal->primitive_cell_volume * a0 * a0 * a0 * au_to_eV * au_to_eV / hbar / hbar / pi;
         for (j = 0; j < npoints; j++) {
          energy = energy_range[0] + range * (double) j / double(npoints);
           for (l = 0; l < job->spin_dim; l++) {
            for (j3 = 0; j3 < job->field_dirs; j3++) {
             spectrum1[l][j3 + 2 * job->field_dirs][j] = spectrum1[l][j3 + 2 * job->field_dirs][j - 1] + \
             energy * spectrum1[l][2 * j3][j] * increment * sum_rule_fac;
            }
           }
          }
         }

     if (job->taskid == 0) {
      spect = fopen("optical_spectrum.dat", "w");
       if (spect == NULL) { fprintf(file.out, "cannot open file optical_spectrum.dat\n"); exit(1); }
        for (j = 0; j < npoints; j++) {
         energy = energy_range[0] + range * (double) j / double(npoints);
          fprintf(spect, "%10.4e   ", energy * au_to_eV);
           for (l = 0; l < job->spin_dim; l++) {
            for (j3 = 0; j3 < 3 * job->field_dirs; j3++) {
            //for (j3 = 0; j3 < 2 * job->field_dirs; j3++) {
            fprintf(spect, "%12.4e", spectrum1[l][j3][j]);
            //fprintf(spect, "%12.4e", spectrum1[l][j]);
            //fprintf(spect, "%12.4e", spectrum[l][j]);
           }
          }
         fprintf(spect,"\n");
        }
       fflush(spect);;
       fclose(spect);
      } // end if job->taskid

  DestroyComplexMatrix(&M_k,job);
  DestroyComplexMatrix(&M_x,job);
  DestroyComplexMatrix(&tmp,job);
  DestroyComplexMatrix(&eigvec0,job);
  DestroyComplexMatrix(&eigvec1,job);
  DestroyComplexMatrix(&eigvec2,job);
  free_INT_1E(&one_ints, Function, job, file);
  free(eigval);
  free_k_points(&knet,job);
  //free(knet);
  //free(pair_p);
  //free(pair_p);

  MPI_File_close(&fh) ;
  MPI_File_close(&gh) ;

}

void dielectric_function(int is[3], int bands[2], int npoints, double energy_range[2], ATOM *atoms,  ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int begin_k1[job->numtasks], end_k1[job->numtasks];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int receive_s[job->numtasks], offset_s[job->numtasks];
  int nbands = bands[1] - bands[0] + 1;
  int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
  int block_size = vector_size * sizeof(Complex);
  int value_size = atoms->number_of_sh_bfns_in_unit_cell;
  int i, j, k, l, m, p, q, s, v, n, nq;
  int i1, j1, j3;
  int dim, dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim3, dimf;
  int dim5;
  int dimp, dimg;
  int Function[8] ;
  int count, count_i, count_j;
  double energy, range = energy_range[1] - energy_range[0];
  double time1, time2;
  char buf2[110], buf3[110];
  char zz1[10] = "/evalfile", yy1[10] = "/datafile";
  char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
  FILE *spect;
  PAIR_TRAN pair_p;
  TRIPLE_TRAN Triple;
  KPOINT_TRAN knet;
  INT_1E one_ints;
  FERMI fermi;
  MPI_File fh ;
  MPI_File gh ;
  Complex alpha, beta;
  alpha.real() = k_one;
  alpha.imag() = k_zero;
  beta.real() = k_zero;
  beta.imag() = k_zero;

  // ******************************************************************************************
  // * Routine calculates RPA dielectric function                                             *
  // ******************************************************************************************

        if (job->taskid == 0) {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                 RPA DIELECTRIC FUNCTION CALCULATION                                     |\n");
        fprintf(file.out,"===========================================================================================================\n");
        //fprintf(file.out,"| E  %9.3e to %9.3e | BAND ENERGY RANGE      %9.3e to %9.3e eV | FERMI ENERGY   %8.5lf |\n", \
        energy_range[0] * au_to_eV,energy_range[1] * au_to_eV,*eigval * au_to_eV,*(eigval + nbands - 1) * au_to_eV,job->fermi_energy);
        //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| BAND RANGE   %4d to %4d | OCCUPIED STATES         | TARGET STATES           | TOTAL ELECTRONS         |\n", \
        bands[0],bands[1]);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        //fprintf(file.out,"| MONKHORST-PACK   %2d %2d %2d | UNIQUE PAIRS     %6d | TOTAL PAIRS      %6d | INTERACTION RANGE %5.2lf |\n", \
        is[0],is[1],is[2],pair_p.nump,pair_p.tot,R->cutoff * bohr_to_AA);
        //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        //fprintf(file.out,"|                                     APPLIED FIELD DIRECTION COSINES                                     |\n");
        //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        //for (i = 0; i < job->field_dirs; i++) {
        //fprintf(file.out,"| %3d                       |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        i + 1,field[i].comp1,field[i].comp2,field[i].comp3);
        //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
       //}
      fflush(file.out);
     }

  // ******************************************************************************************
  // * Generate k points                                                                      *
  // ******************************************************************************************

  count_k_points(&knet,is,crystal,symmetry,job,file);
  allocate_k_points(&knet,crystal,job,file);
  if (job->C09 == 1)
  read_XCBD_crystal_09(&knet,crystal,job,file);
  generate_k_points(&knet,is,crystal,symmetry,job,file);
  fermi.nkunique = knet.unique;
  allocate_fermi(&fermi,atoms,job,file);
  fermi.knet = &knet;
  fermi.nktot = knet.nktot;
  fermi.bands[0] = bands[0];
  fermi.bands[1] = bands[1];

  // ******************************************************************************************
  // * open MPI files and load eigenvalues                                                    *
  // ******************************************************************************************

/*
  FILE *scf_evalues;
  char zz4[24] = "scf_evalues";
  char xx1[4];

  AllocateDoubleArray(&eigval,&dim1,job);
  ResetDoubleArray(eigval,&dim1);

  if (crystal->type[0] == 'M') {
  scf_evalues = fopen(zz4, "r");
  fseek(scf_evalues, SEEK_SET, 0);
  fread(eigval, sizeof(double), value_size, scf_evalues);
 }

  if (job->taskid == 0 && job->verbosity >= 1) {
  fprintf(file.out,"Eigenvalues\n");
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < knet.unique; i++) {
      for (j = 0; j < nbands; j++) {
        fprintf(file.out,"%5d %5d %10.4lf\n",i,j,*(eigval + s * knet.unique * nbands + i * nbands + j));
       }
      }
     }
    }
*/

  // ******************************************************************************************
  // * Calculate Fermi level                                                                  *
  // ******************************************************************************************

  //calculate_fermi_level_metal_old(&fermi,eigval,bands,&knet,knet.unique,knet.nktot,atoms,job,file);
 
  // ******************************************************************************************
  // * Read parameters in new_density_matrix                                                  *
  // ******************************************************************************************

  double *P0;
  int dimp_read = 0;
  int dimf_read = 0;
  read_density_matrix(&fermi,&P0,&dimp_read,&dimf_read,atoms,job,file);
  int dimp_read_spin = dimp_read * job->spin_dim;
  int dimf_read_spin = dimf_read * job->spin_dim;
  AllocateDoubleArray(&P0,&dimp_read_spin,job);
  read_density_matrix(&fermi,&P0,&dimp_read,&dimf_read,atoms,job,file);
  DestroyDoubleArray(&P0,&dimp_read_spin,job);

  // ******************************************************************************************
  // * Count and generate atom pairs for DIELECTRIC FUNCTION calculation                      *
  // * Count size of reduced (dimp) and full (dimf) atom pair array                           *
  // ******************************************************************************************
 
  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity > 1) print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Count size of reduced (dimp) and full (dimf) atom pair array                           *
  // ******************************************************************************************

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,&pair_p,atoms,job,file);
  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file);

  // ******************************************************************************************
  // * Calculate Overlap and Coulomb matrices                                                 *
  // ******************************************************************************************
 
  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 1 ;

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  allocate_INT_1E(&one_ints, dim, Function, job, file);

  INT_1E one_ints_buffer;
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  time1 = MPI_Wtime();
  //fock_element_1e(&one_ints, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  fock_element_1e1(&one_ints_buffer, dim, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  time2 = MPI_Wtime() - time1;
  MPI_Allreduce(&one_ints_buffer.Overlap[0],&one_ints.Overlap[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&one_ints_buffer.Coulomb[0],&one_ints.Coulomb[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);

  // ******************************************************************************************
  // * Calculate Polarisability and Dielectric Function                                       *
  // ******************************************************************************************

  int dimtp, dimtf;
  int dimp1 = nbfn * nbfn;
  double *overlap_integrals;
  double *Polarisability;

  count_triples1(&Triple, atoms, atom_p, symmetry, R, R_tables, job, file);
  allocate_TRIPLE_TRAN(&Triple, job, file);
  generate_triples1(&Triple, atoms, atom_p, symmetry, R, R_tables, job, file);

  sh_triple_array_dimensions(&dimtp,&dimtf,&Triple,atoms,job,file);
  AllocateDoubleArray(&overlap_integrals,&dimtp,job);
  AllocateDoubleArray(&Polarisability,&dimp1,job);
  ResetDoubleArray(overlap_integrals,&dimtp);
  ResetDoubleArray(Polarisability,&dimp1);

  three_centre_overlap2(overlap_integrals, &Triple, R, atoms, shells, gaussians, crystal, job, file);
  //contract_three_centre_overlap(Polarisability,overlap_integrals,&fermi,&pair_p,&Triple,crystal,atoms,shells,symmetry,R,job,file);
 
  // ******************************************************************************************
  // * Loop over k points                                                                     *
  // ******************************************************************************************
 
  int nk[2];
  int begin_k, end_k;
  int dim2;
  int info = 0;
  char uplo = 'U';
  char jobz = 'V';

  begin_k = 0;
  end_k = 1;

  for (k = begin_k; k < end_k; k++) {

  nk[0] = knet.ibz[k];
  nk[1] = knet.ibz[k];

  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  dim2 = dim1 * dim1;

  // ******************************************************************************************
  // * Find transformation to orthogonalise basis at this k point and store in xtrn1          *
  // ******************************************************************************************
 
  double *eigenvalues;
  ComplexMatrix *S_k, *eigenvectors;
  ComplexMatrix *xtrn1, *xtmp1, *xtmp11;

  AllocateComplexMatrix(&S_k,&dim1,&dim1,job);
  AllocateComplexMatrix(&eigenvectors,&dim1,&dim1,job);
  AllocateDoubleArray(&eigenvalues,&dim1,job);
  ResetDoubleArray(eigenvalues,&dim1);
  ResetComplexMatrix(eigenvectors);
  ResetComplexMatrix(S_k);

  fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
  DiagonaliseHermitian(&S_k, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);

  int n = 0;
  for (int ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
    if (eigenvalues[ii] < 1.0e-04 && job->taskid == 0)  fprintf(file.out,"small eigenvalue %d %d %e\n",k,ii,eigenvalues[ii]);
    if (eigenvalues[ii] < 1.0e-04)  n++;
  }
  int dim11 = dim1 - n;
  printf("%3d %3d %3d\n",dim11,dim1,k);

  AllocateComplexMatrix(&xtrn1,&dim11,&dim1,job);
  AllocateComplexMatrix(&xtmp1,&dim1,&dim1,job);
  AllocateComplexMatrix(&xtmp11,&dim11,&dim1,job);

  //if (job->kpoints == 0 && dim11 < fermi.occupied[k]) {
  //if (job->taskid == 0)
  //fprintf(file.out,"Limited basis set is also linearly dependent at k-point %3d. Stopping.\n",k);
  //MPI_Finalize();
  //exit(1);
  //}

  for (i = n; i < dim1; i++) {
    for (j = 0; j < dim1; j++) {
      xtrn1->a[i - n][j] = eigenvectors->a[i][j] / sqrt(*(eigenvalues + i));
     }
    }

  DestroyDoubleArray(&eigenvalues,&dim1,job);
  DestroyComplexMatrix(&eigenvectors,job);

  // ******************************************************************************************
  // * Check that overlap matrix is correctly Fourier transformed                             *
  // ******************************************************************************************

  if (job->taskid == 0 && job->verbosity > 1) {
    ComplexMatrix *S_k1, *S_k2;
    AllocateComplexMatrix(&S_k1,&dim1,&dim1,job);
    AllocateComplexMatrix(&S_k2,&dim11,&dim11,job);
    fourier_transform(&one_ints.Overlap[0], &S_k1->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
    ResetComplexMatrix(xtmp11);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &S_k1, &beta, &xtmp11);
    ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp11, &xtrn1, &beta, &S_k2);
    fprintf(file.out,"S_k1 orthogonalised %d\n",k);
    print_complex_matrix(S_k2, file);
    DestroyComplexMatrix(&S_k1,job);
    DestroyComplexMatrix(&S_k2,job);
  }

  // ******************************************************************************************
  // * Fourier Transform V and P and Calculate Dielectric Matrix D_k                          *
  // ******************************************************************************************

  ComplexMatrix *D_k, *P_k, *V_k, *VSP_k;

  AllocateComplexMatrix(&D_k,&dim1,&dim1,job);
  AllocateComplexMatrix(&P_k,&dim1,&dim1,job);
  AllocateComplexMatrix(&V_k,&dim1,&dim1,job);
  AllocateComplexMatrix(&VSP_k,&dim1,&dim1,job);
  ResetComplexMatrix(D_k);
  ResetComplexMatrix(V_k);
  ResetComplexMatrix(S_k);
  ResetComplexMatrix(VSP_k);

  fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
  fourier_transform(&one_ints.Coulomb[0], &V_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
  //fourier_transform(&Polarisability[0], &P_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
 
  for (i = 0; i < dim1; i++) {
    for (j = 0; j < dim1; j++) {
      P_k->a[i][j] = Polarisability[i * atoms->number_of_sh_bfns_in_unit_cell + j];
     }
    }

  fprintf(file.out,"P_k\n");
  print_complex_matrix(P_k, file);

  ResetComplexMatrix(xtmp1);
  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &V_k, &S_k, &beta, &xtmp1);
  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &P_k, &beta, &VSP_k);

  for (i = 0; i < dim1; i++) {
    for (j = 0; j < dim1; j++) {
      D_k->a[i][j] = S_k->a[i][j] - VSP_k->a[i][j];
     }
    }

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"V_k\n");
  print_complex_matrix(V_k, file);
  fprintf(file.out,"P_k\n");
  print_complex_matrix(P_k, file);
  fprintf(file.out,"D_k\n");
  print_complex_matrix(D_k, file);
 }

  //DestroyComplexMatrix(&V_k, job);
  DestroyComplexMatrix(&S_k, job);
  DestroyComplexMatrix(&VSP_k, job);

  // ******************************************************************************************
  // * Transform Dielectric Matrix to Orthogonal Basis and Invert                             *
  // ******************************************************************************************

  ComplexMatrix *D_k1, *D_k1_inv;
  AllocateComplexMatrix(&D_k1,&dim11,&dim11,job);
  AllocateComplexMatrix(&D_k1_inv,&dim11,&dim11,job);
  ResetComplexMatrix(xtmp11);
  ResetComplexMatrix(D_k1);
  ResetComplexMatrix(D_k1_inv);

  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &D_k, &beta, &xtmp11);
  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp11, &xtrn1, &beta, &D_k1);

  DestroyComplexMatrix(&D_k, job);

  int ss;
  gsl_matrix_complex *M = gsl_matrix_complex_alloc (dim11,dim11), *I = gsl_matrix_complex_alloc (dim11,dim11);
  gsl_permutation *pp = gsl_permutation_alloc (dim11);

  for (i = 0; i < dim11; i++) {
    for (j = 0; j < dim11; j++) {
      gsl_matrix_complex_set(M,i,j,gsl_complex_rect((D_k1->a[i][j]).real(), (D_k1->a[i][j]).imag()));
      gsl_matrix_complex_set(I,i,j,gsl_complex_rect(k_zero, k_zero));
     }
    }

  gsl_linalg_complex_LU_decomp (M, pp, &ss);
  gsl_linalg_complex_LU_invert (M, pp, I);

  for (i = 0; i < dim11; i++) {
    for (j = 0; j < dim11; j++) {
      (D_k1_inv->a[i][j]).real() = GSL_REAL(gsl_matrix_complex_get(I,i,j));
      (D_k1_inv->a[i][j]).imag() = GSL_IMAG(gsl_matrix_complex_get(I,i,j));
     }
    }

  gsl_permutation_free (pp);
  gsl_matrix_complex_free (I);
  gsl_matrix_complex_free (M);

  if (job->taskid == 0 && job->verbosity > 1) {
  ComplexMatrix *D_k4;
  AllocateComplexMatrix(&D_k4,&dim11,&dim11,job);
  ResetComplexMatrix(D_k4);
  for (i = 0; i < dim11; i++) {
    for (j = 0; j < dim11; j++) {
      for (int k1 = 0; k1 < dim11; k1++) {
        D_k4->a[i][k1] += D_k1->a[i][j] * D_k1_inv->a[j][k1];
       }
      }
     }
  fprintf(file.out,"D_k1\n");
  print_complex_matrix(D_k1, file);
  fprintf(file.out,"D_k1_inv\n");
  print_complex_matrix(D_k1_inv, file);
  fprintf(file.out,"D_k1 * D_1_inv\n");
  print_complex_matrix(D_k4, file);
  DestroyComplexMatrix(&D_k4, job);
 }
       
  DestroyComplexMatrix(&D_k1, job);

  // ******************************************************************************************
  // * Transform P_k and V_k Matrices to Orthogonal Basis                                     *
  // ******************************************************************************************

  ComplexMatrix *P_k1, *V_k1;
  AllocateComplexMatrix(&P_k1,&dim11,&dim11,job);
  AllocateComplexMatrix(&V_k1,&dim11,&dim11,job);
  ResetComplexMatrix(xtmp11);
  ResetComplexMatrix(P_k1);
  ResetComplexMatrix(V_k1);

  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &P_k, &beta, &xtmp11);
  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp11, &xtrn1, &beta, &P_k1);

  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &V_k, &beta, &xtmp11);
  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp11, &xtrn1, &beta, &V_k1);

  if (job->taskid == 0 && job->verbosity >= 1) {
  fprintf(file.out,"P_k\n");
  print_complex_matrix(P_k, file);
  fprintf(file.out,"P_k1\n");
  print_complex_matrix(P_k1, file);
 }
  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"V_k\n");
  print_complex_matrix(V_k, file);
  fprintf(file.out,"V_k1\n");
  print_complex_matrix(V_k1, file);
 }

  DestroyComplexMatrix(&P_k, job);

  // move S_k freeing three_centre_coulomb_exchange(V_k,S_k,xtrn1,overlap_integrals,&fermi,&pair_p,&Triple,crystal,atoms,shells,symmetry,R,job,file);

  // ******************************************************************************************
  // * Calculate and Diagonalise Response Matrix R_k = D_k1_inv * P_k                         *
  // ******************************************************************************************

  double *eigval;
  ComplexMatrix *R_k1, *W_k1, *eigenvectors1, *eigenvectors2;

  AllocateComplexMatrix(&R_k1,&dim11,&dim11,job);
  AllocateComplexMatrix(&W_k1,&dim11,&dim11,job);
  AllocateComplexMatrix(&eigenvectors1,&dim11,&dim11,job);
  AllocateComplexMatrix(&eigenvectors2,&dim1,&dim1,job);
  AllocateDoubleArray(&eigval,&dim1,job);

  ResetComplexMatrix(R_k1);
  ResetComplexMatrix(W_k1);
  ResetComplexMatrix(xtmp11);
  ResetComplexMatrix(eigenvectors1);
  ResetComplexMatrix(eigenvectors2);
  ResetDoubleArray(eigval,&dim1);

  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &P_k1, &D_k1_inv, &beta, &R_k1);

  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &V_k1, &R_k1, &beta, &xtmp11);
  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp11, &V_k1, &beta, &W_k1);

  if (job->taskid == 0 && job->verbosity >= 1) {
  fprintf(file.out,"R_k1\n");
 }
  if (job->taskid == 0 && job->verbosity > 1) {
  print_complex_matrix(R_k1, file);
  fprintf(file.out,"V_k1\n");
  print_complex_matrix(V_k1, file);
  fprintf(file.out,"W_k1\n");
  print_complex_matrix(W_k1, file);
 }

  DestroyComplexMatrix(&V_k, job);
  DestroyComplexMatrix(&V_k1, job);
  DestroyComplexMatrix(&P_k1, job);

  DiagonaliseHermitian(&R_k1, &eigval, &eigenvectors1, &jobz, &uplo, &info);

  // ******************************************************************************************
  // * Back-transform Diagonalised Response Matrix                                            *
  // ******************************************************************************************

  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors1, &xtrn1, &beta, &eigenvectors2);
  print_complex_eigenvector_matrix2(eigenvectors2, eigval, 5, dim11, k_one, file);

  DestroyComplexMatrix(&R_k1, job);
  DestroyComplexMatrix(&xtmp1, job);
  DestroyComplexMatrix(&xtmp11, job);
  DestroyComplexMatrix(&xtrn1, job);
  DestroyComplexMatrix(&D_k1_inv, job);
  DestroyComplexMatrix(&eigenvectors1, job);
  DestroyComplexMatrix(&eigenvectors2, job);
  DestroyDoubleArray(&eigval,&dim11,job);

  } // close loop on k

  DestroyDoubleArray(&Polarisability,&dimp1,job);
  DestroyDoubleArray(&overlap_integrals,&dimtp,job);
  free_INT_1E(&one_ints, Function, job, file);
  free_k_points(&knet,job);
  free_PAIR_TRAN(&pair_p,job);
  free_TRIPLE_TRAN(&Triple,job);
  free_fermi(&fermi,job);

  //if (crystal->type[0] == 'M' && job->taskid == 0) {
  //fclose(scf_evalues);
 //}

}

