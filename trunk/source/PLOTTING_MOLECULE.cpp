
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <cstring>
#include <cstdlib>
#include <mpi.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "SETUP_SYMMETRY.h"
#include "SETUP_RECIPROCAL_LATTICE.h"
#include "TOOLS.h"
#include "PRINT_UTIL.h"
#include "PARALLEL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "PAIRS_QUADS.h"
#include "DENSITY_MATRIX_MOLECULE.h"
#include "GW_BSE_MOLECULE.h"
#include "PLOTTING_MOLECULE.h"

using namespace std;

void plot_correlated_electron_hole(int grid_par[3], VECTOR_DOUBLE *points, int num_points, VECTOR_DOUBLE *points1, int bse_state[2], FERMI* fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  int i, j, m, n, p, s;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;
  double *bse_eigenvalues;
  DoubleMatrix *scf_eigenvectors, *bse_eigenvectors;
  PAIR_TRAN pair_p;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity > 1)
  print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Allocate memory and generate grid                                                      *
  // ******************************************************************************************

  int grid_size;
  int radmx = 1;
  int number_Rvec;
  int ngrid_points;
  VECTOR_DOUBLE *grid, *Rvec;

  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);
  //spin_grid_size = job->spin_dim * grid_size;

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  calc_grid1(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  // ******************************************************************************************
  // * Read SCF and BSE eigenvectors                                                          *
  // ******************************************************************************************

  int scf_vector_size, bse_vector_size;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim2 = dim1 * dim1;
  int num_bse_states, ntransitions;
  char buf2[120], buf4[120];
  char xx[20] = "/bse_eigenvectors", yy[10] = "/scf_evec";
  MPI_File fh, gh;

  strcpy(buf4,file.scf_eigvec);
  strcat(buf4,xx);
  MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,yy);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&gh) ;

  ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in plot_correlated_electron_hole. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }
  scf_vector_size = job->spin_dim * nbands * dim1;
  num_bse_states = bse_state[1] - bse_state[0] + 1;
  bse_vector_size = num_bse_states * ntransitions;

  AllocateDoubleMatrix(&scf_eigenvectors,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateDoubleMatrix(&bse_eigenvectors,&num_bse_states,&ntransitions,job);

  MPI_File_seek(fh, (bse_state[0] - 1) * ntransitions * sizeof(double), MPI_SEEK_SET);
  MPI_File_read(fh, &bse_eigenvectors->a[0][0], bse_vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE); // read BSE eigenvectors

  MPI_File_seek(gh, dim1 * (fermi->bands[0] - 1) * sizeof(double), MPI_SEEK_SET);
  MPI_File_read(gh, &scf_eigenvectors->a[0][0], scf_vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE); // read SCF eigenvectors

  if (job->verbosity > 1) {
  for (i = 0; i < num_bse_states; i++) {
    for (j = 0; j < ntransitions; j++) {
      fprintf(file.out,"%3d %3d %10.4f\n",i,j,bse_eigenvectors->a[i][j]);
     }
    fprintf(file.out,"\n");
    }

  for (i = 0; i < nbands; i++) {
    for (j = 0; j < dim1; j++) {
      fprintf(file.out,"%3d %3d %10.4f\n",i,j,scf_eigenvectors->a[i][j]);
     }
    fprintf(file.out,"\n");
    }
   }

  // ******************************************************************************************
  // Calculate wavefunctions on grid and electron hole correlations                           *
  // ******************************************************************************************
 
  int count;
  double *psiamp1_buffer, *psiamp2_buffer;
  double *psiamp1, *psiamp2;
  double *Psiamp;
  char buf[22] = "bse_isosurface.xsf";

  VECTOR_DOUBLE Rvec1[radmx], Rvec2[radmx];

  //printf("num_bse %3d nocc %3d nvir %3d bands %3d %3d HOMO %3d bse_state[0] %3d ntrans %3d\n",\
  num_bse_states,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],bse_state[0],ntransitions);

  Rvec1[0].comp1 = 0.0;
  Rvec1[0].comp2 = 0.0;
  Rvec1[0].comp3 = 0.0;
  Rvec2[0].comp1 = 0.0;
  Rvec2[0].comp2 = 0.0;
  Rvec2[0].comp3 = 0.0;

  Psiamp = (double *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(double));
  if (Psiamp == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp \n");
  MPI_Finalize(); exit(1); }

  psiamp1 = (double *) malloc(job->spin_dim * nocc * grid_size * sizeof(double));
  if (psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp2 = (double *) malloc(job->spin_dim * nvir * num_points * sizeof(double));
  if (psiamp2 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2 \n");
  MPI_Finalize(); exit(1); }

  psiamp1_buffer = (double *) malloc(job->spin_dim * nocc * grid_size * sizeof(double));
  if (psiamp1_buffer == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1_buffer \n");
  MPI_Finalize(); exit(1); }

  psiamp2_buffer = (double *) malloc(job->spin_dim * nvir * num_points * sizeof(double));
  if (psiamp2_buffer == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2_buffer \n");
  MPI_Finalize(); exit(1); }

  //fprintf(file.out,"particle grids\n");
  for (n = 0; n < num_points; n++) {
    wavefunction_gridpoint(&psiamp2[n * nvir],&points1[n],Rvec2,&scf_eigenvectors->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
   }

  //fprintf(file.out,"hole grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(psiamp1,&grid[p],Rvec1,&scf_eigenvectors->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
    for (n = 0; n < num_points; n++) {
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < num_bse_states; m++) {
          count = 0;
          for (i = 0; i < nocc; i++) {
            for (j = 0; j < nvir; j++) {
              Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
              bse_eigenvectors->a[m][count] * psiamp1[s * nocc + i] * psiamp2[s * num_points * nvir + n * nvir + j];
              //bse_eigenvector->a[m][count] * psiamp1[s * nocc + i] ;
              //fprintf(file.out,"   %3d %3d %3d %3d %3d %12.6lf %10.4lf %10.4lf %14.8lf   %10.4lf %10.4lf %10.4lf\n",\
              p,n,m,i,j,bse_eigenvectors->a[m][count],psiamp1[i],psiamp2[n * nvir + j],\
              Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p],\
	      grid[p].comp1,grid[p].comp2,grid[p].comp3);
              count++;
             }
            }
           //fprintf(file.out,"\n");
           }
          }
         }
        }

  if (job->taskid == 0)\
  write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,num_bse_states,num_points,crystal,atoms,shells,gaussians,job,file);

  MPI_File_close(&fh);
  MPI_File_close(&gh);

  free(Psiamp);
  free(psiamp1);
  free(psiamp2);
  free(psiamp1_buffer);
  free(psiamp2_buffer);
  //DestroyDoubleArray(&bse_eigenvalues,&ntransitions,job);
  DestroyDoubleMatrix(&bse_eigenvectors,job);
  DestroyDoubleMatrix(&scf_eigenvectors,job);
  free(Rvec);
  free(grid);
  free_PAIR_TRAN(&pair_p,job);

}

void plot_electron_hole_molecule(int grid_par[3], VECTOR_DOUBLE *points, FERMI* fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  int i, j, m, n, p, s;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;
  int ntransitions;
  PAIR_TRAN pair_p;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity > 1)
  print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Allocate memory and generate grid                                                      *
  // ******************************************************************************************

  int grid_size;
  int radmx = 1;
  int number_Rvec;
  int ngrid_points;
  VECTOR_DOUBLE *grid, *Rvec;

  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  calc_grid1(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  // ******************************************************************************************
  // * Read SCF and BSE eigenvectors                                                          *
  // ******************************************************************************************

  int scf_vector_size;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  char buf4[120];
  char xx[10] = "/scf_evec";
  DoubleMatrix *scf_eigenvectors;
  MPI_File fh;

  strcpy(buf4,file.scf_eigvec);
  strcat(buf4,xx);
  MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  //printf("%s\n",buf4);

  ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in plot_electron_hole_molecule. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }
  scf_vector_size = job->spin_dim * nbands * dim1;

  AllocateDoubleMatrix(&scf_eigenvectors,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);

  MPI_File_seek(fh, dim1 * (fermi->bands[0] - 1) * sizeof(double), MPI_SEEK_SET);
  MPI_File_read(fh, &scf_eigenvectors->a[0][0], scf_vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE); // read SCF eigenvectors

  if (job->verbosity > 1) {
  for (i = 0; i < nbands; i++) {
    for (j = 0; j < dim1; j++) {
      fprintf(file.out,"%3d %3d %10.4f\n",i,j,scf_eigenvectors->a[i][j]);
     }
    fprintf(file.out,"\n");
    }
   }

  // ******************************************************************************************
  // Calculate wavefunctions on grid and electron hole correlations                           *
  // ******************************************************************************************
 
  double *psiamp1_buffer, *psiamp2_buffer;
  double *psiamp1, *psiamp2;
  double *Psiamp;
  char yy[16] = "isosurface.xsf";
  char buf1[21] = "hole_";
  char buf2[25] = "electron_";
  char buf3[35] = "transition_density_";
  strcat(buf1,yy);
  strcat(buf2,yy);
  strcat(buf3,yy);

  VECTOR_DOUBLE Rvec1[radmx], Rvec2[radmx];

  Rvec1[0].comp1 = 0.0;
  Rvec1[0].comp2 = 0.0;
  Rvec1[0].comp3 = 0.0;
  Rvec2[0].comp1 = 0.0;
  Rvec2[0].comp2 = 0.0;
  Rvec2[0].comp3 = 0.0;

  Psiamp = (double *) calloc(job->spin_dim * nocc * nvir * grid_size, sizeof(double));
  if (Psiamp == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp \n");
  MPI_Finalize(); exit(1); }

  psiamp1 = (double *) malloc(job->spin_dim * nocc * grid_size * sizeof(double));
  if (psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp2 = (double *) malloc(job->spin_dim * nvir * grid_size * sizeof(double));
  if (psiamp2 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2 \n");
  MPI_Finalize(); exit(1); }

  psiamp1_buffer = (double *) malloc(job->spin_dim * nocc * grid_size * sizeof(double));
  if (psiamp1_buffer == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1_buffer \n");
  MPI_Finalize(); exit(1); }

  psiamp2_buffer = (double *) malloc(job->spin_dim * nvir * grid_size * sizeof(double));
  if (psiamp2_buffer == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2_buffer \n");
  MPI_Finalize(); exit(1); }

  //fprintf(file.out,"hole grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(&psiamp1[p * nocc],&grid[p],Rvec1,&scf_eigenvectors->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
   }

  //fprintf(file.out,"particle grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(&psiamp2[p * nvir],&grid[p],Rvec2,&scf_eigenvectors->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
   }

  for (p = 0; p < grid_size; p++) {
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < nocc; i++) {
        Psiamp[s * nocc * grid_size + i * grid_size + p] = psiamp1[s * grid_size * nocc + p * nocc + i];
       }
      }
     }

  if (job->taskid == 0)
  write_isosurface_xsf(buf1, Psiamp, grid_size, grid_par, points, nocc, 1, crystal, atoms, shells, gaussians, job,file);

  for (p = 0; p < grid_size; p++) {
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < nvir; i++) {
        Psiamp[s * nvir * grid_size + i * grid_size + p] = psiamp2[s * grid_size * nvir + p * nvir + i];
       }
      }
     }

  if (job->taskid == 0)
  write_isosurface_xsf(buf2, Psiamp, grid_size, grid_par, points, nvir, 1, crystal, atoms, shells, gaussians, job,file);

  for (p = 0; p < grid_size; p++) {
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < nocc; i++) {
        for (j = 0; j < nvir; j++) {
          //fprintf(file.out,"   %3d %3d %3d %14.8lf %14.8lf %16.10lf   %10.4lf %10.4lf %10.4lf\n",\
          p,i,j,psiamp1[p * nocc + i],psiamp2[p * nvir + j],\
          Psiamp[i * nvir * grid_size + j * grid_size + p],grid[p].comp1,grid[p].comp2,grid[p].comp3);
          Psiamp[s * nocc * nvir * grid_size + i * nvir * grid_size + j * grid_size + p] = \
          psiamp1[s * grid_size * nocc + p * nocc + i] * psiamp2[s * grid_size  * nvir + p * nvir + j];
         }
        }
       }
      //fprintf(file.out,"\n");
      }

  if (job->taskid == 0)
  write_isosurface_xsf(buf3, Psiamp, grid_size, grid_par, points, nocc, nvir, crystal, atoms, shells, gaussians, job,file);

  free(grid);
  free(Rvec);
  free(Psiamp);
  free(psiamp1);
  free(psiamp2);
  free(psiamp1_buffer);
  free(psiamp2_buffer);
  free_PAIR_TRAN(&pair_p,job);
  DestroyDoubleMatrix(&scf_eigenvectors,job);

  MPI_File_close(&fh);

}

void Rvec_grid_size(int *number_Rvec, int radmx, CRYSTAL *crystal, FILES file)

{

      switch (crystal->type[0]) {

        case 'C':
       
        *number_Rvec = (2 * radmx + 1) * (2 * radmx + 1) * (2 * radmx + 1);

          break;

        case 'S':

        *number_Rvec = (2 * radmx + 1) * (2 * radmx + 1);

          break;

        case 'P':

        *number_Rvec = (2 * radmx + 1);

          break;

        case 'M':
   
        *number_Rvec = 1;

          break;

      } // close switch (crystal->type

}

void generate_Rvec_grid(int number_Rvec, int radmx, VECTOR_DOUBLE *Rvec, CRYSTAL *crystal, FILES file)

{

  int i, j, k;
  VECTOR_DOUBLE *p_Rvec;

      switch (crystal->type[0]) {

        case 'C':

  p_Rvec = Rvec;
  for (i = -radmx; i < (radmx + 1); i++) {
    for (j = -radmx; j < (radmx + 1); j++) {
      for (k = -radmx; k < (radmx + 1); k++) {
        p_Rvec->comp1 = i * crystal->primitive_cell[0].comp1 + j * crystal->primitive_cell[1].comp1 +
                        k * crystal->primitive_cell[2].comp1;
        p_Rvec->comp2 = i * crystal->primitive_cell[0].comp2 + j * crystal->primitive_cell[1].comp2 +
                        k * crystal->primitive_cell[2].comp2;
        p_Rvec->comp3 = i * crystal->primitive_cell[0].comp3 + j * crystal->primitive_cell[1].comp3 +
                        k * crystal->primitive_cell[2].comp3;
        p_Rvec++;
      }
    }
  } // close loops on i, j, k

          break;

        case 'S':

  p_Rvec = Rvec;
  for (i = -radmx; i < (radmx + 1); i++) {
    for (j = -radmx; j < (radmx + 1); j++) {
        p_Rvec->comp1 = i * crystal->primitive_cell[0].comp1 + j * crystal->primitive_cell[1].comp1;
        p_Rvec->comp2 = i * crystal->primitive_cell[0].comp2 + j * crystal->primitive_cell[1].comp2;
        p_Rvec->comp3 = k_zero;
        p_Rvec++;
      }
    } // close loops on i, j

          break;

        case 'P':

  p_Rvec = Rvec;
  for (i = -radmx; i < (radmx + 1); i++) {
        p_Rvec->comp1 = k_zero;
        p_Rvec->comp2 = k_zero;
        p_Rvec->comp3 = i * crystal->primitive_cell[2].comp3;
        p_Rvec++;
       } // close loop on i

          break;

        case 'M':
   
  p_Rvec = Rvec;
        p_Rvec->comp1 = k_zero;
        p_Rvec->comp2 = k_zero;
        p_Rvec->comp3 = k_zero;

          break;

      } // close switch (crystal->type

}

void calc_grid1(int *grid_par, VECTOR_DOUBLE *grid, int grid_size, int *number_of_grid_points, VECTOR_DOUBLE *points, JOB_PARAM *job, FILES file)

{

  int i, j, k;
  VECTOR_DOUBLE vec[3];
  double incrementX[3], incrementY[3], incrementZ[3];
  unsigned char plot_dim = 3; //DEFAULT is a 3D plot

  // define three vectors for plotting
  // NOTE: VECTOR_DOUBLE point[0] is the origin

  for (i = 0; i < 3; i++) {
    switch (grid_par[i]) {
      case 1:
        plot_dim -= 1;
        vec[i].comp1 = 0;
        vec[i].comp2 = 0;
        vec[i].comp3 = 0;
        break;
      default:
        vec[i].comp1 = points[i + 1].comp1;
        vec[i].comp2 = points[i + 1].comp2;
        vec[i].comp3 = points[i + 1].comp3;
    } // end switch
  }

  // Check orthogonalisation of the vectors

  for (i = 0; i < 2; i++) {
    for (j = i + 1; j < 3; j++) {
      if (double_vec_dot(&vec[i], &vec[j]) > 1e-4) {
        fprintf(file.out, "WARNING: Vectors for plotting not orthogonal\n");

      }
    }
  }

  // Calc increments


  for (i = 0; i < 3; i++) {
    switch (grid_par[i]) {
      case 1:
        incrementX[i] = (double) 0;
        incrementY[i] = (double) 0;
        incrementZ[i] = (double) 0;
        break;
      default:
        incrementX[i] = vec[i].comp1 / ((double) (grid_par[i] - 1));
        incrementY[i] = vec[i].comp2 / ((double) (grid_par[i] - 1));
        incrementZ[i] = vec[i].comp3 / ((double) (grid_par[i] - 1));
        if (job->verbosity > 1) {
        fprintf(file.out,"GRID INCREMENTS %3d %10.4lf %10.4lf %10.4lf\n",i,incrementX[i],incrementY[i],incrementZ[i]);
       }
    }
  }

  // Calc grid

  int count = 0;
  //for (i = 0; i < grid_par[0]; i++) {
    //for (j = 0; j < grid_par[1]; j++) {
      //for (k = 0; k < grid_par[2]; k++) {
      for (k = 0; k < grid_par[2]; k++) {
    for (j = 0; j < grid_par[1]; j++) {
  for (i = 0; i < grid_par[0]; i++) {
        grid[count].comp1 = points[0].comp1 + incrementX[0] * i + incrementX[1] * j + incrementX[2] * k;
        grid[count].comp2 = points[0].comp2 + incrementY[0] * i + incrementY[1] * j + incrementY[2] * k;
        grid[count].comp3 = points[0].comp3 + incrementZ[0] * i + incrementZ[1] * j + incrementZ[2] * k;
        if (job->verbosity > 1)
          fprintf(file.out, "Grid components %3d %3d %3d %10.4lf %10.4lf %10.4lf\n", i, j, k, grid[count].comp1,
              grid[count].comp2, grid[count].comp3);
        count++;
        if (count > grid_size) {
          fprintf(file.out, "grid size %d is larger than the buffer %d \n", count, grid_size);
          (*number_of_grid_points) = count;
          return;
        }

      }
    }
  }
  (*number_of_grid_points) = count;

  return;

}

void wavefunction_gridpoint(double *psiamp, VECTOR_DOUBLE *rpoint, VECTOR_DOUBLE *Rvec, double *eigvec, int nstates, int nbands, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  // get psi for nbands and spins on rpoint

  int m, s;
  int index_i, sheli, shelposi, bfposi, gausposi;
  int i1, i2, j2, i3, j3, i4, j4;
  int ip;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  double psiexp, rsqrd;
  double prefac[7];
  double r[3];

  for (i1 = 0; i1 < job->spin_dim * nstates; i1++)
  psiamp[i1] = k_zero;

  //bfposi = 0;
  for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
    r[0] = atoms->cell_vector[ip].comp1 + Rvec->comp1;
    r[1] = atoms->cell_vector[ip].comp2 + Rvec->comp2;
    r[2] = atoms->cell_vector[ip].comp3 + Rvec->comp3;
    rsqrd = (r[0] - rpoint->comp1) * (r[0] - rpoint->comp1) + \
            (r[1] - rpoint->comp2) * (r[1] - rpoint->comp2) + \
            (r[2] - rpoint->comp3) * (r[2] - rpoint->comp3);
    //fprintf(file.out,"rsq %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n", \
    rsqrd,r[0],r[1],r[2],rpoint->comp1,rpoint->comp2,rpoint->comp3);
    //if (rsqrd > 150.0) continue;
    bfposi   = atoms->bfnposn_sh[ip];
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli = shells->type_sh[index_i];
      switch (sheli) {

        case 1:
          prefac[0] =   k_one;
          break;
        case 4:
          prefac[0] =   k_one;
          prefac[1] =   rpoint->comp1 - r[0];
          prefac[2] =   rpoint->comp2 - r[1];
          prefac[3] =   rpoint->comp3 - r[2];
          break;
        case 3:
          prefac[0] =   rpoint->comp1 - r[0];
          prefac[1] =   rpoint->comp2 - r[1];
          prefac[2] =   rpoint->comp3 - r[2];
          break;
        case 5:
          prefac[0] = ((rpoint->comp3 - r[2]) * (rpoint->comp3 - r[2]) -
                       (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) / two -
                       (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) / two) / sqrt(three); // 2 zz - xx - yy
          prefac[1] =  (rpoint->comp1 - r[0]) * (rpoint->comp3 - r[2]); // xz
          prefac[2] =  (rpoint->comp2 - r[1]) * (rpoint->comp3 - r[2]); // yz
          prefac[3] = ((rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) -
                       (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1])) / two; // xx - yy
          prefac[4] =  (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]); // xy
          break;
        case 7:
          prefac[0] =  (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]) *  three/two/sqrt(six) +
                       (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * -k_one/two/sqrt(six);
          prefac[1] =  (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]) * (rpoint->comp3 - r[2]);
          prefac[2] =  (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]) * -k_one/two/sqrt(ten) +
                       (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * -k_one/two/sqrt(ten) +
                       (rpoint->comp2 - r[1]) * (rpoint->comp3 - r[2]) * (rpoint->comp3 - r[2]) *  four/two/sqrt(ten);
          prefac[3] =  (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * (rpoint->comp3 - r[2]) * -three/two/sqrt(fifteen) +
                       (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * (rpoint->comp3 - r[2]) * -three/two/sqrt(fifteen) +
                       (rpoint->comp3 - r[2]) * (rpoint->comp3 - r[2]) * (rpoint->comp3 - r[2]) *  k_one/sqrt(fifteen);
          prefac[4] =  (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * -k_one/two/sqrt(ten) +
                       (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * -k_one/two/sqrt(ten) +
                       (rpoint->comp1 - r[0]) * (rpoint->comp3 - r[2]) * (rpoint->comp3 - r[2]) *  four/two/sqrt(ten);
          prefac[5] =  (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * (rpoint->comp3 - r[2]) *  k_one/two +
                       (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * (rpoint->comp3 - r[2]) * -k_one/two;
          prefac[6] =  (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) *  k_one/two/sqrt(six) +
                       (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1]) * -three/two/sqrt(six);
          break;
   
      } // close switch
    psiexp = k_zero;
    for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      psiexp += gaussians->c_sh[gausposi + i4] * exp(- gaussians->expo_sh[gausposi + i4] * rsqrd);
      //fprintf(file.out,"psiexp %3d %3d %10.4lf %10.4lf   %12.3e\n",i4,gausposi+i4,gaussians->c_sh[gausposi+i4],\
      gaussians->expo_sh[gausposi+i4],psiexp);
     }
    for (s = 0; s < job->spin_dim; s++) {
      for (m = 0; m < nstates; m++) {
        for (i3 = 0; i3 < sheli; i3++) {
         //psiamp[s * nbands + m] += eigvec->a[s * nbands + m][bfposi + i3] * psiexp * prefac[i3];
         //fprintf(file.out,"psiamp %3d %3d %3d %10.4lf %10.4lf  %10.4lf\n",s,m,i3,eigvec->a[s * nbands + m][bfposi + i3], \
         psiexp, prefac[i3]);
         psiamp[s * nstates + m] += eigvec[s * nbands * dim1 + m * dim1 + bfposi + i3] * psiexp * prefac[i3];
         //fprintf(file.out,"psiamp %3d %3d %3d %3d %3d %3d %10.4lf %10.4lf  %10.4lf %10.4lf\n",ip,bfposi,index_i,s,m,i3, \
         eigvec[s * nbands * dim1 + m * dim1 + bfposi + i3], psiexp, prefac[i3],psiamp[s * nstates + m]);
        } // close loop over i3
       } // close loop on m
      } // close loop on s
     gausposi += shells->ng_sh[index_i];
     bfposi += shells->type_sh[index_i];
    } // close loop on index_i
   } // close loop on ip

    if (job->taskid == 0 && job->verbosity > 1) {
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < nstates; m++) {
          fprintf(file.out,"%3d %3d %10.4lf\n",s,m,psiamp[s * nstates + m]);
         }
        }
       }

}

void write_isosurface_xsf(char *buf, double *grid_data, int grid_size, int grid_par[3], VECTOR_DOUBLE *points, int num_states, int num_points, CRYSTAL *crystal, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i, j, k, j1, j2;
int rem, num;
double *p_grid_data;
FILE *file_prt;

  // ******************************************************************************************
  // * Open file for output                                                                   *
  // ******************************************************************************************

  if (job->taskid == 0) {
  file_prt = fopen(buf, "w");
  if (file_prt == NULL) { fprintf(file.out, "CANNOT OPEN FILE FOR PLOTTING IN write_isosurface_xsf\n"); MPI_Finalize(); exit(1); }
 }

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

   fprintf(file_prt,"#                                     E X C I T O N\n");
   fprintf(file_prt,"#  \n");
   fprintf(file_prt,"#                            Isosurface plotting using XCrysDen\n");
   fprintf(file_prt,"#  \n");

   switch (crystal->type[0]) {

        case 'C':

        fprintf(file_prt,"CRYSTAL\n");                         //Structure begin protocol
        fprintf(file_prt,"PRIMVEC\n");
          for (i = 0; i < 3; i++)
             {
              fprintf(file_prt,"%8.4f%8.4f%8.4f\n", crystal->primitive_cell[i].comp1 * bohr_to_AA, \
              crystal->primitive_cell[i].comp2 * bohr_to_AA, crystal->primitive_cell[i].comp3 * bohr_to_AA);
             }
              fprintf(file_prt,"PRIMCOORD\n");
              fprintf(file_prt,"%4d     1\n",atoms->number_of_atoms_in_unit_cell );

             break;

        case 'S':

        fprintf(file_prt,"SLAB\n");
        fprintf(file_prt,"PRIMVEC\n");
          for (i = 0; i < 2; i++)
            {
             fprintf(file_prt,"%8.4f%8.4f%8.4f\n", crystal->primitive_cell[i].comp1 * bohr_to_AA,\
             crystal->primitive_cell[i].comp2 * bohr_to_AA, crystal->primitive_cell[i].comp3 * bohr_to_AA);
            }
             fprintf(file_prt,"%8.4f%8.4f%8.4f\n", k_zero, k_zero, k_one);
             fprintf(file_prt,"PRIMCOORD\n");
             fprintf(file_prt,"%4d     1\n",atoms->number_of_atoms_in_unit_cell );

            break;

        case 'P':

        fprintf(file_prt,"ATOMS\n");

         break;

        case 'M':

         fprintf(file_prt,"ATOMS\n");

          break;

        }

       for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)
          {
        fprintf(file_prt,"%4d %16.8lf %16.8lf %16.8lf\n",atoms->atomic_number[i], \
        atoms->cell_vector[i].comp1 * bohr_to_AA, atoms->cell_vector[i].comp2 * bohr_to_AA, atoms->cell_vector[i].comp3 * bohr_to_AA);
       }

        fprintf(file_prt,"BEGIN_BLOCK_DATAGRID_3D\n3D_isosurface\n");                       // Begin WF data protocol

    rem = grid_size / 6;
    num = grid_size - rem * 6;
    p_grid_data = grid_data;
    for (j = 0; j < job->spin_dim; j++) {
      for (j1 = 0; j1 < num_states; j1++) {
        for (j2 = 0; j2 < num_points; j2++) {
      //fprintf(file_prt, "# SPIN %d BANDS %d to %d\n",j, bands[0], bands[1]);
      fprintf(file_prt, "BEGIN_DATAGRID_3D_BAND\n");
      fprintf(file_prt, "%8d%5d%5d\n", grid_par[0], grid_par[1], grid_par[2]);
      //fprintf(file_prt, "%8d%5d%5d\n", grid_par[2], grid_par[1], grid_par[0]);
      fprintf(file_prt, "%16.8lf%16.8lf%16.8lf\n", points[0].comp1 * bohr_to_AA, points[0].comp2 * bohr_to_AA, \
      points[0].comp3 * bohr_to_AA);
      for (i = 1; i <= 3 ; i++) {
        fprintf(file_prt, "%16.8lf%16.8lf%16.8lf\n",  points[i].comp1 * bohr_to_AA, \
        points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
       }
        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_prt, " %12.5e", *p_grid_data);
            p_grid_data++;
          }
          fprintf(file_prt, "\n");
        }
        for (i = 0; i < num; i++) {
            fprintf(file_prt, " %12.5e", *p_grid_data);
            p_grid_data++;
        }
        if (num > 0)
          fprintf(file_prt, "\n");
          fprintf(file_prt, "END_DATAGRID_3D\n");   
        } // close loop on j2
        } // close loop on j1
        } // close loop on j

  fprintf(file_prt, "END_BLOCK_DATAGRID_3D\n");  //End protocol print
  fclose(file_prt);

}

