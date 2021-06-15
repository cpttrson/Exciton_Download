#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
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
#include "SYMMETRY.h"
#include "TOOLS.h"
#include "PRINT_UTIL.h"
#include "KPOINTS.h"
#include "CRYSTAL09.h"
#include "PARALLEL.h"
#include "ALLOCATE_MEMORY.h"
#include "INTEGRALS1.h"
#include "FOURIER_TRANSFORM.h"
//#include "INTEGRALS_TWO_CENTRE.h"
#include "INTEGRALS_2C_CRYSTAL.h"
#include "PAIRS_QUADS.h"
#include "DENSITY_MATRIX.h"
#include "GW_BSE_MOLECULE.h"
#include "TDHF_CRYSTAL.h"
#include "PLOTTING.h"

using namespace std;

int init_Complex_array(int dim, Complex *Array)

{

  int i;
  Complex *p_Array;

  p_Array = Array;
  for (i = 0; i < dim; i++) {
    *p_Array = Complex(0, 0);
    p_Array++;
  }

  return 0;

}

void eigenvec_plot(int *is, int *bands, int nkpoints, KPOINT_TRAN *knet, ComplexMatrix *eigvec, VECTOR_INT *kpoints,
int *grid_par, VECTOR_DOUBLE *points, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, p, q, s;
  int nbands = bands[1] - bands[0] + 1, plot_dim = 0;
  int kplot[nkpoints];

  for (i = 0; i < nkpoints; i++)
  kplot[i] = kpoints[i].comp1 * is[1] * is[2] + kpoints[i].comp2 * is[2] + kpoints[i].comp3;
 
  for (i = 0; i < 3; i++) {
    if (grid_par[i] > 1)
      plot_dim++;
  }

  fprintf(file.out, "WAVEFUNCTION PLOT\n");
  fprintf(file.out, "K POINTS IN OBLIQUE COORDINATES\n");
  for (i = 0; i < nkpoints; i++) {
    fprintf(file.out, "%d/%d %d/%d %d/%d\n", kpoints[i].comp1, is[0], kpoints[i].comp2, is[1], kpoints[i].comp3, is[2]);
  }
  fprintf(file.out, "BANDS PLOTTED %d TO %d\n", bands[0], bands[1]);
  fprintf(file.out, "GRID PARAMETERS FOR %d-D PLOT", plot_dim);
  for (i = 0; i < plot_dim; i++)
    fprintf(file.out, " %d ", grid_par[i]);
  fprintf(file.out, "\n");

  if (nkpoints > 10) {
    fprintf(file.out, "NUMBER OF K POINTS MUST BE 10 OR FEWER\n");
    exit(1);
  }

  struct filename {
    char name[15];
  };
  struct filename *filearray;
  filearray = (struct filename *) malloc(10 * sizeof(struct filename));
  if (filearray == NULL) {
    fprintf(file.out, "CANNOT OPEN MEMORY FOR FILEARRAY IN EIGENVEC_PLOT\n");
    exit(1);
  }

  strcpy(filearray[0].name, "eeigvec01.dat");
  strcpy(filearray[1].name, "eeigvec02.dat");
  strcpy(filearray[2].name, "eeigvec03.dat");
  strcpy(filearray[3].name, "eeigvec04.dat");
  strcpy(filearray[4].name, "eeigvec05.dat");
  strcpy(filearray[5].name, "eeigvec06.dat");
  strcpy(filearray[6].name, "eeigvec07.dat");
  strcpy(filearray[7].name, "eeigvec08.dat");
  strcpy(filearray[8].name, "eeigvec09.dat");
  strcpy(filearray[9].name, "eeigvec10.dat");

  FILE *file_prt[10];
  for (k = 0; k < nkpoints; k++) {
    file_prt[k] = fopen(filearray[k].name, "w");
    if (file_prt[k] == NULL) {
      fprintf(file.out, "CANNOT OPEN FILES FOR PLOTTING IN eigenvec_plot\n");
      exit(1);
    }
  }

  int number_Rvec, radmx = 3;
  VECTOR_DOUBLE *Rvec;

  Rvec_grid_size(&number_Rvec, radmx, crystal, file);

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) {
    fprintf(file.out, "There is not enough memory for Rvec array\n");
    exit(1);
  }

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  int PSIdim;
  Complex *Psirkn, *p_Psirkn;

  Psirkn = (Complex *) malloc(nbands * nkpoints * sizeof(Complex));
  if (Psirkn == NULL) {
    fprintf(file.out, "There is not enough memory for Psirkn array");
    exit(1);
  }

  PSIdim = nbands * nkpoints;

  // Generate plot                                                                          *

  int i3, k3;
  int ngrid_points, grid_size;
  VECTOR_DOUBLE *grid, *p_grid;

  grid_size = (grid_par[0] + 1) * (grid_par[1] + 1) * (grid_par[2] + 1);

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) {
    fprintf(file.out, "Cannot open memory for grid array\n");
    exit(1);
  }

  calc_grid(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  p_grid = grid;

  for (i = 0; i < grid_par[0]; i++) {
    for (j = 0; j < grid_par[1]; j++) {
      for (k = 0; k < grid_par[2]; k++) {

        init_Complex_array(PSIdim, Psirkn);
        plot_points(nkpoints, kplot, knet, bands, p_grid, number_Rvec, Rvec, Psirkn, eigvec, atoms, shells, gaussians, job, file);

        p_Psirkn = Psirkn;
        for (k3 = 0; k3 < nkpoints; k3++) {
          fprintf(file_prt[k3],"%10.5lf %10.5lf %10.5lf",p_grid->comp1*bohr_to_AA, p_grid->comp2 * bohr_to_AA, p_grid->comp3 * bohr_to_AA);
          for (i3 = bands[0] - 1; i3 < bands[1]; i3++) {
            fprintf(file_prt[k3], "%15.5lf %15.5lf", p_Psirkn->real(), p_Psirkn->imag());
            p_Psirkn++;
          }
          fprintf(file_prt[k3], "\n");
        }
        p_grid++;
      }
    } // close loops over i, j
    for (k3 = 0; k3 < nkpoints; k3++) {
      fprintf(file_prt[k3], "\n");
    }
  } // close loop on k}
  for (k3 = 0; k3 < nkpoints; k3++)
    fclose(file_prt[k3]);

  free(Psirkn);
  free(filearray);
  free(Rvec);
  free(grid);

  return;

}

void eigenvec_isosurface(int *is, int *bands,  int nkpoints, KPOINT_TRAN *knet, VECTOR_INT *kplot, int *grid_par, \
VECTOR_DOUBLE *points, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, \
SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurfaces of wave functions at selected k points                  *
  // ******************************************************************************************

int i, j, k, p, q, s;
int nbands = bands[1] - bands[0] + 1, plot_dim = 0;
int kpoints[nkpoints];
int number_Rvec, radmx = 2;
int PSIdim;
int i3, j3, k3;
int ngrid_points, grid_size;
int rem, num, count;
int dim0 = job->spin_dim * nbands;
int dim1 = job->spin_dim * nbands * nkpoints;
int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
int block_size = vector_size * sizeof(Complex);
char buf[110];
char xx[10] = "/datafile";
struct filename { char name[15]; } ;
struct filename *filearray;
VECTOR_DOUBLE *Rvec, *grid, *p_grid;
Complex *wavefn, *p_wavefn;
Complex *Psirkn, *p_Psirkn;
FILE *file_prt[10];
ComplexMatrix *eigvec0, *eigvec1;
strcpy(buf,file.directory1);
strcat(buf,xx);
MPI_File fh ;

  // ******************************************************************************************
  // * Open eigenvector and eigenvalue MPI files and read eigenvalues                         *
  // ******************************************************************************************

  AllocateComplexMatrix(&eigvec0,&dim0,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateComplexMatrix(&eigvec1,&dim1,&atoms->number_of_sh_bfns_in_unit_cell,job);
  MPI_File_open(MPI_COMM_WORLD,buf,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  for (s = 0; s < job->spin_dim; s++) {
    for (k = 0; k < nkpoints; k++) {
      p = kplot[k].comp1 * is[1] * is[2] + kplot[k].comp2 * is[2] + kplot[k].comp3;
      q = knet->fbz[p];
      printf("spin %3d k %3d k_unique %3d %10.4lf %10.4lf %10.4lf\n",s,p,q,knet->cart[p].comp1,knet->cart[p].comp2,\
      knet->cart[p].comp3);
      MPI_File_seek(fh, (s * knet->unique + q) * block_size, MPI_SEEK_SET) ;
      MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
      rotate_psi(&eigvec0->a[0][0],&eigvec1->a[(s * nkpoints + k)*nbands][0],nbands,p,knet,atom_p,atoms,R,shells,symmetry,job,file);
     }
    }
    MPI_File_close(&fh);

/*

 // wannier orbital plots

FILE *scf_evectors;
int dummy;
char dummy_string[150], dummy_string1[150];

     scf_evectors = fopen("boys_adatom3", "r");
     read_line(scf_evectors, dummy_string, 124);
     read_line(scf_evectors, dummy_string, 124);

   
     for (j = 0; j < 11; j++) {
     for (i = 0; i < 3; i++) {
     read_line(scf_evectors, dummy_string, 124);
     sscanf(dummy_string,"%s",dummy_string1);
     printf("%s\n",dummy_string);
    }
   for (i = 0; i < atoms->number_of_sh_bfns_in_unit_cell; i++) {
   fscanf(scf_evectors, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&dummy,&eigvec1->a[j * 10 + 0][i].real(),&eigvec1->a[j * 10 + 1][i].real(),&eigvec1->a[j * 10 + 2][i].real(), \
        &eigvec1->a[j * 10 + 3][i].real(),&eigvec1->a[j * 10 + 4][i].real(),&eigvec1->a[j * 10 + 5][i].real(),&eigvec1->a[j * 10 + 6][i].real(),&eigvec1->a[j * 10 + 7][i].real(), \
        &eigvec1->a[j * 10 + 8][i].real(),&eigvec1->a[j * 10 + 9][i].real());
       }
      }

     for (i = 0; i < atoms->number_of_sh_bfns_in_unit_cell; i++) {
     for (j = 0; j < 10; j++) {
     printf("%12.3e ",eigvec1->a[j][i].real());
     eigvec1->a[j][i].imag() = k_zero;
    }
     printf("\n");
   }
*/

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"BANDS FOR PLOTTING\n");
  for (s = 0; s < job->spin_dim; s++) {
    for (k = 0; k < nkpoints; k++) {
      for (i = 0; i < nbands; i++)     {
        for (j = 0; j < atoms->number_of_sh_bfns_in_unit_cell; j++) {
          fprintf(file.out,"%3d %3d %3d   %10.4lf   %10.4lf\n",k,i,j, eigvec1->a[(s * nkpoints + k) * nbands + i][j].real(), \
          eigvec1->a[(s * nkpoints + k) * nbands + i][j].imag());
         }
        }
       }
      }
     }

  // ******************************************************************************************
  // * Open files for output                                                                  *
  // ******************************************************************************************

  filearray = (struct filename *) malloc(10 * sizeof(struct filename));
  if (filearray == NULL) { fprintf(file.out, "CANNOT OPEN MEMORY FOR FILEARRAY IN EIGENVEC_PLOT\n"); exit(1); }

  strcpy(filearray[0].name, "eeigvec01.xsf");
  strcpy(filearray[1].name, "eeigvec02.xsf");
  strcpy(filearray[2].name, "eeigvec03.xsf");
  strcpy(filearray[3].name, "eeigvec04.xsf");
  strcpy(filearray[4].name, "eeigvec05.xsf");
  strcpy(filearray[5].name, "eeigvec06.xsf");
  strcpy(filearray[6].name, "eeigvec07.xsf");
  strcpy(filearray[7].name, "eeigvec08.xsf");
  strcpy(filearray[8].name, "eeigvec09.xsf");
  strcpy(filearray[9].name, "eeigvec10.xsf");

  for (k = 0; k < nkpoints; k++) {
  file_prt[k] = fopen(filearray[k].name, "w");
  if (file_prt[k] == NULL) { fprintf(file.out, "CANNOT OPEN FILES FOR PLOTTING IN eigenvec_plot\n"); exit(1); }
 }

  for (i = 0; i < nkpoints; i++)
  kpoints[i] = kplot[i].comp1 * is[1] * is[2] + kplot[i].comp2 * is[2] + kplot[i].comp3;

  // ******************************************************************************************
  // * Allocate memory                                                                        *
  // ******************************************************************************************

  PSIdim = job->spin_dim * nkpoints * nbands;
  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);

  Psirkn = (Complex *) malloc(PSIdim * sizeof(Complex));
  if (Psirkn == NULL) { fprintf(file.out, "There is not enough memory for Psirkn array"); exit(1); }

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  wavefn = (Complex*) malloc(PSIdim * grid_size * sizeof(Complex));
  if (wavefn == NULL) { fprintf(file.out, "Cannot open memory for wavefn array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  calc_grid(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Rvec grid in eigenvec_isosurface\n");
  for(i = 0; i < number_Rvec; i++)
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);
 }

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

 for (k3 = 0; k3 < nkpoints; k3++) {
 
   fprintf(file_prt[k3],"#                                     E X C I T O N\n");
   fprintf(file_prt[k3],"#  \n");
   fprintf(file_prt[k3],"#                    Wavefunction isosurface plotting using XCrysDen\n");
   fprintf(file_prt[k3],"#  \n");

   switch (crystal->type[0]) {

        case 'C':

        fprintf(file_prt[k3],"CRYSTAL\n");                         //Structure begin protocol
        fprintf(file_prt[k3],"PRIMVEC\n");
          for (i = 0;i < 3; i++)
             {
              fprintf(file_prt[k3],"%8.4f%8.4f%8.4f\n", crystal->primitive_cell[i].comp1 * bohr_to_AA, \
              crystal->primitive_cell[i].comp2 * bohr_to_AA, crystal->primitive_cell[i].comp3 * bohr_to_AA);
             }
              fprintf(file_prt[k3],"PRIMCOORD\n");
              fprintf(file_prt[k3],"%4d     1\n",atoms->number_of_atoms_in_unit_cell );

             break;

        case 'S':

        fprintf(file_prt[k3],"SLAB\n");
        fprintf(file_prt[k3],"PRIMVEC\n");
          for (i = 0; i < 2; i++)
            {
             fprintf(file_prt[k3],"%8.4f%8.4f%8.4f\n", crystal->primitive_cell[i].comp1 * bohr_to_AA,\
             crystal->primitive_cell[i].comp2 * bohr_to_AA, crystal->primitive_cell[i].comp3 * bohr_to_AA);
            }
             fprintf(file_prt[k3],"%8.4f%8.4f%8.4f\n", k_zero, k_zero, k_one);
             fprintf(file_prt[k3],"PRIMCOORD\n");
             fprintf(file_prt[k3],"%4d     1\n",atoms->number_of_atoms_in_unit_cell );

            break;

        case 'P':

        fprintf(file_prt[k3],"ATOMS\n");

         break;

        case 'M':

         fprintf(file_prt[k3],"ATOMS\n");

          break;

        }

       for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)
          {
        fprintf(file_prt[k3],"%4d %16.8lf %16.8lf %16.8lf\n",atoms->atomic_number[i], \
        atoms->cell_vector[i].comp1 * bohr_to_AA, atoms->cell_vector[i].comp2 * bohr_to_AA, atoms->cell_vector[i].comp3 * bohr_to_AA);
       }

        fprintf(file_prt[k3],"BEGIN_BLOCK_DATAGRID_3D\n3D_isosurface\n");                       // Begin WF data protocol

          } // close loop on k3

  count = 0;
  p_grid = grid;

  for (i = 0; i < grid_par[0]; i++) {
    for (j = 0; j < grid_par[1]; j++) {
      for (k = 0; k < grid_par[2]; k++) {
        init_Complex_array(PSIdim, Psirkn);
        plot_points(nkpoints, kpoints, knet, bands, p_grid, number_Rvec, Rvec, Psirkn, eigvec1, atoms, shells, gaussians, job, file);
        p_Psirkn = Psirkn;
        p_wavefn = wavefn + count;
        for (k3 = 0; k3 < nkpoints; k3++) {
          for (j3 = 0; j3 < job->spin_dim; j3++) {
            for (i3 = 0; i3 < nbands; i3++) {
            *p_wavefn = *p_Psirkn;
            p_wavefn += grid_size;
            p_Psirkn++;
          }
         }
        }
        count++;
        p_grid++;
      }
    }
  }

  rem = grid_size / 6;
  num = grid_size - rem * 6;

  for (k3 = 0; k3 < nkpoints; k3++) {
    for (j3 = 0; j3 < job->spin_dim; j3++) {
      for (i3 = 0; i3 < nbands; i3++) {
        fprintf(file_prt[k3], "# KPOINT %d SPIN %d BAND %d\n",k3, j3, bands[0] + i3);
        fprintf(file_prt[k3], "BEGIN_DATAGRID_3D_BAND\n");
        fprintf(file_prt[k3], "%8d%5d%5d\n", grid_par[2], grid_par[1], grid_par[0]);
        fprintf(file_prt[k3], "%16.8lf%16.8lf%16.8lf\n", points[0].comp1 * bohr_to_AA, points[0].comp2 * bohr_to_AA, \
        points[0].comp3 * bohr_to_AA);
        for (i = 3; i >=1 ; i--) {
          fprintf(file_prt[k3], "%16.8lf%16.8lf%16.8lf\n",  points[i].comp1 * bohr_to_AA, \
          points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
         }
        p_wavefn = wavefn + k3 * job->spin_dim * nbands * grid_size + j3 * nbands * grid_size + i3 * grid_size;
        ////p_wavefn = wavefn + k3 * nbands * grid_size + i3 * grid_size;
        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_prt[k3], " %12.5e", p_wavefn->real() * p_wavefn->real() + p_wavefn->imag() * p_wavefn->imag());
            //fprintf(file_prt[k3], " %12.5e", p_wavefn->real());
            p_wavefn++;
          }
          fprintf(file_prt[k3], "\n");
        }
        for (i = 0; i < num; i++) {
          //fprintf(file_prt[k3], " %12.5e", p_wavefn->real());
          fprintf(file_prt[k3], " %12.5e", p_wavefn->real() * p_wavefn->real() + p_wavefn->imag() * p_wavefn->imag());
          p_wavefn++;
        }
        if (num > 0)
          fprintf(file_prt[k3], "\n");
          fprintf(file_prt[k3], "END_DATAGRID_3D\n");   
         }
        }
       } // close loop on i3, j3, k3

  for (k3 = 0; k3 < nkpoints; k3++) {
  fprintf(file_prt[k3], "END_BLOCK_DATAGRID_3D\n");  //End protocol print
  fclose(file_prt[k3]);
 }

  free(Psirkn);
  free(filearray);
  free(Rvec);
  free(grid);
  free(wavefn);
  DestroyComplexMatrix(&eigvec0,job);
  DestroyComplexMatrix(&eigvec1,job);

  return;

}

void density_isosurface(double *grid_sum, int *is, int *bands,  int *grid_par, VECTOR_DOUBLE *points, KPOINT_TRAN *knet, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurface of charge density for selected range of bands            *
  // ******************************************************************************************

int i, j, k, p, q, s;
int nbands = bands[1] - bands[0] + 1, plot_dim = 0;
int number_Rvec, radmx = 3;
int i3, j3, k3;
int ngrid_points, grid_size, spin_grid_size;
int rem, num, count;
int dimp, dimf, dimp_spin, dimf_spin;
int total_tasks;
int begin_p[job->numtasks], end_p[job->numtasks];
double *P, *F;
double *rho_grid, *p_rhogrid, *rho_grid_sum;
//char buf[110];
//char xx[10] = "/datafile";
char filename[15] = "/datafile";
FILE *file_prt;
VECTOR_DOUBLE *Rvec, *grid, *p_grid;
FERMI fermi;
PAIR_TRAN pair_p;
//strcpy(buf,file.directory1);
//strcat(buf,xx);
//MPI_File fh ;

  // ******************************************************************************************
  // * Open file for output                                                                   *
  // ******************************************************************************************

  if (job->taskid == 0) {
  file_prt = fopen("density_isosurface.xsf", "w");
  if (file_prt == NULL) { fprintf(file.out, "CANNOT OPEN FILE FOR PLOTTING IN density_isosurface\n"); exit(1); }
 }

  // ******************************************************************************************
  // * Allocate memory                                                                        *
  // ******************************************************************************************

  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);
  spin_grid_size = job->spin_dim * grid_size;

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  // ******************************************************************************************
  // * Generate pairs and Overlap Matrix                                                      *
  // ******************************************************************************************

  fermi.is[0] = is[0];
  fermi.is[1] = is[1];
  fermi.is[2] = is[2];
  fermi.bands[0] = bands[0];
  fermi.bands[1] = bands[1];
  if (job->spin_polarisation == 1) {
  fermi.bands[2] = bands[2];
  fermi.bands[3] = bands[3];
 }
  //printf("bands %3d %3d %3d %3d\n",bands[0],bands[1],bands[2],bands[3]);
  fermi.nkunique = knet->unique;
  fermi.nktot = knet->nktot;

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //print_pairs(&pair_p,atoms,R,file);
  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file);
  dimp_spin = job->spin_dim * job->dimp;
  dimf_spin = job->spin_dim * job->dimf;

  AllocateDoubleArray(&P,&dimp_spin,job);
  AllocateDoubleArray(&F,&dimf_spin,job);
  AllocateDoubleArray(&rho_grid,&spin_grid_size,job);
  if (job->taskid == 0)
  AllocateDoubleArray(&rho_grid_sum,&spin_grid_size,job);
  job->guess_type = 0;
  job->density = 2;
  job->fix_occ = 1;
  density_matrix_crystal2(&fermi,P,F,knet,filename,&fermi.nkunique,R,R_tables,&pair_p,atom_p,atoms,shells,symmetry,job,file);

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  calc_grid(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Rvec grid in density_isosurface\n");
  for(i = 0; i < number_Rvec; i++)
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);
 }

  total_tasks = grid_size;
  mpi_begin_end(begin_p,end_p,total_tasks,job->numtasks,job,file);

  ResetDoubleArray(rho_grid,&spin_grid_size);
  if (job->taskid == 0)
  ResetDoubleArray(rho_grid_sum,&spin_grid_size);
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    p_grid = grid + p;
    p_rhogrid = rho_grid + p;
    plot_points1(p_rhogrid, p_grid, grid_size, Rvec, number_Rvec, F, &pair_p, R, atoms, shells, gaussians, job, file);
    if (job->taskid == 0) printf("proc %3d plot_points1 %3d %12.3e %12.3e\n",job->taskid,p,*p_rhogrid,*(p_rhogrid + grid_size));
   }

    MPI_Reduce(rho_grid, rho_grid_sum, spin_grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   if (job->taskid == 0) {

  // ******************************************************************************************
  // * Calculate integrated charge density                                                    *
  // ******************************************************************************************

   count = 0;
   for (i = 0; i < job->spin_dim; i++) {
     for (j = 0; j < grid_size; j++) {
       grid_sum[i] += rho_grid[j + count];
       count++;
      }
     }

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

   fprintf(file_prt,"#                                     E X C I T O N\n");
   fprintf(file_prt,"#  \n");
   fprintf(file_prt,"#                       Density isosurface plotting using XCrysDen\n");
   fprintf(file_prt,"#  \n");

   switch (crystal->type[0]) {

        case 'C':

        fprintf(file_prt,"CRYSTAL\n");                         //Structure begin protocol
        fprintf(file_prt,"PRIMVEC\n");
          for (i = 0;i < 3; i++)
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
    p_rhogrid = rho_grid_sum;
    for (j3 = 0; j3 < job->spin_dim; j3++) {
      fprintf(file_prt, "# SPIN %d BANDS %d to %d\n",j3, bands[0], bands[1]);
      fprintf(file_prt, "BEGIN_DATAGRID_3D_BAND\n");
      fprintf(file_prt, "%8d%5d%5d\n", grid_par[2], grid_par[1], grid_par[0]);
      fprintf(file_prt, "%16.8lf%16.8lf%16.8lf\n", points[0].comp1 * bohr_to_AA, points[0].comp2 * bohr_to_AA, \
      points[0].comp3 * bohr_to_AA);
      for (i = 3; i >=1 ; i--) {
        fprintf(file_prt, "%16.8lf%16.8lf%16.8lf\n",  points[i].comp1 * bohr_to_AA, \
        points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
       }
        p_rhogrid = rho_grid_sum;
        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_prt, " %12.5e", *p_rhogrid);
            p_rhogrid++;
          }
          fprintf(file_prt, "\n");
        }
        for (i = 0; i < num; i++) {
            fprintf(file_prt, " %12.5e", *p_rhogrid);
            p_rhogrid++;
        }
        if (num > 0)
          fprintf(file_prt, "\n");
          fprintf(file_prt, "END_DATAGRID_3D\n");   
        } // close loop on j3

  fprintf(file_prt, "END_BLOCK_DATAGRID_3D\n");  //End protocol print
  fclose(file_prt);

 } // close if (job->taskid == 0)

  free(Rvec);
  free(grid);
  free_PAIR_TRAN(&pair_p,job);
  DestroyDoubleArray(&rho_grid,&spin_grid_size,job);
  if (job->taskid == 0)
  DestroyDoubleArray(&rho_grid_sum,&spin_grid_size,job);
  DestroyDoubleArray(&P,&dimp_spin,job);
  DestroyDoubleArray(&F,&dimf_spin,job);

}

void density_isosurface1(double *grid_sum, int *is, int *bands,  int *grid_par, VECTOR_DOUBLE *points, KPOINT_TRAN *knet, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurface of charge density for selected range of bands            *
  // ******************************************************************************************

int i, j, k, p, q, s;
int nbands = bands[1] - bands[0] + 1, plot_dim = 0;
int number_Rvec, radmx = 1;
int i3, j3, k3;
int ngrid_points, grid_size, spin_grid_size;
int rem, num, count;
int dimp, dimf, dimp_spin, dimf_spin;
int total_tasks;
int begin_p[job->numtasks], end_p[job->numtasks];
double *P, *F;
double *rho_grid, *p_rhogrid, *rho_grid_sum;
//char buf[110];
//char xx[10] = "/datafile";
char filename[15] = "/datafile";
FILE *file_prt;
VECTOR_DOUBLE *Rvec, *grid, *p_grid;
FERMI fermi;
PAIR_TRAN pair_p;
//strcpy(buf,file.directory1);
//strcat(buf,xx);
//MPI_File fh ;

  // ******************************************************************************************
  // * Open file for output                                                                   *
  // ******************************************************************************************

  if (job->taskid == 0) {
  file_prt = fopen("density_isosurface.xsf", "w");
  if (file_prt == NULL) { fprintf(file.out, "CANNOT OPEN FILE FOR PLOTTING IN density_isosurface\n"); exit(1); }
 }

  // ******************************************************************************************
  // * Allocate memory                                                                        *
  // ******************************************************************************************

  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);
  spin_grid_size = job->spin_dim * grid_size;

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  // ******************************************************************************************
  // * Generate pairs and Overlap Matrix                                                      *
  // ******************************************************************************************

  fermi.is[0] = is[0];
  fermi.is[1] = is[1];
  fermi.is[2] = is[2];
  fermi.bands[0] = bands[0];
  fermi.bands[1] = bands[1];
  if (job->spin_polarisation == 1) {
  fermi.bands[2] = bands[2];
  fermi.bands[3] = bands[3];
 }
  //printf("bands %3d %3d %3d %3d\n",bands[0],bands[1],bands[2],bands[3]);
  fermi.nkunique = knet->unique;
  fermi.nktot = knet->nktot;
  //CHANGES2014
  allocate_fermi(&fermi,atoms,job,file);
  for (i = 0; i < job->spin_dim * fermi.nkunique * atoms->number_of_sh_bfns_in_unit_cell; i++) {
  fermi.occupation[i] = k_one;
 }
  //CHANGES2014

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //print_pairs(&pair_p,atoms,R,file);
  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file);
  dimp_spin = job->spin_dim * job->dimp;
  dimf_spin = job->spin_dim * job->dimf;

  AllocateDoubleArray(&P,&dimp_spin,job);
  AllocateDoubleArray(&F,&dimf_spin,job);
  AllocateDoubleArray(&rho_grid,&spin_grid_size,job);
  if (job->taskid == 0)
  AllocateDoubleArray(&rho_grid_sum,&spin_grid_size,job);
  job->guess_type = 0;
  job->density = 2;
  job->fix_occ = 1;
  density_matrix_crystal2(&fermi,P,F,knet,filename,&fermi.nkunique,R,R_tables,&pair_p,atom_p,atoms,shells,symmetry,job,file);

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  calc_grid1(grid_par, grid, grid_size, &ngrid_points, points, job, file);
  //calc_grid(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Rvec grid in density_isosurface\n");
  for(i = 0; i < number_Rvec; i++)
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);
 }

  total_tasks = grid_size;
  mpi_begin_end(begin_p,end_p,total_tasks,job->numtasks,job,file);

  ResetDoubleArray(rho_grid,&spin_grid_size);
  if (job->taskid == 0)
  ResetDoubleArray(rho_grid_sum,&spin_grid_size);
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    p_grid = grid + p;
    p_rhogrid = rho_grid + p;
    plot_points1(p_rhogrid, p_grid, grid_size, Rvec, number_Rvec, F, &pair_p, R, atoms, shells, gaussians, job, file);
    if (job->taskid == 0) printf("proc %3d plot_points1 %3d %12.3e %12.3e\n",job->taskid,p,*p_rhogrid,*(p_rhogrid + grid_size));
   }

    MPI_Reduce(rho_grid, rho_grid_sum, spin_grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   if (job->taskid == 0) {

  // ******************************************************************************************
  // * Calculate integrated charge density                                                    *
  // ******************************************************************************************

   count = 0;
   for (i = 0; i < job->spin_dim; i++) {
     for (j = 0; j < grid_size; j++) {
       grid_sum[i] += rho_grid[j + count];
       count++;
      }
     }

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

   fprintf(file_prt,"#                                     E X C I T O N\n");
   fprintf(file_prt,"#  \n");
   fprintf(file_prt,"#                       Density isosurface plotting using XCrysDen\n");
   fprintf(file_prt,"#  \n");

   switch (crystal->type[0]) {

        case 'C':

        fprintf(file_prt,"CRYSTAL\n");                         //Structure begin protocol
        fprintf(file_prt,"PRIMVEC\n");
          for (i = 0;i < 3; i++)
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
    p_rhogrid = rho_grid_sum;
    for (j3 = 0; j3 < job->spin_dim; j3++) {
      //fprintf(file_prt, "# SPIN %d BANDS %d to %d\n",j3, bands[0], bands[1]);
      fprintf(file_prt, "BEGIN_DATAGRID_3D_BAND\n");
      fprintf(file_prt, "%8d%5d%5d\n", grid_par[0], grid_par[1], grid_par[2]);
      //fprintf(file_prt, "%8d%5d%5d\n", grid_par[2], grid_par[1], grid_par[0]);
      fprintf(file_prt, "%16.8lf%16.8lf%16.8lf\n", points[0].comp1 * bohr_to_AA, points[0].comp2 * bohr_to_AA, \
      points[0].comp3 * bohr_to_AA);
      for (i = 1; i <= 3 ; i++) {
      //for (i = 3; i >=1 ; i--) {
        fprintf(file_prt, "%16.8lf%16.8lf%16.8lf\n",  points[i].comp1 * bohr_to_AA, \
        points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
       }
        ////p_rhogrid = rho_grid_sum;
        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_prt, " %12.5e", *p_rhogrid);
            p_rhogrid++;
          }
          fprintf(file_prt, "\n");
        }
        for (i = 0; i < num; i++) {
            fprintf(file_prt, " %12.5e", *p_rhogrid);
            p_rhogrid++;
        }
        if (num > 0)
          fprintf(file_prt, "\n");
          fprintf(file_prt, "END_DATAGRID_3D\n");   
        } // close loop on j3

  fprintf(file_prt, "END_BLOCK_DATAGRID_3D\n");  //End protocol print
  fclose(file_prt);

 } // close if (job->taskid == 0)

  free(Rvec);
  free(grid);
  free_PAIR_TRAN(&pair_p,job);
  DestroyDoubleArray(&rho_grid,&spin_grid_size,job);
  if (job->taskid == 0)
  DestroyDoubleArray(&rho_grid_sum,&spin_grid_size,job);
  DestroyDoubleArray(&P,&dimp_spin,job);
  DestroyDoubleArray(&F,&dimf_spin,job);

}

/*
void eigenvec_isosurface(int *is, int *bands,  int nkpoints, KPOINT_TRAN *knet, ComplexMatrix *eigvec, VECTOR_INT *kpoints, \
int *grid_par, VECTOR_DOUBLE *points, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, \
FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurfaces of wave functions at selected k points                  *
  // ******************************************************************************************

  int i, j, k, p, q, s;
  int nbands = bands[1] - bands[0] + 1, plot_dim = 0;
  int kplot[nkpoints];
  int number_Rvec, radmx = 2;
  int PSIdim;
  int i3, j3, k3;
  int ngrid_points, grid_size;
  int rem, num, count;
  struct filename {
  char name[15];
  };
  struct filename *filearray;
  VECTOR_DOUBLE *Rvec, *grid, *p_grid, *grid_sum, *p_grid_sum;
  Complex *wavefn, *p_wavefn;
  Complex *Psirkn, *p_Psirkn;
  FILE *file_prt[10];

  // ******************************************************************************************
  // * Open files for output                                                                  *
  // ******************************************************************************************

  filearray = (struct filename *) malloc(10 * sizeof(struct filename));
  if (filearray == NULL) { fprintf(file.out, "CANNOT OPEN MEMORY FOR FILEARRAY IN EIGENVEC_PLOT\n"); exit(1); }

  strcpy(filearray[0].name, "eeigvec01.xsf");
  strcpy(filearray[1].name, "eeigvec02.xsf");
  strcpy(filearray[2].name, "eeigvec03.xsf");
  strcpy(filearray[3].name, "eeigvec04.xsf");
  strcpy(filearray[4].name, "eeigvec05.xsf");
  strcpy(filearray[5].name, "eeigvec06.xsf");
  strcpy(filearray[6].name, "eeigvec07.xsf");
  strcpy(filearray[7].name, "eeigvec08.xsf");
  strcpy(filearray[8].name, "eeigvec09.xsf");
  strcpy(filearray[9].name, "eeigvec10.xsf");

  for (k = 0; k < nkpoints; k++) {
  file_prt[k] = fopen(filearray[k].name, "w");
  if (file_prt[k] == NULL) { fprintf(file.out, "CANNOT OPEN FILES FOR PLOTTING IN eigenvec_plot\n"); exit(1); }
 }

  for (i = 0; i < nkpoints; i++)
  kplot[i] = kpoints[i].comp1 * is[1] * is[2] + kpoints[i].comp2 * is[2] + kpoints[i].comp3;

  //for (i = 0; i < 3; i++) {
  //if (grid_par[i] > 1)
  //plot_dim++;
 //}

  // ******************************************************************************************
  // * Allocate memory                                                                        *
  // ******************************************************************************************

  PSIdim = job->spin_dim * nkpoints * nbands;
  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);

  Psirkn = (Complex *) malloc(PSIdim * sizeof(Complex));
  if (Psirkn == NULL) { fprintf(file.out, "There is not enough memory for Psirkn array"); exit(1); }

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  grid_sum = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid_sum == NULL) { fprintf(file.out, "Cannot open memory for grid_sum array\n"); exit(1); }

  wavefn = (Complex*) malloc(PSIdim * grid_size * sizeof(Complex));
  if (wavefn == NULL) { fprintf(file.out, "Cannot open memory for wavefn array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);

  calc_grid(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Rvec grid in eigenvec_isosurface\n");
  for(i = 0; i < number_Rvec; i++)
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);
 }

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************


 for (k3 = 0; k3 < nkpoints; k3++) {
 
   fprintf(file_prt[k3],"#                                     E X C I T O N\n");
   fprintf(file_prt[k3],"#  \n");
   fprintf(file_prt[k3],"#                    Wavefunction isosurface plotting using XCrysDen\n");
   fprintf(file_prt[k3],"#  \n");

   switch (crystal->type[0]) {

        case 'C':

        fprintf(file_prt[k3],"CRYSTAL\n");                         //Structure begin protocol
        fprintf(file_prt[k3],"PRIMVEC\n");
          for (i = 0;i < 3; i++)
             {
              fprintf(file_prt[k3],"%8.4f%8.4f%8.4f\n", crystal->primitive_cell[i].comp1 * bohr_to_AA, \
              crystal->primitive_cell[i].comp2 * bohr_to_AA, crystal->primitive_cell[i].comp3 * bohr_to_AA);
             }
              fprintf(file_prt[k3],"PRIMCOORD\n");
              fprintf(file_prt[k3],"%4d     1\n",atoms->number_of_atoms_in_unit_cell );

             break;

        case 'S':

        fprintf(file_prt[k3],"SLAB\n");
        fprintf(file_prt[k3],"PRIMVEC\n");
          for (i = 0; i < 2; i++)
            {
             fprintf(file_prt[k3],"%8.4f%8.4f%8.4f\n", crystal->primitive_cell[i].comp1 * bohr_to_AA,\
             crystal->primitive_cell[i].comp2 * bohr_to_AA, crystal->primitive_cell[i].comp3 * bohr_to_AA);
            }
             fprintf(file_prt[k3],"%8.4f%8.4f%8.4f\n", k_zero, k_zero, k_one);
             fprintf(file_prt[k3],"PRIMCOORD\n");
             fprintf(file_prt[k3],"%4d     1\n",atoms->number_of_atoms_in_unit_cell );

            break;

        case 'P':

        fprintf(file_prt[k3],"ATOMS\n");

         break;

        case 'M':

         fprintf(file_prt[k3],"ATOMS\n");

          break;

        }

       for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)
          {
        fprintf(file_prt[k3],"%d%16.8lf%16.8lf%16.8lf\n",atoms->atomic_number[i], \
        atoms->cell_vector[i].comp1 * bohr_to_AA, atoms->cell_vector[i].comp2 * bohr_to_AA, atoms->cell_vector[i].comp3 * bohr_to_AA);
       }

        fprintf(file_prt[k3],"BEGIN_BLOCK_DATAGRID_3D\n3D_isosurface\n");                       // Begin WF data protocol

          } // close loop on k3

  count = 0;
  p_grid = grid;

  for (i = 0; i < grid_par[0]; i++) {
    for (j = 0; j < grid_par[1]; j++) {
      for (k = 0; k < grid_par[2]; k++) {
        init_Complex_array(PSIdim, Psirkn);
        plot_points(nkpoints, kplot, knet, bands, p_grid, number_Rvec, Rvec, Psirkn, eigvec, atoms, shells, gaussians, job, file);
        p_Psirkn = Psirkn;
        p_wavefn = wavefn + count;
        for (k3 = 0; k3 < nkpoints; k3++) {
          for (j3 = 0; j3 < job->spin_dim; j3++) {
            for (i3 = 0; i3 < nbands; i3++) {
            *p_wavefn = *p_Psirkn;
            p_wavefn += grid_size;
            p_Psirkn++;
          }
         }
        }
        count++;
        p_grid++;
      }
    }
  }


  rem = grid_size / 6;
  num = grid_size - rem * 6;

  for (k3 = 0; k3 < nkpoints; k3++) {
    for (j3 = 0; j3 < job->spin_dim; j3++) {
      for (i3 = 0; i3 < nbands; i3++) {
        fprintf(file_prt[k3], "# KPOINT %d SPIN %d BAND %d\n",k3, j3, bands[0] + i3);
        fprintf(file_prt[k3], "BEGIN_DATAGRID_3D_BAND\n");
        fprintf(file_prt[k3], "%8d%5d%5d\n", grid_par[2], grid_par[1], grid_par[0]);
        fprintf(file_prt[k3], "%16.8lf%16.8lf%16.8lf\n", points[0].comp1 * bohr_to_AA, points[0].comp2 * bohr_to_AA, \
        points[0].comp3 * bohr_to_AA);
        for (i = 3; i >=1 ; i--) {
          fprintf(file_prt[k3], "%16.8lf%16.8lf%16.8lf\n",  points[i].comp1 * bohr_to_AA, \
          points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
         }
        p_wavefn = wavefn + k3 * job->spin_dim * nbands * grid_size + j3 * nbands * grid_size + i3 * grid_size;
        ////p_wavefn = wavefn + k3 * nbands * grid_size + i3 * grid_size;
        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_prt[k3], " %12.5e", p_wavefn->real() * p_wavefn->real() + p_wavefn->imag() * p_wavefn->imag());
            //fprintf(file_prt[k3], " %12.5e", p_wavefn->real());
            p_wavefn++;
          }
          fprintf(file_prt[k3], "\n");
        }
        for (i = 0; i < num; i++) {
          //fprintf(file_prt[k3], " %12.5e", p_wavefn->real());
          fprintf(file_prt[k3], " %12.5e", p_wavefn->real() * p_wavefn->real() + p_wavefn->imag() * p_wavefn->imag());
          p_wavefn++;
        }
        if (num > 0)
          fprintf(file_prt[k3], "\n");
          fprintf(file_prt[k3], "END_DATAGRID_3D\n");   
         }
        }
       } // close loop on i3, j3, k3


        fprintf(file_prt[k3], "# KPOINT %d SPIN %d BAND %d\n",k3, j3, bands[0] + i3);
        fprintf(file_prt[k3], "BEGIN_DATAGRID_3D_BAND\n");
        fprintf(file_prt[k3], "%8d%5d%5d\n", grid_par[2], grid_par[1], grid_par[0]);
        fprintf(file_prt[k3], "%16.8lf%16.8lf%16.8lf\n", points[0].comp1 * bohr_to_AA, points[0].comp2 * bohr_to_AA, \
        points[0].comp3 * bohr_to_AA);
        for (i = 3; i >=1 ; i--) {
          fprintf(file_prt[k3], "%16.8lf%16.8lf%16.8lf\n",  points[i].comp1 * bohr_to_AA, \
          points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
         }

  for (k3 = 0; k3 < nkpoints; k3++) {
    for (j3 = 0; j3 < job->spin_dim; j3++) {
      for (i3 = 0; i3 < nbands; i3++) {
        p_wavefn = wavefn + k3 * job->spin_dim * nbands * grid_size + j3 * nbands * grid_size + i3 * grid_size;
        for (k = 0; k < grid_size; k++) {
          grid_sum[k] += p_wavefn->real() * p_wavefn->real() + p_wavefn->imag() * p_wavefn->imag());
            p_wavefn++;
          }
         }
        }
       } // close loop on i3, j3, k3

        for (k = 0; k < rem; k++) {
          for (i = 0; i < 6; i++) {
            fprintf(file_prt[k3], " %12.5e", *p_grid_sum);
            p_grid_sum++;
          }
         }

        for (k = 0; k < num; k++) {
            fprintf(file_prt[k3], " %12.5e", *p_grid_sum);
            p_grid_sum++;
          }
        if (num > 0)
          fprintf(file_prt[k3], "\n");

          fprintf(file_prt[k3], "END_DATAGRID_3D\n");   


  for (k3 = 0; k3 < nkpoints; k3++) {
  fprintf(file_prt[k3], "END_BLOCK_DATAGRID_3D\n");  //End protocol print
  fclose(file_prt[k3]);
 }

  free(Psirkn);
  free(filearray);
  free(Rvec);
  free(grid);
  free(grid_sum);
  free(wavefn);

  return;

}
*/

/*
void band_plot(int *is, int *bands, int nkpoints, double *atom_proj, KPOINT_TRAN *knet, ComplexMatrix *eigvec, double *eigval, VECTOR_INT *kpoints, ATOM_TRAN *atom_p, ATOM_TRAN *atom_i, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, p;
  int nbands = bands[1] - bands[0] + 1, plot_dim = 0;
  int nk[2], kplot[nkpoints];
  int Function[8];
  int dim, dimf, dimg, dimp, num_p;
  int max_f[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  double projection_weight[nkpoints][nbands];
  double weight[nkpoints][nbands];
  double pointsize;
  FILE *band_structure_0, *band_structure_1, *band_gnu;
  ComplexMatrix *S_k;
  PAIR_TRAN pair_p;
  INT_1E one_ints;

  if (job->taskid == 0) {
  band_structure_0 = fopen("band_structure_0.dat", "w");
  if (band_structure_0 == NULL) { fprintf(file.out, "cannot open file band_structure_0.dat\n"); MPI_Finalize(); exit(1); }
 }

  AllocateComplexMatrix(&S_k,&atoms->number_of_sh_bfns_in_unit_cell,&atoms->number_of_sh_bfns_in_unit_cell,job);

  for (i = 0; i < nkpoints; i++) {
  kplot[i] = kpoints[i].comp1 * is[1] * is[2] + kpoints[i].comp2 * is[2] + kpoints[i].comp3;
  //fprintf(file.out,"%3d %3d\n",i,kplot[i]);
 }

  for (i = 0; i < nkpoints; i++) {
    for (j = 0; j < nbands; j++) {
      projection_weight[i][j] = k_zero;
      weight[i][j] = k_zero;
     }
    }

  int num_pairs[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int pair_limits[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int max_pairs;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  print_pairs(&pair_p,atoms,R,file);

  //count_pairs(atoms,atom_p,atom_i,symmetry,R,R_tables,&num_p,num_pairs,&max_pairs,pair_limits,job,file);
  //allocate_density_pairs(&pair_p,num_p,max_pairs,symmetry,atoms,R_tables,file);
  //generate_pairs(atoms,atom_p,atom_i,symmetry,R,R_tables,&pair_p,&num_p,&max_pairs,pair_limits,job,file);
  //print_pairs(&num_p, &pair_p, file);


  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  //mpi_begin_end(begin_p,end_p,num_p,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  //array_dimensions(&dim, &dimg, &num_p, &pair_p, atoms, job, file); // don't change
  //printf("dim %d %d %d\n",dim,dimg,num_p);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e(&one_ints, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //fock_element_1e(&one_ints, dim, &pair_p, num_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);

  //printf("overlap %10.4lf \n",one_ints.Overlap[0]);
  //sh_array_dimensions(&dimp, &dimf, &num_p, &pair_p, atoms, job, file); // don't change
  //printf("dimp %d %d %d\n",num_p,dimp,dimf);
  //printf("overlap %10.4lf \n",one_ints.Overlap[0]);

  MPI_Allgatherv(&one_ints.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);

  //printf("overlap %10.4lf \n",one_ints.Overlap[0]);

  for (i = 0; i < nkpoints; i++) {
    nk[0] = kplot[i];
    nk[1] = kplot[i];
      fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
    //fourier_transform1(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, num_p, R, atoms, shells, symmetry, job, file);
    //fourier_transform2(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, num_p, R, atoms, shells, symmetry, job, file);
    fprintf(file.out,"KPOINT %3d %3d %3d\n",i,kplot[i],knet->fbz[nk[0]]);
    //print_complex_matrix(S_k, file);
    for (j = 0; j < nbands; j++) {
      for (k = 0; k < atoms->number_of_sh_bfns_in_unit_cell; k++) {
        for (l = 0; l < atoms->number_of_sh_bfns_in_unit_cell; l++) {
          projection_weight[i][j] += (conj(eigvec->a[i * nbands + j][k]) * S_k->a[k][l] * eigvec->a[i * nbands + j][l]).real() * atom_proj[k];
//fprintf(file.out,"%3d %3d %3d %3d %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf \n",i,j,k,l,projection_weight[i][j],(conj(eigvec->a[i * nbands + j][k])).real() , (conj(eigvec->a[i * nbands + j][k])).imag(),(S_k->a[k][l]).real() , (S_k->a[k][l]).imag(),(eigvec->a[i * nbands + j][l]).real(),(eigvec->a[i * nbands + j][l]).imag());
          weight[i][j] += (conj(eigvec->a[i * nbands + j][k]) * S_k->a[k][l] * eigvec->a[i * nbands + j][l]).real();
         }
        }
       }
      }


  for (i = 0; i < nbands; i++) {
    for (j = 0; j < nkpoints; j++) {
      fprintf(file.out,"%3d %3d %20.15lf %20.15lf    %10.4lf %10.4lf    %10.4lf %10.4lf\n",i,j,projection_weight[j][i],weight[j][i],eigvec->a[j * nbands + i][0].real(),eigvec->a[j * nbands + i][0].imag(),eigvec->a[j * nbands + i][1].real(),eigvec->a[j * nbands + i][1].imag());
      fprintf(band_structure_0,"%3d %10.4lf %10.4e\n",j, eigval[i + nbands * knet->fbz[kplot[j]]], 6*projection_weight[j][i]);
      //fprintf(band_structure_0,"%3d %10.4lf %10.4e\n",j, (eigval[i + nbands * knet->fbz[kplot[j]]] - 0*job->fermi_energy) * au_to_eV, projection_weight[j][i]);
     }
    fprintf(band_structure_0,"\n");
   }
    fflush(band_structure_0);

    DestroyComplexMatrix(&S_k,job);

}

void band_plot1(int *is, int *bands, int nkpoints, int nproj, KPOINT_TRAN *knet, ComplexMatrix *eigvec, double *eigval, VECTOR_INT *kpoints, ATOM_TRAN *atom_p, ATOM_TRAN *atom_i, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, n, p, s;
  int nbands = bands[1] - bands[0] + 1, plot_dim = 0, nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int nk[2], kplot[nkpoints];
  int Function[8];
  int dim, dimf, dimg, dimp, num_p;
  int dim2;
  int max_f[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  double projection_weight[nproj + 1][nkpoints][nbands];
  double weight[nkpoints][nbands];
  double pointsize;
  char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
  FILE *band_structure_0, *band_structure_1, *band_gnu;
  Complex alpha, beta;
  alpha.real() = k_one;
  alpha.imag() = k_zero;
  beta.real() = k_zero;
  beta.imag() = k_zero;
  ComplexMatrix *S_k, *S, *tmp;
  PAIR_TRAN pair_p;
  INT_1E one_ints;

  if (job->taskid == 0) {
  band_structure_0 = fopen("band_structure_0.dat", "w");
  if (band_structure_0 == NULL) { fprintf(file.out, "cannot open file band_structure_0.dat\n"); MPI_Finalize(); exit(1); }
 }

  dim2 = (nproj + 1) * job->spin_dim * nkpoints * nbands;

  AllocateComplexMatrix(&S_k,&nbfn,&nbfn,job);
  AllocateComplexMatrix(&S,&dim2,&nbands,job);
  AllocateComplexMatrix(&tmp,&nbfn,&nbands,job);
  ResetComplexMatrix(S);
  ResetComplexMatrix(tmp);

  for (i = 0; i < nkpoints; i++) {
  kplot[i] = kpoints[i].comp1 * is[1] * is[2] + kpoints[i].comp2 * is[2] + kpoints[i].comp3;
  //fprintf(file.out,"%3d %3d\n",i,kplot[i]);
 }

  for (n = 0; n <= nproj; n++) {
    for (i = 0; i < nkpoints; i++) {
      for (j = 0; j < nbands; j++) {
        projection_weight[n][i][j] = k_zero;
        weight[i][j] = k_zero;
       }
      }
     }

  int num_pairs[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int pair_limits[atoms->number_of_unique_atoms * atoms->number_of_unique_atoms];
  int max_pairs;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  print_pairs(&pair_p,atoms,R,file);

  //count_pairs(atoms,atom_p,atom_i,symmetry,R,R_tables,&num_p,num_pairs,&max_pairs,pair_limits,job,file);
  //allocate_density_pairs(&pair_p,num_p,max_pairs,symmetry,atoms,R_tables,file);
  //generate_pairs(atoms,atom_p,atom_i,symmetry,R,R_tables,&pair_p,&num_p,&max_pairs,pair_limits,job,file);
  //print_pairs(&num_p, &pair_p, file);


  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  //mpi_begin_end(begin_p,end_p,num_p,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  //array_dimensions(&dim, &dimg, &num_p, &pair_p, atoms, job, file); // don't change
  //printf("dim %d %d %d\n",dim,dimg,num_p);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e(&one_ints, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //fock_element_1e(&one_ints, dim, &pair_p, num_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);

  //printf("overlap %10.4lf \n",one_ints.Overlap[0]);
  //sh_array_dimensions(&dimp, &dimf, &num_p, &pair_p, atoms, job, file); // don't change
  //printf("dimp %d %d %d\n",num_p,dimp,dimf);
  //printf("overlap %10.4lf \n",one_ints.Overlap[0]);

  MPI_Allgatherv(&one_ints.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);

  //printf("overlap %10.4lf \n",one_ints.Overlap[0]);
      ComplexMatrix *eigvec0;
      AllocateComplexMatrix(&eigvec0,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);

      int dim5 = nkpoints * nbands;
      int dim6 = job->spin_dim * dim5;
      for (n = 0; n <= nproj; n++) {
        for (s = 0; s < job->spin_dim; s++) {
          for (i = 0; i < nkpoints; i++) {
            for (k = 0; k < nbands; k++) {
              for (l = 0; l < atoms->number_of_sh_bfns_in_unit_cell; l++) {
                eigvec0->a[k][l] = eigvec->a[n * dim6 + s * dim5 + i * nbands + k][l];
               }
              }
        nk[0] = kplot[i];
        nk[1] = kplot[i];
         fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
         //fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
        //fourier_transform1(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, num_p, R, atoms, shells, symmetry, job, file);
        //fourier_transform2(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, num_p, R, atoms, shells, symmetry, job, file);
        ResetComplexMatrix(S);
        fprintf(file.out,"KPOINT %3d %3d %3d\n",i,kplot[i],knet->fbz[nk[0]]);
        ComplexGEMM1(&NoTrans, &Trans, &alpha, &S_k, &eigvec0, &beta, &tmp);
        for (j = 0; j < nbfn; j++) {
          for (p = 0; p < nbfn;p++) {
            (eigvec0->a[j][p]).imag() *= -k_one;
           }
          }
        ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigvec0, &tmp, &beta, &S);
        for (p = 0; p < nbands; p++) {
          projection_weight[n][i][p] = S->a[p][p].real();
          //fprintf(file.out,"%3d %3d %16.10lf %16.10lf\n",j,p,S->a[p][p].real(),S->a[p][p].imag());
         }
             }
        }
       }

      for (i = 0; i < nbands; i++) {
        for (j = 0; j < nkpoints; j++) {
          //fprintf(band_structure_0,"%3d %10.4lf %10.4e\n",j, eigval[i + nbands * knet->fbz[kplot[j]]], projection_weight[n][j][i]);
          fprintf(band_structure_0,"%3d %10.4lf ",j, eigval[i + nbands * knet->fbz[kplot[j]]]);
          for (n = 0; n <= nproj; n++) {
      fprintf(band_structure_0," %10.4e ",projection_weight[n][j][i]);
     }
    fprintf(band_structure_0,"\n");
     }
    fprintf(band_structure_0,"\n");
   }
    fflush(band_structure_0);

    DestroyComplexMatrix(&S_k,job);
    DestroyComplexMatrix(&tmp,job);


(eigvec0->a[0][0]).real() = 1.0;
(eigvec0->a[0][0]).imag() = 0.1;
(eigvec0->a[0][1]).real() = 2.0;
(eigvec0->a[0][1]).imag() = 0.2;
(eigvec0->a[1][0]).real() = 3.0;
(eigvec0->a[1][0]).imag() = 0.3;
(eigvec0->a[1][1]).real() = 4.0;
(eigvec0->a[1][1]).imag() = 0.4;

(S_k->a[0][0]).real() = 2.0;
(S_k->a[0][0]).imag() = 0.2;
(S_k->a[0][1]).real() = 1.0;
(S_k->a[0][1]).imag() = 0.1;
(S_k->a[1][0]).real() = 4.0;
(S_k->a[1][0]).imag() = 0.4;
(S_k->a[1][1]).real() = 3.0;
(S_k->a[1][1]).imag() = 0.3;
    //fprintf(file.out,"TMP %3d\n",i);
    //print_complex_matrix(tmp, file);
*/
/*
    for (j = 0; j < nkpoints; j++) {
  for (i = 0; i < nbands; i++) {
fprintf(file.out,"BAND %3d %3d %5d\n",j,i,&(eigvec->a[j * nbands]));
      for (k = 0; k < nbfn; k++) {
fprintf(file.out,"%3d %3d %3d %10.4lf %10.4lf \n",j,i,k, eigvec->a[j * nbands + i][k].real(),eigvec->a[j * nbands + i][k].imag()); }}}

}
*/

void band_plot(int *is, int *bands, int nkpoints, int nproj, double atom_proj[][5], KPOINT_TRAN *knet, VECTOR_INT *kpoints, VECTOR_INT *kpoints1, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int b, i, ii, i1, j, k, l, m, n, p, q, s;
  int nbands = bands[1] - bands[0] + 1;
  int nk[2];
  int Function[8];
  int dim, dimf, dimg, dimp;
  int dim0 = job->spin_dim * nbands;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  int dim2 = dim1 * knet->unique;
  int dim5 = nkpoints * nbands;
  int dim6 = job->spin_dim * dim5;
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int info = 0;
  int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
  int block_size = vector_size * sizeof(Complex);
  double kmag;
  double projection_weight[nproj + 1][nkpoints][nbands];
  double *eigval, *eigenvalues;
  char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
  char uplo = 'U';
  char jobz = 'V';
  char buf[110], buf2[110];
  char xx[10] = "/datafile", yy[10] = "/evalfile";
  VECTOR_DOUBLE k_temp;
  FILE *band_structure_0, *band_structure_1, *band_gnu;
  Complex alpha = Complex(k_one, k_zero);
  Complex  beta = Complex(k_zero, k_zero);
  //Complex alpha, beta;
  //alpha.real() = k_one;
  //alpha.imag() = k_zero;
  //beta.real() = k_zero;
  //beta.imag() = k_zero;
  ComplexMatrix *S_k;
  ComplexMatrix *xtmp, *S_x;
  ComplexMatrix *eigvec0, *eigvec1, *eigvec2;
  PAIR_TRAN pair_p;
  INT_1E one_ints, one_ints_buffer;
  MPI_File fh, gh;

  strcpy(buf,file.directory1);
  strcat(buf,xx);
  strcpy(buf2,file.directory1);
  strcat(buf2,yy);

  if (job->taskid == 0) {
  band_structure_0 = fopen("band_structure_0.dat", "w");
  if (band_structure_0 == NULL) { fprintf(file.out, "cannot open file band_structure_0.dat\n"); MPI_Finalize(); exit(1); }
 }

  for (m = 0; m <= nproj; m++) {
    for (k = 0; k < nkpoints; k++) {
      for (b = 0; b < nbands; b++) {
        projection_weight[m][k][b] = k_zero;
       }
      }
     }

  // ******************************************************************************************
  // * Generate pairs and Overlap Matrix                                                      *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //print_pairs(&pair_p,atoms,R,file);

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p, end_p, receive_p, offset_p, &pair_p, atoms, job, file);
  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  //fock_element_1e(&one_ints, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //MPI_Allgatherv(&one_ints.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  fock_element_1e(&one_ints_buffer, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allgatherv(&one_ints_buffer.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);

  // ******************************************************************************************
  // * Open eigenvector and eigenvalue MPI files and read eigenvalues                         *
  // ******************************************************************************************

  AllocateComplexMatrix(&eigvec0,&dim0,&dim1,job);
  AllocateComplexMatrix(&eigvec1,&dim0,&dim1,job);
  AllocateComplexMatrix(&eigvec2,&dim0,&dim1,job);
  AllocateComplexMatrix(&S_k,&dim1,&dim1,job);
  AllocateComplexMatrix(&S_x,&dim0,&dim0,job);
  AllocateComplexMatrix(&xtmp,&dim0,&dim1,job);
  AllocateDoubleArray(&eigval,&dim2,job);

  ResetDoubleArray(eigval,&dim2);
  MPI_File_open(MPI_COMM_WORLD,buf ,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
  MPI_File_read(gh, eigval, job->spin_dim * knet->unique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"EIGENVALUES FOR PLOTTING\n");
  for (k = 0; k < knet->unique; k++) {
    for (i = 0; i < nbands; i++)     {
      fprintf(file.out,"%3d %3d   %10.4lf\n",k,i, eigval[k * nbands + i]);
     }
    }
   }


  ResetComplexMatrix(eigvec0);
  ResetComplexMatrix(eigvec1);
  for (s = 0; s < job->spin_dim; s++) {
    for (k = 0; k < nkpoints; k++) {
      p = kpoints[k].comp1 * is[1] * is[2] + kpoints[k].comp2 * is[2] + kpoints[k].comp3;
      q = knet->fbz[p];
      printf("spin %3d k %3d k_unique %3d %10.4lf %10.4lf %10.4lf\n",s,p,q,knet->cart[p].comp1,knet->cart[p].comp2,knet->cart[p].comp3);
      MPI_File_seek(fh, (s * knet->unique + q) * block_size, MPI_SEEK_SET) ;
      MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
      rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],nbands,p,knet,atom_p,atoms,R,shells,symmetry,job,file);
      //CHANGE2014rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],bands,p,knet,atom_p,atoms,R,shells,symmetry,job,file);
      if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"eigvec1 array %3d \n",k);
      print_complex_matrix(eigvec0, file);
     }
      nk[0] = p;
      nk[1] = p;
      fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

  // ******************************************************************************************
  // * Check that overlap matrix is correctly Fourier transformed                             *
  // ******************************************************************************************

      if (job->taskid == 0 && job->verbosity > 1) {

      fprintf(file.out,"S_k2 orthogonalised %d\n",k);
      ComplexMatrix *xtrn1, *S_k1, *S_k2, *xtmp1, *eigenvectors;
      AllocateComplexMatrix(&S_k1,&dim1,&dim1,job);
      AllocateComplexMatrix(&S_k2,&dim1,&dim1,job);
      fourier_transform(&one_ints.Overlap[0], &S_k1->a[0][0], knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      fourier_transform(&one_ints.Overlap[0], &S_k2->a[0][0], knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      AllocateComplexMatrix(&eigenvectors,&dim1,&dim1,job);
      AllocateDoubleArray(&eigenvalues,&dim1,job);
      ResetDoubleArray(eigenvalues,&dim1);
      DiagonaliseHermitian(&S_k1, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
      //fprintf(file.out,"Eigenvectors of Fourier transformed overlap matrix at k point %3d\n",k);
      //print_complex_matrix(S_k1, file);
      n = 0;
      for (ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
      if (eigenvalues[ii] < 1.0e-04 && job->taskid == 0)  fprintf(file.out,"small eigenvalue %d %d %e\n",k,ii,eigenvalues[ii]);
      if (eigenvalues[ii] < 1.0e-04)  n++;
     }
      int dim11 = dim1 - n;
      AllocateComplexMatrix(&xtrn1,&dim11,&dim1,job);
      AllocateComplexMatrix(&xtmp1,&dim11,&dim1,job);
      for (i1 = n; i1 < dim1; i1++) {
      //fprintf(file.out,"%3d %10.4lf \n",k,*(eigenvalues + k));
      for (j = 0; j < dim1; j++) {
      xtrn1->a[i1 - n][j] = eigenvectors->a[i1][j] / sqrt(*(eigenvalues + i1));
      //fprintf(file.out,"%10.4lf ", (eigenvectors->a[k][j]).real());
     }
      //fprintf(file.out,"\n");
     }
      print_complex_matrix(xtrn1, file);
      ResetComplexMatrix(xtmp1);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &S_k2, &beta, &xtmp1);
      ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &xtrn1, &beta, &S_k2);
      fprintf(file.out,"S_k2 orthogonalised %d\n",k);
      print_complex_matrix(S_k2, file);
      DestroyComplexMatrix(&S_k1,job);
      DestroyComplexMatrix(&S_k2,job);
      DestroyComplexMatrix(&xtrn1,job);
      DestroyComplexMatrix(&xtmp1,job);
      DestroyComplexMatrix(&eigenvectors,job);
      DestroyDoubleArray(&eigenvalues,&dim1,job);

    }

  // ******************************************************************************************
  // * Calculate projection weights                                                           *
  // ******************************************************************************************

      for (m = 0; m <= nproj; m++) {
        for (j = 0; j < nbands; j++)     {
          for (l = 0; l < atoms->number_of_sh_bfns_in_unit_cell; l++) {
            eigvec2->a[j][l] = eigvec1->a[j][l] * atom_proj[l][m];
            //fprintf(file.out,"%3d %3d %3d   %10.4lf   %10.4lf\n",m,j,l,eigvec2->a[j][l].real(),eigvec2->a[j][l].imag());
           }
          }
           ResetComplexMatrix(S_x);
           ResetComplexMatrix(xtmp);
           ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigvec1, &S_k, &beta, &xtmp);
           ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp, &eigvec2, &beta, &S_x);
           //fprintf(file.out,"S_x %d\n",k);
           //print_complex_matrix(S_x, file);
           for (b = 0; b < nbands; b++) {
             projection_weight[m][k][b] = S_x->a[b][b].real();
             fprintf(file.out,"s %3d k %3d proj %3d band %3d %16.10lf %16.10lf %16.4e\n",\
             s,k,m,b,S_x->a[b][b].real(),S_x->a[b][b].imag(),fabs(k_one - S_x->a[b][b].real()));
            } // close loop over b
           } // close loop over m
          } // close loop over k
         } // close loop over s

  // ******************************************************************************************
  // * Write eigenvalues and projection weights                                               *
  // ******************************************************************************************

      for (b = 0; b < nbands; b++) {
        for (k = 0; k < nkpoints; k++) {
          p = kpoints[k].comp1 * is[1] * is[2] + kpoints[k].comp2 * is[2] + kpoints[k].comp3;
          q = knet->fbz[p];
          k_temp.comp1 = kpoints1[k].comp1 / double (is[0]) * crystal->reciprocal_cell[0].comp1 + \
                         kpoints1[k].comp2 / double (is[1]) * crystal->reciprocal_cell[1].comp1 + \
                         kpoints1[k].comp3 / double (is[2]) * crystal->reciprocal_cell[2].comp1;
          k_temp.comp2 = kpoints1[k].comp1 / double (is[0]) * crystal->reciprocal_cell[0].comp2 + \
                         kpoints1[k].comp2 / double (is[1]) * crystal->reciprocal_cell[1].comp2 + \
                         kpoints1[k].comp3 / double (is[2]) * crystal->reciprocal_cell[2].comp2;
          k_temp.comp3 = kpoints1[k].comp1 / double (is[0]) * crystal->reciprocal_cell[0].comp3 + \
                         kpoints1[k].comp2 / double (is[1]) * crystal->reciprocal_cell[1].comp3 + \
                         kpoints1[k].comp3 / double (is[2]) * crystal->reciprocal_cell[2].comp3;
          kmag = sqrt(double_vec_dot(&k_temp,&k_temp)) / bohr_to_AA ;
          fprintf(band_structure_0,"%3d %10.4lf %10.4lf ",k, kmag, eigval[q * nbands + b]);
          for (m = 0; m <= nproj; m++) {
            fprintf(band_structure_0," %11.4e ",projection_weight[m][k][b]);
           }
          fprintf(band_structure_0,"\n");
         }
        fprintf(band_structure_0,"\n");
       }
      fflush(band_structure_0);

    DestroyDoubleArray(&eigval,&dim2,job);
    DestroyComplexMatrix(&S_k,job);
    DestroyComplexMatrix(&S_x,job);
    DestroyComplexMatrix(&xtmp,job);
    DestroyComplexMatrix(&eigvec0,job);
    DestroyComplexMatrix(&eigvec1,job);
    DestroyComplexMatrix(&eigvec2,job);
    MPI_File_close(&fh);
    MPI_File_close(&gh);

}

void bulk_band_projection(int *is, int *is_proj, int *bands, KPOINT_TRAN *knet, double *eigval, ATOM_TRAN *atom_p, ATOM_TRAN *atom_i, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, m, n, p, q, s;
  int count, count1, k_start;
  int nbands = bands[1] - bands[0] + 1;
  int not_orthogonal;
  int G_surf_num = 5;
  double kdotn, temp_k, ae_bd, a, b, c, d, e, f, surface_proj_mag;
  double lambda, mu;
  FILE *band_project_0, *band_project_1;
  VECTOR_INT surface_norm, surface_proj;
  VECTOR_DOUBLE G_vec;

  if (job->taskid == 0) {
  band_project_0 = fopen("band_project_0.dat", "w");
  if (band_project_0 == NULL) { fprintf(file.out, "cannot open file band_project_0.dat\n"); MPI_Finalize(); exit(1); }
 }
 
      surface_norm.comp1 = 1;
      surface_norm.comp2 = 1;
      surface_norm.comp3 = 1;

      surface_proj.comp1 = -1;
      surface_proj.comp2 =  1;
      surface_proj.comp3 =  0;

      a = surface_norm.comp1 * crystal->reciprocal_cell[0].comp1 + surface_norm.comp2 * crystal->reciprocal_cell[1].comp1 + \
          surface_norm.comp3 * crystal->reciprocal_cell[2].comp1;
      b = surface_norm.comp1 * crystal->reciprocal_cell[0].comp2 + surface_norm.comp2 * crystal->reciprocal_cell[1].comp2 + \
          surface_norm.comp3 * crystal->reciprocal_cell[2].comp2;
      c = surface_norm.comp1 * crystal->reciprocal_cell[0].comp3 + surface_norm.comp2 * crystal->reciprocal_cell[1].comp3 + \
          surface_norm.comp3 * crystal->reciprocal_cell[2].comp3;
     
      d = surface_proj.comp1 * crystal->primitive_cell[0].comp1 + surface_proj.comp2 * crystal->primitive_cell[1].comp1 + \
          surface_proj.comp3 * crystal->primitive_cell[2].comp1;
      e = surface_proj.comp1 * crystal->primitive_cell[0].comp2 + surface_proj.comp2 * crystal->primitive_cell[1].comp2 + \
          surface_proj.comp3 * crystal->primitive_cell[2].comp2;
      f = surface_proj.comp1 * crystal->primitive_cell[0].comp3 + surface_proj.comp2 * crystal->primitive_cell[1].comp3 + \
          surface_proj.comp3 * crystal->primitive_cell[2].comp3;

      not_orthogonal = a * d + b * e + c * f;
 
      if (not_orthogonal) {
      if (job->taskid == 0)
      fprintf(file.out,"Surface normal and surface projection are not orthogonal\n");
      MPI_Finalize();
      exit(1);
     }

      ae_bd = k_one / (a * e - b * d);
   
      printf("%10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf   %10.4lf\n",a,b,c,d,e,f,ae_bd);

      a /= ae_bd;
      e /= ae_bd;
      b /= ae_bd;
      d /= ae_bd;

VECTOR_DOUBLE G_surf[G_surf_num];
G_surf[0].comp1 =  0.0;
G_surf[0].comp2 =  0.0;
G_surf[0].comp3 =  0.0;
G_surf[1].comp1 = -0.076555;
G_surf[1].comp2 =  0.535890;
G_surf[1].comp3 =  0.127593;
G_surf[2].comp1 =  0.076555;
G_surf[2].comp2 = -0.535890;
G_surf[2].comp3 = -0.127593;
G_surf[3].comp1 =  0.127593;
G_surf[3].comp2 =  0.127593;
G_surf[3].comp3 = -0.255186;
G_surf[4].comp1 = -0.127593;
G_surf[4].comp2 = -0.127593;
G_surf[4].comp3 =  0.255186;

G_surf_num = 1;

      count = 0;
      for (s = 0; s < job->spin_dim; s++) {
        for (i = 0; i < is[0]; i++) {
          for (j = 0; j < is[1]; j++) {
            for (k = 0; k < is[2]; k++) {
              p = i * is[1] * is[2] + j * is[2] + k;
              q = knet->fbz[p];
fprintf(file.out,"%3d %3d %3d  %3d %3d    %12.6lf %12.6lf %12.6lf    %12.6lf %12.6lf %12.6lf     %10.4lf\n",i,j,k,p,G->last_vector,G->vec_bi[146].comp1,G->vec_bi[146].comp2,G->vec_bi[146].comp3,knet->cart[p].comp1,knet->cart[p].comp2,knet->cart[p].comp3, eigval[nbands * q + 3] * au_to_eV + 3.322);
              for (l = 0; l < nbands; l++) {
                for (m = 0; m < G_surf_num; m++) {
                  for (n = 0; n < G->last_vector; n++) {
                    G_vec.comp1 = G->vec_bi[n].comp1 + knet->cart[p].comp1;
                    G_vec.comp2 = G->vec_bi[n].comp2 + knet->cart[p].comp2;
                    G_vec.comp3 = G->vec_bi[n].comp3 + knet->cart[p].comp3;
                    lambda = e * G_vec.comp1 - d * G_vec.comp2;
                    mu    = -b * G_vec.comp1 + a * G_vec.comp2;
                    kdotn = G_vec.comp2;
                    //kdotn = G_vec.comp1;
                    //kdotn = (- G_vec.comp1 + G_vec.comp2) / sqrt(two);
                    //if (fabs(c * lambda + f * mu - G_vec.comp3) < 0.0001 && fabs(kdotn) < 0.95)
                    //if (fabs(G_vec.comp1 + G_vec.comp2 - two * G_vec.comp3) < 0.0001 && fabs(kdotn) < 0.95)
                   //this oneif (fabs(G_vec.comp1 + G_vec.comp2 - two * G_vec.comp3 + G_surf[m].comp1 - G_surf[m].comp2 + two * G_surf[m].comp3) < 0.0001 && fabs(kdotn) < 0.95)
                    ////if (fabs(G_vec.comp2 - 0.087611) < 0.0001 && fabs(kdotn) < 0.44)
                    //if (fabs(G_vec.comp2) < 0.0001 && fabs(kdotn) < 0.44)
                    //if (fabs(G_vec.comp1 - 0.219026) < 0.0001 && (kdotn < 0.44 && kdotn > -0.001))
                    if (fabs(G_vec.comp1 - 0.219026) < 0.0001 && (kdotn < 1.5 * 0.87916 && kdotn > 0.87914))
                    //if (fabs(G_vec.comp1 - 0.43805) < 0.0001 && fabs(kdotn) < 0.11)
                    count++;
                   }
                  }
                 }
                }
               }
              }
             }

      int k_break[count];
      double k_value[count], eigval1[count];

      count = 0;
      for (s = 0; s < job->spin_dim; s++) {
        for (i = 0; i < is[0]; i++) {
          for (j = 0; j < is[1]; j++) {
            for (k = 0; k < is[2]; k++) {
              for (l = 0; l < nbands; l++) {
                p = i * is[1] * is[2] + j * is[2] + k;
                q = knet->fbz[p];
                for (m = 0; m < G_surf_num; m++) {
                  for (n = 0; n < G->last_vector; n++) {
                    G_vec.comp1 = G->vec_bi[n].comp1 + knet->cart[p].comp1;
                    G_vec.comp2 = G->vec_bi[n].comp2 + knet->cart[p].comp2;
                    G_vec.comp3 = G->vec_bi[n].comp3 + knet->cart[p].comp3;
                    lambda = e * G_vec.comp1 - d * G_vec.comp2;
                    mu    = -b * G_vec.comp1 + a * G_vec.comp2;
                    //kdotn = (- G_vec.comp1 + G_vec.comp2) / sqrt(two);
                    kdotn = G_vec.comp2;
                    //kdotn = G_vec.comp1;
                    //if (fabs(c * lambda + f * mu - G_vec.comp3) < 0.0001 && fabs(kdotn) < 0.95) {
                    //if (fabs(G_vec.comp1 + G_vec.comp2 - two * G_vec.comp3) < 0.0001 && fabs(kdotn) < 0.95) {
                  //this oneif (fabs(G_vec.comp1 + G_vec.comp2 - two * G_vec.comp3 + G_surf[m].comp1 - G_surf[m].comp2 + two * G_surf[m].comp3) < 0.0001 && fabs(kdotn) < 0.95) {
                    //if (fabs(G_vec.comp2) < 0.0001 && fabs(kdotn) < 0.95) {
                    //if (fabs(G_vec.comp2 - 0.087611) < 0.0001 && fabs(kdotn) < 0.44) {
                    //if (fabs(G_vec.comp2) < 0.0001 && fabs(G_vec.comp1) < 0.44) {
                    //if (fabs(G_vec.comp1) < 0.0001 && fabs(kdotn) < 0.11) {
                    //if (fabs(G_vec.comp1) < 0.0001 && (kdotn < 0.44 && kdotn > -0.001))  {
                    //if (fabs(G_vec.comp1 - 0.219026) < 0.0001 && (kdotn < 0.44 && kdotn > -0.001)) {
                    if (fabs(G_vec.comp1 - 0.219026) < 0.0001 && (kdotn < 1.5 * 0.87916 && kdotn > 0.87914)) {
                    //if (fabs(G_vec.comp1 - 0.43805) < 0.0001 && fabs(kdotn) < 0.11) {
                    k_value[count] = kdotn;
                    eigval1[count] = eigval[nbands * q + l];
                    count++;
                   }
                  }
                 }
                }
               }
              }
             }
            }

  for (i = 1; i < count; ++i) {
    for (j = count - 1; j >= i; --j) {
      if ((k_value[j - 1]) > (k_value[j])) {
     
      temp_k = k_value[j - 1];
      k_value[j - 1] = k_value[j];
      k_value[j] = temp_k;

      temp_k = eigval1[j - 1];
      eigval1[j - 1] = eigval1[j];
      eigval1[j] = temp_k;

     }
    }
   }

   for (i=0;i<G->last_vector;i++) {
     double gmag = sqrt(double_vec_dot(&G->vec_bi[i],&G->vec_bi[i]));
     printf("%3d %10.4lf %10.4lf %10.4lf   %10.4lf\n",i,G->vec_bi[i].comp1,G->vec_bi[i].comp2,G->vec_bi[i].comp3,gmag);
    }

   //for (i=0;i<count;i++) {
     //fprintf(band_project_0,"%3d %10.4lf\n",i,k_value[i]);
    //}

  count1 = 0;
  for (i = 1; i < count; ++i) {
    if (fabs(k_value[i] - k_value[i - 1]) > 0.0001) {
      k_break[count1] = i;
      //printf("%3d %3d\n",count1,k_break[count1]);
      count1++;
    }
   }

   k_start = 0;
   for (i = 0; i < count1; i++) {
     for (j = k_start; j < k_break[i]; j++) {
       for (k = k_break[i] - 1; k >= j; k--) {
         if ((eigval1[k - 1] > eigval1[k])) {
            temp_k = eigval1[k - 1];
            eigval1[k - 1] = eigval1[k];
            eigval1[k] = temp_k;
           }
          }
         }
        k_start = k_break[i];
       }

   for (i = 0; i < count1; i++) {
     for (j = k_start; j < k_break[i]; j++) {
       fprintf(band_project_0,"%10.4lf %10.4lf\n",k_value[j],eigval1[j] * au_to_eV + 6.8226);
      }
      k_start = k_break[i];
     }

}

void plot_points(int nkpoints, int *kplot, KPOINT_TRAN *kpoint, int *bands, VECTOR_DOUBLE *rpoint, int number_Rvec,
VECTOR_DOUBLE *Rvec, Complex *Psirkn, ComplexMatrix *eigvec, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, k;
  int index_i, sheli, shelposi, bfposi, gausposi;
  int i1, i4, j1, j2, j3, k3;
  int nbands = bands[1] - bands[0] + 1;
  int Fdim = atoms->number_of_sh_bfns_in_unit_cell * nkpoints;
  double rsqrd, kdotr;
  double prefac[7], psiexp;
  double r[3];
  double *fiamp, *p_fiamp;
  Complex *Fiamp, *p_Fiamp;

  fiamp = (double *) malloc(atoms->number_of_sh_bfns_in_unit_cell * number_Rvec * sizeof(double));
  if (fiamp == NULL) {
    fprintf(stderr, "ERROR: There is not enough memory for fiamp array \n");
    exit(1);
  }

  Fiamp = (Complex *) malloc(atoms->number_of_sh_bfns_in_unit_cell * nkpoints * sizeof(Complex));
  if (Fiamp == NULL) {
    fprintf(stderr, "ERROR: There is not enough memory for Fiamp array \n");
    exit(1);
  }

   bfposi = 0;
   for (i1 = 0; i1 < atoms->number_of_atoms_in_unit_cell; i1++) {
     shelposi = atoms->shelposn_sh[i1];
     gausposi = atoms->gausposn_sh[i1];
       for (index_i = 0; index_i < atoms->nshel_sh[i1]; index_i++) {
         sheli = shells->type_sh[shelposi + index_i];
           //fprintf(file.out,"i1 %d %d %d %d %d\n",i1,shelposi,gausposi,index_i,sheli);
           for (k3 = 0; k3 < number_Rvec; k3++) {
             r[0] = atoms->cell_vector[i1].comp1 + Rvec[k3].comp1;
             r[1] = atoms->cell_vector[i1].comp2 + Rvec[k3].comp2;
             r[2] = atoms->cell_vector[i1].comp3 + Rvec[k3].comp3;
               rsqrd = (r[0] - rpoint->comp1) * (r[0] - rpoint->comp1) + (r[1] - rpoint->comp2) * (r[1] - rpoint->comp2) +
                       (r[2] - rpoint->comp3) * (r[2] - rpoint->comp3);
                   //fprintf(file.out,"rsq %lf %lf %lf %lf %lf %lf %lf\n",rsqrd,r[0],r[1],r[2],rpoint->comp1,rpoint->comp2,rpoint->comp3);
                     switch (sheli) {
                       case 1:
                         prefac[0] = k_one;
                         break;
                       case 4:
                         prefac[0] = k_one;
                         prefac[1] = rpoint->comp1 - r[0];
                         prefac[2] = rpoint->comp2 - r[1];
                         prefac[3] = rpoint->comp3 - r[2];
                         break;
                       case 3:
                         prefac[0] = rpoint->comp1 - r[0];
                         prefac[1] = rpoint->comp2 - r[1];
                         prefac[2] = rpoint->comp3 - r[2];
                         break;
                       case 5:
                         prefac[0] = ((rpoint->comp3 - r[2]) * (rpoint->comp3 - r[2]) -
                                      (rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) / two -
                                      (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1])/two) / sqrt(three); // 2 zz - xx - yy
                         prefac[1] =  (rpoint->comp1 - r[0]) * (rpoint->comp3 - r[2]); // xz
                         prefac[2] =  (rpoint->comp2 - r[1]) * (rpoint->comp3 - r[2]); // yz
                         prefac[3] = ((rpoint->comp1 - r[0]) * (rpoint->comp1 - r[0]) -
                                      (rpoint->comp2 - r[1]) * (rpoint->comp2 - r[1])) / two; // xx - yy
                         prefac[4] =  (rpoint->comp1 - r[0]) * (rpoint->comp2 - r[1]); // xy
                         break;
                       case 7:
                         fprintf(file.out,"f coefficients not entered yet\n");
                         exit(1);
                         break;
                     } // close switch

      psiexp = k_zero;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
        psiexp += gaussians->c_sh[gausposi + i4] * exp(- gaussians->expo_sh[gausposi + i4] * rsqrd);
        //fprintf(file.out,"psiexp %d %d %lf %lf\n",i4,gausposi +i4,gaussians->c_sh[gausposi+i4],gaussians->expo_sh[gausposi+i4]);
       }
      //fprintf(file.out,"bfposi %5d %5d %5d %5d %5d %5d %5d\n",i1,shelposi,gausposi,index_i,sheli,k3,bfposi);
      for (j3 = 0; j3 < sheli; j3++) {
        p_fiamp = fiamp + (bfposi + j3) * number_Rvec + k3;
        (*p_fiamp) = psiexp * prefac[j3];
        //fprintf(file.out, "bfpos %d psiexp %e prefac %d fiamp %e\n", bfposi+j3, psiexp, j3, (*p_fiamp)) ;
      } // close loop over j3
    } // close loop over k3
    gausposi += shells->ng_sh[shelposi + index_i];
    bfposi += shells->type_sh[shelposi + index_i];
  } // close loop over index_i
  } // close loop over i1

  init_Complex_array(Fdim, Fiamp);
  Complex *p_Psirkn, temp;

  p_Fiamp = Fiamp;
  for (k = 0; k < nkpoints; k++) {
    //fprintf(file.out,"KPOINT %3d %10.4lf %10.4lf %10.4lf\n",k,kpoint->cart[kplot[k]].comp1,kpoint->cart[kplot[k]].comp2, \
    kpoint->cart[kplot[k]].comp3);
    p_fiamp = fiamp;
      for (i = 0; i < atoms->number_of_sh_bfns_in_unit_cell; i++) {
      temp = Complex(k_zero, k_zero);
      for (k3 = 0; k3 < number_Rvec; k3++) {
        kdotr = double_vec_dot(&Rvec[k3], &kpoint->cart[kplot[k]]);
        //kdotr = double_vec_dot(&Rvec[k3], &kpoint->cart[k]);
        //kdotr = double_vec_dot(&Rvec[k3], &kpoint[k].cart);
        temp += (*p_fiamp) * Complex(cos(kdotr), sin(kdotr));
        p_fiamp++;
      }
      (*p_Fiamp) = temp;
      //fprintf(file.out, "%e , %e\n", p_Fiamp->real(), p_Fiamp->imag()) ;
      p_Fiamp++;
    } // close loop over i
  } // close loop over k

  p_Psirkn = Psirkn;
  for (k = 0; k < nkpoints; k++) {
 for (int s = 0; s < job->spin_dim; s++) {
    for (j = 0; j < nbands; j++) {
      p_Fiamp = Fiamp + k * atoms->number_of_sh_bfns_in_unit_cell;
      temp = Complex(k_zero, k_zero);
        for (i = 0; i < atoms->number_of_sh_bfns_in_unit_cell; i++) {
        temp += eigvec->a[(s * nkpoints + k) * nbands + j][i] * (*p_Fiamp);
        ////temp += eigvec->a[k * nbands + j][i] * (*p_Fiamp);
        //fprintf(file.out, "%d %d %d %lf %lf\n",k,j,i,eigvec->a[k*nbands+j][i].real(),eigvec->a[k*nbands+j][i].imag()) ;
        //fprintf(file.out, " %lf , %lf ", p_Fiamp->real(), p_Fiamp->imag()) ;
        p_Fiamp++;
       } // close loop over i
      (*p_Psirkn) = temp;
      p_Psirkn++;
    } // close loop over j
   } // close loop over s
  } // close loop over k

  free(fiamp);
  free(Fiamp);

  return;
}

void plot_points1(double *rho_grid, VECTOR_DOUBLE *rpoint, int grid_size, VECTOR_DOUBLE *Rvec, int number_Rvec, double *F, PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, p, s;
  int index_i, sheli, index_j, shelj, shelposi, shelposj, bfposi, bfposj, gausposi, gausposj;
  int i1, j1, i2, j2, i3, j3, i4, j4;
  int ip, jp, gj;
  int nd1, nd2;
  int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
  int dim3 = atoms->number_of_sh_bfns_in_unit_cell, dim4 = dim3 * dim3;
  int D_ptr, S_ptr;
  double rsqrd;
  double prefac[7], rhoexp;
  double r[3];
  double *rhoamp, *p_rhoamp;

  rhoamp = (double *) malloc(job->spin_dim * dim4 * number_Rvec * sizeof(double));
  if (rhoamp == NULL) { fprintf(file.out, "ERROR: There is not enough memory for rhoamp array \n"); MPI_Finalize(); exit(1); }

  for (i = 0; i < job->spin_dim * dim4 * number_Rvec; i++)
  rhoamp[i] = k_zero;

  for (p = 0; p < pair_p->tot; p++) {
    ip = pair_p->cell1[p];
    jp = pair_p->cell2[p];
    gj = pair_p->latt2[p];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    bfposi = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli = shells->type_sh[index_i];
      bfposj = 0;
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      //fprintf(file.out,"pair_p %3d  %3d %3d   %3d     %3d %3d  %10.4lf %10.4lf %10.4lf   %10.4lf %10.4lf %10.4lf\n",\
      p,ip,jp,gj,shelposj,gausposj,rpoint->comp1,rpoint->comp2,rpoint->comp3,\
      atoms->cell_vector[jp].comp1,atoms->cell_vector[jp].comp2,atoms->cell_vector[jp].comp3);
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj = shells->type_sh[index_j];
        //fprintf(file.out,"jp %3d %3d %3d %3d %3d\n",jp,shelposj,gausposj,index_j,shelj);
        for (j1 = 0; j1 < number_Rvec; j1++) {
          r[0] = atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1 + Rvec[j1].comp1;
          r[1] = atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2 + Rvec[j1].comp2;
          r[2] = atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3 + Rvec[j1].comp3;
          rsqrd = (r[0] - rpoint->comp1) * (r[0] - rpoint->comp1) + (r[1] - rpoint->comp2) * (r[1] - rpoint->comp2) +
                  (r[2] - rpoint->comp3) * (r[2] - rpoint->comp3);
                  //fprintf(file.out,"rsq %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n", \
                  rsqrd,r[0],r[1],r[2],rpoint->comp1,rpoint->comp2,rpoint->comp3);
         if (rsqrd > 100.0) continue;

         switch (shelj) {

             case 1:
               prefac[0] = k_one;
               break;
             case 4:
               prefac[0] = k_one;
               prefac[1] = rpoint->comp1 - r[0];
               prefac[2] = rpoint->comp2 - r[1];
               prefac[3] = rpoint->comp3 - r[2];
               break;
             case 3:
               prefac[0] = rpoint->comp1 - r[0];
               prefac[1] = rpoint->comp2 - r[1];
               prefac[2] = rpoint->comp3 - r[2];
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
               fprintf(file.out,"f coefficients not entered yet\n");
               exit(1);
               break;

           } // close switch

  rhoexp = k_zero;
  for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
    rhoexp += gaussians->c_sh[gausposj + j4] * exp(- gaussians->expo_sh[gausposj + j4] * rsqrd);
    //fprintf(file.out,"rhoexp %3d %3d %10.4lf %10.4lf   %12.3e\n",j4,gausposj +j4,gaussians->c_sh[gausposj+j4],gaussians->expo_sh[gausposj+j4],rhoexp);
   }
   if (rhoexp < 1.0e-6) continue;
   //fprintf(file.out,"bfposj %5d %5d %5d %5d %5d %5d\n",shelposj,gausposj,index_j,j1,shelj,bfposj);
   for (s = 0; s < job->spin_dim; s++) {
     S_ptr = s * job->dimf;
     for (i3 = 0; i3 < sheli; i3++) {
       for (j3 = 0; j3 < shelj; j3++) {
         p_rhoamp = rhoamp + s * dim4 * number_Rvec + (atoms->bfnposn_sh[ip] + bfposi + i3) * dim3 * number_Rvec + \
         (atoms->bfnposn_sh[jp] + bfposj + j3) * number_Rvec + j1;
         D_ptr = S_ptr + pair_p->off[pair_p->ptr[gj * dim2 + ip * dim1 + jp]];
         //fprintf(file.out,"%3d %3d %3d   %3d   %3d %3d   %3d    %3d %3d   %10.4lf  %12.3e %12.3e %12.3e\n", \
         s,bfposj+j3,j1,s * atoms->number_of_sh_bfns_in_unit_cell * number_Rvec + (atoms->bfnposn_sh[jp] + bfposj + j3) * number_Rvec + j1, \
         ip, jp, gj, pair_p->ptr[gj * dim2 + ip * dim1 + jp],D_ptr,rsqrd,rhoexp,prefac[j3],F[D_ptr + bfposj + j3]);
         (*p_rhoamp) += rhoexp * prefac[j3] * F[D_ptr + (bfposi + i3) * nd2 + bfposj + j3];
         if ((bfposi+i3)*nd2+bfposj+j3 > nd1 * nd2) fprintf(file.out,"ERROR  %3d %3d %3d   %3d\n",\
         atoms->bfnposn_sh[jp], bfposj, j3, atoms->number_of_sh_bfns_in_unit_cell);
        } // close loop over j3
       } // close loop over i3
      } // close loop over s
     } // close loop over j1
    gausposj += shells->ng_sh[index_j];
    bfposj += shells->type_sh[index_j];
   } // close loop over index_j
   gausposi += shells->ng_sh[index_i];
   bfposi += shells->type_sh[index_i];
  } // close loop over index_i
 } // close loop over p

   //for (i = 0; i < dim4 * number_Rvec; i++)
   //fprintf(file.out,"%5d %12.3E\n",i,rhoamp[i]);

   for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
     bfposi = 0;
     shelposi = atoms->shelposn_sh[ip];
     gausposi = atoms->gausposn_sh[ip];
     for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
       sheli = shells->type_sh[index_i];
       for (i1 = 0; i1 < number_Rvec; i1++) {
         r[0] = atoms->cell_vector[ip].comp1 + Rvec[i1].comp1;
         r[1] = atoms->cell_vector[ip].comp2 + Rvec[i1].comp2;
         r[2] = atoms->cell_vector[ip].comp3 + Rvec[i1].comp3;
         rsqrd = (r[0] - rpoint->comp1) * (r[0] - rpoint->comp1) + (r[1] - rpoint->comp2) * (r[1] - rpoint->comp2) +
                 (r[2] - rpoint->comp3) * (r[2] - rpoint->comp3);
                 //fprintf(file.out,"rsq %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n", \
                 rsqrd,r[0],r[1],r[2],rpoint->comp1,rpoint->comp2,rpoint->comp3);
         if (rsqrd > 100.0) continue;

           switch (sheli) {

             case 1:
               prefac[0] = k_one;
               break;
             case 4:
               prefac[0] = k_one;
               prefac[1] = rpoint->comp1 - r[0];
               prefac[2] = rpoint->comp2 - r[1];
               prefac[3] = rpoint->comp3 - r[2];
               break;
             case 3:
               prefac[0] = rpoint->comp1 - r[0];
               prefac[1] = rpoint->comp2 - r[1];
               prefac[2] = rpoint->comp3 - r[2];
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
               fprintf(file.out,"f coefficients not entered yet\n");
               exit(1);
               break;

           } // close switch

  rhoexp = k_zero;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    rhoexp += gaussians->c_sh[gausposi + i4] * exp(- gaussians->expo_sh[gausposi + i4] * rsqrd);
    //fprintf(file.out,"rhoexp %3d %3d %10.4lf %10.4lf   %12.3e\n",i4,gausposi +i4,gaussians->c_sh[gausposi+i4],gaussians->expo_sh[gausposi+i4],rhoexp);
   }
   if (rhoexp < 1.0e-6) continue;
   //fprintf(file.out,"bfposi %5d %5d %5d %5d %5d %5d\n",shelposi,gausposi,index_i,i1,sheli,bfposi);
   for (s = 0; s < job->spin_dim; s++) {
     for (i3 = 0; i3 < sheli; i3++) {
       for (jp = 0; jp < atoms->number_of_atoms_in_unit_cell; jp++) {
         bfposj = 0;
         shelposj = atoms->shelposn_sh[jp];
         gausposj = atoms->gausposn_sh[jp];
         for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
           shelj = shells->type_sh[index_j];
           for (j1 = 0; j1 < number_Rvec; j1++) {
             for (j3 = 0; j3 < shelj; j3++) {
               if (atoms->bfnposn_sh[jp] + bfposj + j3 > atoms->number_of_sh_bfns_in_unit_cell) fprintf(file.out,"ERROR  %3d %3d %3d   %3d\n",\
               atoms->bfnposn_sh[jp], bfposj, j3, atoms->number_of_sh_bfns_in_unit_cell);
               p_rhoamp = rhoamp + s * dim4 * number_Rvec + (atoms->bfnposn_sh[ip] + bfposi + i3) * dim3 * number_Rvec + \
               (atoms->bfnposn_sh[jp] + bfposj + j3) * number_Rvec + j1;
               //p_rhoamp = rhoamp + s * atoms->number_of_sh_bfns_in_unit_cell * number_Rvec + (atoms->bfnposn_sh[jp] + bfposj + j3) * number_Rvec + j1;
               //fprintf(file.out,"%3d %3d %3d   %3d  %12.3e %12.3e %12.3e\n", \
               s,atoms->bfnposn_sh[jp]+bfposj+j3,j1,s*atoms->number_of_sh_bfns_in_unit_cell*number_Rvec+(atoms->bfnposn_sh[jp]+bfposj+j3)*number_Rvec+j1,\
               rhoexp,prefac[i3],*p_rhoamp);
               *(rho_grid + s * grid_size) += rhoexp * prefac[i3] * *p_rhoamp;
               // need spin offset here
              } // close loop over j3
             } // close loop over i3
            gausposj += shells->ng_sh[index_j];
            bfposj += shells->type_sh[index_j];
           } // close loop over index_j
          } // close loop over jp
         } // close loop over i3
        } // close loop over s
       } // close loop over i1
      gausposi += shells->ng_sh[index_i];
      bfposi += shells->type_sh[index_i];
     } // close loop over index_i
    } // close loop over ip

  free(rhoamp);

}

void plot_correlated_electron_hole(int grid_par[3], VECTOR_DOUBLE *points, int num_points, VECTOR_DOUBLE *points1, int bse_state[2], FERMI* fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  int i, j, m, n, p, s;
  int dimk = 1;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;
  PAIR_TRAN pair_p;
  //double *bse_eigenvalue;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity >= 1)
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
  // * Read and broadcast scf eigenvectors                                                    *
  // ******************************************************************************************

  int vector_size;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim2 = dim1 * dim1;
  char zz2[24] = "scf_evectors";
  size_t result;
  FILE *scf_evectors;
  DoubleMatrix *eigvec;
  ComplexMatrix *scf_eigenvectors;
  dimk = job->spin_dim * nbands;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&scf_eigenvectors,&dimk,&dim1,job);
  AllocateDoubleMatrix(&eigvec,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);

  //strcpy(buf2,file.directory1);
  //strcat(buf2,yy);
  //MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;

  if (job->taskid == 0) {
  scf_evectors = fopen(zz2, "rb");
  fseek(scf_evectors, dim1 * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
  result = fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
  fclose(scf_evectors);

//ResetDoubleMatrix(eigvec);
//eigvec->a[0][10]  = 1.0;
//eigvec->a[1][0]  = 1.0;

  for (i = 0; i < nbands; i++) {
    for (j = 0; j < dim1; j++) {
      eigvec->a[i][j] = (scf_eigenvectors->a[i][j]).real();
      //fprintf(file.out,"eigvec %3d %3d %10.4lf\n",i,j,eigvec->a[i][j]);
     }
    }
   }

  MPI_Bcast(&eigvec->a[0][0],nbands * dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Read and broadcast bse eigenvectors                                                    *
  // ******************************************************************************************

  int ntransitions;
  int num_bse_states;

  ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  //ntransitions = (nbands - fermi->occupied[0]) * fermi->occupied[0];
  num_bse_states = bse_state[1] - bse_state[0] + 1;
  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in RPA_molecule. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }

  double *bse_eigenvalues;
  DoubleMatrix *bse_eigenvector;
  AllocateDoubleMatrix(&bse_eigenvector,&num_bse_states,&ntransitions,job);
  ResetDoubleMatrix(bse_eigenvector);

  if (job->taskid == 0) {
  char xx[20] = "./bse_eigenvectors";
  char yy[20] = "./bse_eigenvalues";
  size_t result;
  FILE *bse_evectors, *bse_evalues;
  bse_evectors = fopen(xx, "rb");
  bse_evalues  = fopen(yy, "rb");
  fseek(bse_evectors, (bse_state[0] - 1) * ntransitions * sizeof(double), SEEK_SET);
  result = fread(&bse_eigenvector->a[0][0], sizeof(double), num_bse_states * ntransitions, bse_evectors);
  fseek(bse_evalues, (bse_state[0] - 1) * sizeof(double), SEEK_SET);
  //fread(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
  fclose(bse_evectors);
  fclose(bse_evalues);

//ResetDoubleMatrix(bse_eigenvector);
//bse_eigenvector->a[0][0] = 1.0;

  for (i = 0; i < num_bse_states; i++) {
    for (j = 0; j < ntransitions; j++) {
      //fprintf(file.out,"bse_eigvec %3d %3d %10.4lf\n",i,j,bse_eigenvector->a[i][j]);
     }
    }
 }
  MPI_Bcast(&bse_eigenvector->a[0][0], num_bse_states * ntransitions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(bse_eigenvalues,num_bse_states,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // Calculate wavefunctions on grid and electron hole correlations                           *
  // ******************************************************************************************
 
  int count;
  double *psiamp1_buffer, *psiamp2_buffer;
  double *psiamp1, *psiamp2;
  double *Psiamp;
  char buf[22] = "bse_isosurface.xsf";

  VECTOR_DOUBLE Rvec1[radmx], Rvec2[radmx];

  printf("%3d %3d %3d %3d %3d %3d %3d\n",num_bse_states,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],ntransitions);
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

  fprintf(file.out,"particle grids\n");
  for (n = 0; n < num_points; n++) {
    wavefunction_gridpoint(&psiamp2[n * nvir],&points1[n],Rvec2,&eigvec->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
   }

  fprintf(file.out,"hole grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(psiamp1,&grid[p],Rvec1,&eigvec->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
    for (n = 0; n < num_points; n++) {
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < num_bse_states; m++) {
          count = 0;
          for (i = 0; i < nocc; i++) {
            for (j = 0; j < nvir; j++) {
              Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
              bse_eigenvector->a[m][count] * psiamp1[s * nocc + i] * psiamp2[s * num_points * nvir + n * nvir + j];
              //bse_eigenvector->a[m][count] * psiamp1[s * nocc + i] ;
              //fprintf(file.out,"   %3d %3d %3d %3d %3d %12.6lf %10.4lf %10.4lf %12.6lf   %10.4lf %10.4lf %10.4lf\n",\
              p,n,m,i,j,bse_eigenvector->a[m][count],psiamp1[i],psiamp2[n * nvir + j],\
              Psiamp[m * num_points * grid_size + n * grid_size + p],grid[p].comp1,grid[p].comp2,grid[p].comp3);
              count++;
             }
            }
           //fprintf(file.out,"\n");
           }
          }
         }
        }

  write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,num_bse_states,num_points,crystal,atoms,shells,gaussians,job,file);

  free(grid);
  free(Rvec);
  free(Psiamp);
  free(psiamp1);
  free(psiamp2);
  free(psiamp1_buffer);
  free(psiamp2_buffer);
  free_PAIR_TRAN(&pair_p,job);
  DestroyDoubleMatrix(&bse_eigenvector,job);
  DestroyDoubleMatrix(&eigvec,job);
  DestroyComplexMatrix(&scf_eigenvectors,job);
  //DestroyDoubleArray(&bse_eigenvalue,&ntransitions,job);

}

void plot_bse_transition_density(int grid_par[3], VECTOR_DOUBLE *points, VECTOR_DOUBLE *points1, int bse_state[2], FERMI* fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  int i, j, k, m, n, p, s;
  int dimk = 1;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;
  PAIR_TRAN pair_p;
  //double *bse_eigenvalue;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity >= 1)
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
  // * Read and broadcast scf eigenvectors                                                    *
  // ******************************************************************************************

  int vector_size;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim2 = dim1 * dim1;
  char zz2[24] = "scf_evectors";
  size_t result;
  FILE *scf_evectors;
  DoubleMatrix *eigvec;
  ComplexMatrix *scf_eigenvectors;
  dimk = job->spin_dim * nbands;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&scf_eigenvectors,&dimk,&dim1,job);
  AllocateDoubleMatrix(&eigvec,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);

  if (job->taskid == 0) {
  scf_evectors = fopen(zz2, "rb");
  fseek(scf_evectors, dim1 * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
  result = fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
  fclose(scf_evectors);

  for (i = 0; i < nbands; i++) {
    for (j = 0; j < dim1; j++) {
      eigvec->a[i][j] = (scf_eigenvectors->a[i][j]).real();
      //fprintf(file.out,"eigvec %3d %3d %10.4lf\n",i,j,eigvec->a[i][j]);
     }
    }
   }

  //MPI_Bcast(&eigvec->a[0][0],nbands * dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Read and broadcast bse eigenvectors                                                    *
  // ******************************************************************************************

  int ntransitions;
  int num_bse_states;

  ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  //ntransitions = (nbands - fermi->occupied[0]) * fermi->occupied[0];
  num_bse_states = bse_state[1] - bse_state[0] + 1;
  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in RPA_molecule. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }

  double *bse_eigenvalues;
  DoubleMatrix *bse_eigenvector;
  AllocateDoubleMatrix(&bse_eigenvector,&num_bse_states,&ntransitions,job);
  ResetDoubleMatrix(bse_eigenvector);

  char xc[22] = "/bse_eigenvectors";
  char bufbse[120];
  FILE *bse_evectors;
  strcpy(bufbse,file.bse_eigvec);
  strcat(bufbse,xc);
  //printf("bse_eigenvectors directory: %s %s\n",file.bse_eigvec,bufbse);
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,bufbse,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, (bse_state[0] - 1) * ntransitions * sizeof(double), MPI_SEEK_SET) ;
  //MPI_File_read(fh, &bse_eigenvector->a[0][0], num_bse_states * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  //MPI_File_close(&fh);

  char yy[20] = "bse_eigenvalues";
  AllocateDoubleArray(&bse_eigenvalues,&ntransitions,job);
  ResetDoubleArray(bse_eigenvalues,&ntransitions);
  read_scf_GW_eigenvalues(bse_eigenvalues, 0, ntransitions, yy, job, file);
  //for (i = 0; i < ntransitions; i++) printf("%10.4f\n",bse_eigenvalues[i]);

  job->bse_tda = 0;
  job->bse_lim = ntransitions;
  if (job->taskid == 0 && job->bse_tda == 0) {
  for (i = 0; i < bse_state[1] - bse_state[0] + 1; i++) {
  MPI_File_seek(fh, (job->bse_lim - bse_state[0] - i) * ntransitions * sizeof(double), MPI_SEEK_SET) ;
  MPI_File_read(fh, &bse_eigenvector->a[i][0], ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
 }
  for (j = 0; j < bse_eigenvector->iRows; j++) {
    for (k = 0; k < bse_eigenvector->iCols; k++) {
      //fprintf(file.out,"%3d %3d %10.4f\n",j,k,bse_eigenvector->a[j][k]);
      bse_eigenvector->a[j][k] /= two * sqrt(bse_eigenvalues[j + bse_state[0] - 1]);
     }
    }
   }

  if (job->taskid == 0) {
/*
  char xx[20] = "./bse_eigenvectors";
  char yy[20] = "./bse_eigenvalues";
  FILE *bse_evectors, *bse_evalues;
  bse_evectors = fopen(xx, "rb");
  //bse_evalues  = fopen(yy, "rb");
  fseek(bse_evectors, (bse_state[0] - 1) * ntransitions * sizeof(double), SEEK_SET);
  fread(&bse_eigenvector->a[0][0], sizeof(double), num_bse_states * ntransitions, bse_evectors);
  //fseek(bse_evalues, (bse_state[0] - 1) * sizeof(double), SEEK_SET);
  //fread(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
  fclose(bse_evectors);
  //fclose(bse_evalues);
*/

//ResetDoubleMatrix(bse_eigenvector);
//bse_eigenvector->a[0][0] = 1.0;

  for (i = 0; i < num_bse_states; i++) {
    for (j = 0; j < ntransitions; j++) {
      //fprintf(file.out,"bse_eigvec %3d %3d %10.4lf\n",i,j,bse_eigenvector->a[i][j]);
     }
    }
 }
  //MPI_Bcast(&bse_eigenvector->a[0][0], num_bse_states * ntransitions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(bse_eigenvalues,num_bse_states,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // Calculate wavefunctions on grid and electron hole correlations                           *
  // ******************************************************************************************
 
  int count;
  double *psiamp1_buffer, *psiamp2_buffer;
  double *psiamp1, *psiamp2;
  double *Psiamp;
  char buf[38] = "bse_transition_density_isosurface.xsf";

  VECTOR_DOUBLE Rvec1[radmx], Rvec2[radmx];

  //printf("%3d %3d %3d %3d %3d %3d %3d\n",num_bse_states,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],ntransitions);
  Rvec1[0].comp1 = 0.0;
  Rvec1[0].comp2 = 0.0;
  Rvec1[0].comp3 = 0.0;
  Rvec2[0].comp1 = 0.0;
  Rvec2[0].comp2 = 0.0;
  Rvec2[0].comp3 = 0.0;

  Psiamp = (double *) calloc(job->spin_dim * num_bse_states * grid_size, sizeof(double));
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

  //fprintf(file.out,"particle/hole grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(psiamp1,&grid[p],Rvec1,&eigvec->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
    wavefunction_gridpoint(psiamp2,&grid[p],Rvec1,&eigvec->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
    for (s = 0; s < job->spin_dim; s++) {
      for (m = 0; m < num_bse_states; m++) {
        count = 0;
        for (i = 0; i < nocc; i++) {
          for (j = 0; j < nvir; j++) {
            Psiamp[s * num_bse_states * grid_size + m * grid_size + p] += bse_eigenvector->a[m][count] * \
            psiamp1[s * nocc + i] * psiamp2[s * nvir + j];
            //fprintf(file.out,"   %3d %3d %3d %3d %12.6lf %14.8lf %14.8lf %16.10lf   %10.4lf %10.4lf %10.4lf\n",\
            p,m,i,j,bse_eigenvector->a[m][count],psiamp1[i],psiamp2[j],\
            Psiamp[m * grid_size + p],grid[p].comp1,grid[p].comp2,grid[p].comp3);
            count++;
           }
          }
         //fprintf(file.out,"\n");
         }
        }
       }

  write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,num_bse_states,1,crystal,atoms,shells,gaussians,job,file);

  free(grid);
  free(Rvec);
  free(Psiamp);
  free(psiamp1);
  free(psiamp2);
  free(psiamp1_buffer);
  free(psiamp2_buffer);
  free_PAIR_TRAN(&pair_p,job);
  DestroyDoubleMatrix(&bse_eigenvector,job);
  DestroyDoubleMatrix(&eigvec,job);
  DestroyComplexMatrix(&scf_eigenvectors,job);
  //DestroyDoubleArray(&bse_eigenvalue,&ntransitions,job);

}

void plot_bse_transition_density_crystal(double *grid_sum, int bse_state[2],  int *grid_par, VECTOR_DOUBLE *points, int num_points, VECTOR_DOUBLE *points1, FERMI *fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurface of transition charge density for BSE in crystals         *
  // ******************************************************************************************

  int i, j, k, m, n, p, s;
  int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;
  //PAIR_TRAN pair_p;

  // ******************************************************************************************
  // * Generate knet                                                                          *
  // ******************************************************************************************

  KPOINT_TRAN knet;
  count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  if (job->kss == 0)      fermi->nkunique = knet.nktot;
  else if (job->kss == 1) fermi->nkunique = knet.unique;
  allocate_k_points(&knet,crystal,job,file);
  generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  allocate_fermi(fermi,atoms,job,file);
  fermi->knet = &knet;
  fermi->nktot = knet.nktot;
  fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  //count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  //generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //if (job->taskid == 0 && job->verbosity > 1)
  //print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  int grid_size, spin_grid_size;
  int radmx = 2;
  int plot_dim = 0;
  int number_Rvec;
  int ngrid_points;
  VECTOR_DOUBLE *grid, *Rvec;

  grid_size = grid_par[0] * grid_par[1] * grid_par[2];
  ////spin_grid_size = job->spin_dim * grid_size;
  Rvec_grid_size(&number_Rvec, radmx, crystal, file);

  grid = (VECTOR_DOUBLE *) malloc(grid_size * sizeof(VECTOR_DOUBLE));
  if (grid == NULL) { fprintf(file.out, "Cannot open memory for grid array\n"); exit(1); }

  Rvec = (VECTOR_DOUBLE *) malloc(number_Rvec * sizeof(VECTOR_DOUBLE));
  if (Rvec == NULL) { fprintf(file.out, "There is not enough memory for Rvec array\n"); exit(1); }

  generate_Rvec_grid(number_Rvec, radmx, Rvec, crystal, file);
  calc_grid1(grid_par, grid, grid_size, &ngrid_points, points, job, file);

  //for(i = 0; i < number_Rvec; i++) \
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);

  // ******************************************************************************************
  // * Open scf eigenvector file and allocate memory for eigenvectors                         *
  // ******************************************************************************************

  MPI_File fh;
  int dimk;
  int vector_size = nbands * nbfn;
  int block_size = vector_size * sizeof(Complex);
  char buf2[110], filename[15] = "/scf_evec";
  strcpy(buf2,file.bse_eigvec);
  strcat(buf2,filename);
  //printf("scf eigenvectors directory: %s \n",buf2);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  ComplexMatrix *scf_eigvec_0, *scf_eigvec_1, *scf_eigvec_2;
  dimk = job->spin_dim * nbands;

  AllocateComplexMatrix(&scf_eigvec_0,&nbands,&nbfn,job);
  AllocateComplexMatrix(&scf_eigvec_1,&nbands,&nbfn,job);
  AllocateComplexMatrix(&scf_eigvec_2,&nbands,&nbfn,job);

  // ******************************************************************************************
  // * Open bse eigenvector file and allocate memory for eigenvectors                         *
  // ******************************************************************************************

  int ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1) * fermi->nktot;
  int num_bse_states = bse_state[1] - bse_state[0] + 1;
spin_grid_size = job->spin_dim * num_bse_states * num_points * grid_size;

  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in BSE_PLOT. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }

  char xc[22] = "/bse_eigenvectors";
  char bufbse[120];
  FILE *bse_evectors;
  long memsize = num_bse_states * ntransitions;
  ComplexMatrix *bse_eigenvector;

  AllocateComplexMatrix(&bse_eigenvector,&num_bse_states,&ntransitions,job);
  ResetComplexMatrix(bse_eigenvector);
  strcpy(bufbse,file.bse_eigvec);
  strcat(bufbse,xc);

  MPI_File gh;
  MPI_File_open(MPI_COMM_WORLD,bufbse,MPI_MODE_RDONLY,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, (bse_state[0] - 1) * ntransitions * sizeof(Complex), MPI_SEEK_SET) ;
  MPI_File_read(gh, &bse_eigenvector->a[0][0], 2 * memsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&gh);

  //ResetComplexMatrix(bse_eigenvector);
  //(bse_eigenvector->a[0][2]).real() = 1.0;
  //for (i = 0; i < 64; i++)
  //(bse_eigenvector->a[1][9*i + 2]).real() = 1.0;

  for (j = 0; j < bse_eigenvector->iRows; j++) {
    for (k = 0; k < bse_eigenvector->iCols; k++) {
      //fprintf(file.out,"%3d %3d %10.4f %10.4f\n",j,k,(bse_eigenvector->a[j][k]).real(),(bse_eigenvector->a[j][k]).imag());
      //bse_eigenvector->a[j][k] /= two * sqrt(bse_eigenvalues[j + bse_state[0] - 1]);
     }
    //fprintf(file.out,"\n");
    }

  int begin_k[job->numtasks], end_k[job->numtasks];
  int total_tasks = fermi->knet->unique;
  mpi_begin_end(begin_k,end_k,total_tasks,job->numtasks,job,file);
  //printf("%3d %3d %3d\n",job->taskid,begin_k[job->taskid],end_k[job->taskid]);

  // ******************************************************************************************
  // Calculate transition density on grid                                                     *
  // ******************************************************************************************
 
  int count, offset;
  int q, q1, q3;
  int count_bz, k_bz, kq_bz, fbz1, fbz2;
  Complex *psiamp1_buffer, *psiamp2_buffer;
  Complex *Psiamp1, *Psiamp1_buffer, *psiamp1, *psiamp2;
  double *Psiamp;
  char buf[38] = "bse_transition_density_isosurface.xsf";
  VECTOR_DOUBLE Rvec1[radmx], Rvec2[radmx];

  //printf("%3d %3d %3d %3d %3d %3d %3d\n",num_bse_states,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],ntransitions);
  Rvec1[0].comp1 = 0.0;
  Rvec1[0].comp2 = 0.0;
  Rvec1[0].comp3 = 0.0;
  //Rvec2[0].comp1 = 0.0;
  //Rvec2[0].comp2 = 0.0;
  //Rvec2[0].comp3 = 0.0;

  Psiamp = (double *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(double));
  if (Psiamp == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp \n");
  MPI_Finalize(); exit(1); }

  Psiamp1 = (Complex *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(Complex));
  if (Psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp1 \n");
  MPI_Finalize(); exit(1); }

  Psiamp1_buffer = (Complex *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(Complex));
  if (Psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp1 = (Complex *) calloc(job->spin_dim * grid_size * fermi->nktot * nocc, sizeof(Complex));
  if (psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp2 = (Complex *) calloc(job->spin_dim * num_points * fermi->nktot * nvir, sizeof(Complex));
  if (psiamp2 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2 \n");
  MPI_Finalize(); exit(1); }

  //psiamp1_buffer = (Complex *) calloc(job->spin_dim * grid_size * fermi->nktot * nocc, sizeof(Complex));
  //if (psiamp1_buffer == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1_buffer \n");
  //MPI_Finalize(); exit(1); }

  //psiamp2_buffer = (Complex *) calloc(job->spin_dim * num_points * fermi->nktot * nvir, sizeof(Complex));
  //if (psiamp2_buffer == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2_buffer \n");
  //MPI_Finalize(); exit(1); }

  int ibz, k1, vector_offset;
long count_k;
  int buffer1_size = job->spin_dim * num_bse_states * num_points * grid_size;
  //int buffer2_size = job->spin_dim * grid_size * fermi->nktot * nvir;

/*
  int ip, jp;
  double kdotr, rsqrd;
  double r[3];
  Complex sumr, sumi;
  Complex temp;
  VECTOR_DOUBLE *p_Rvec;
  p = 3431;
  sumr = 0.0; sumi = 0.0;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  for (p = 3429; p < 3440; p++) {
    //sumr = 0.0; sumi = 0.0;
    //count_k = 0;
        p_Rvec = Rvec;
        for (jp = 0; jp < number_Rvec; jp++) {
          sumr = 0.0; sumi = 0.0;
    count_k = 0;
    for (k = 0; k < fermi->knet->unique; k++) {
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        kdotr = double_vec_dot(&knet.cart[k_bz],p_Rvec);
        temp = Complex(cos(kdotr), sin(kdotr));
        //p_Rvec = Rvec;
        //for (jp = 0; jp < number_Rvec; jp++) {
          //kdotr = double_vec_dot(&knet.cart[k_bz],p_Rvec);
          //temp = Complex(cos(kdotr), sin(kdotr));
          //if (jp == 13) fprintf(file.out,"k_bz %3d R %10.4f %10.4f %10.4f  %3d kdotr %10.4f  %10.4f %10.4f\n",\
          k_bz,p_Rvec->comp1,p_Rvec->comp2,p_Rvec->comp3,jp,kdotr,cos(kdotr),sin(kdotr));
          for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
          r[0] = atoms->cell_vector[ip].comp1 + p_Rvec->comp1;
          r[1] = atoms->cell_vector[ip].comp2 + p_Rvec->comp2;
          r[2] = atoms->cell_vector[ip].comp3 + p_Rvec->comp3;
          rsqrd = (r[0] - grid[p].comp1) * (r[0] - grid[p].comp1) + \
                  (r[1] - grid[p].comp2) * (r[1] - grid[p].comp2) + \
                  (r[2] - grid[p].comp3) * (r[2] - grid[p].comp3);
          sumr += exp(-0.2 * rsqrd) * cos(kdotr);
          sumi += exp(-0.2 * rsqrd) * sin(kdotr);
        }
      fprintf(file.out,"k_bz %3d R %10.4f %10.4f %10.4f %3d kdotr %10.4f  %10.4f %10.4f  %10.4f %10.4f %10.4f %10.4f\n",\
      k_bz,p_Rvec->comp1,p_Rvec->comp2,p_Rvec->comp3,jp,kdotr,cos(kdotr),sin(kdotr),grid[p].comp1,grid[p].comp2,grid[p].comp3,rsqrd);
        //p_Rvec++;
        count_k++;
        }
        //count_k++;
        }
        p_Rvec++;
  printf("%3d %3d  %10.4f %10.4f   %10.4f %10.4f %10.4f\n",jp,p,sumr,sumi,grid[p].comp1,grid[p].comp2,grid[p].comp3);
  fprintf(file.out,"%10.4f %10.4f   %10.4f %10.4f %10.4f\n\n",sumr,sumi,grid[p].comp1,grid[p].comp2,grid[p].comp3);
        }
  //printf("%10.4f %10.4f   %10.4f %10.4f %10.4f\n",sumr,sumi,grid[p].comp1,grid[p].comp2,grid[p].comp3);
        }
  MPI_Finalize();
  exit(0);
*/

  ResetComplexArray(Psiamp1,&spin_grid_size);
  ResetComplexArray(Psiamp1_buffer,&spin_grid_size);

  //fprintf(file.out,"particle/hole grids\n");
  //printf("particle/hole grids\n");
  offset = (fermi->bands[0] - 1) * nbfn;
  //printf("offset %3d %3d\n",offset,job->spin_dim);
  for (s = 0; s < job->spin_dim; s++) {
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) {
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; }
        //printf("b %3d %3d %3d %3d\n",begin_k[job->taskid],k,count_k,fermi->knet->num[k]);
       }
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      //printf("0 %3d %3d %3d\n",job->taskid,k,count_k);
    ////for (k = 0; k < fermi->knet->unique; k++) {
      //ibz = fermi->knet->ibz[k];
      MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_2->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
        //for (int jj=0;jj<nbfn;jj++) fprintf(file.out,\
        "EIG %3d %3d %3d %10.4f %10.4f\n",\
        k,k_bz,jj,\
        (scf_eigvec_2->a[0][jj]).real(),(scf_eigvec_2->a[0][jj]).imag());
        //"EIG %3d %3d %3d %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f %10.4f %10.4f %10.4f %10.4f\n",\
        //(scf_eigvec_2->a[0][jj]).real(),(scf_eigvec_2->a[0][jj]).imag(),\
        (scf_eigvec_2->a[1][jj]).real(),(scf_eigvec_2->a[1][jj]).imag(),\
        (scf_eigvec_2->a[2][jj]).real(),(scf_eigvec_2->a[2][jj]).imag(),\
        (scf_eigvec_2->a[3][jj]).real(),(scf_eigvec_2->a[3][jj]).imag(),\
        (scf_eigvec_2->a[4][jj]).real(),(scf_eigvec_2->a[4][jj]).imag(),\
        (scf_eigvec_2->a[5][jj]).real(),(scf_eigvec_2->a[5][jj]).imag());
        //printf("particle grid num_points %3d\n",num_points);
        for (n = 0; n < num_points; n++) {
          wavefunction_gridpoint_crystal(&psiamp2[n * fermi->nktot * nvir + k_bz * nvir],&points1[n],1,Rvec1,k_bz,\
          fermi->knet,&scf_eigvec_2->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
          //for (j = 0; j < nvir; j++) \
          fprintf(file.out,"k %3d k1 %3d k_bz %3d n %3d vir %3d %14.8lf %14.8lf %10.4lf %10.4lf %10.4lf\n",\
          k,k1,k_bz,n,j,(psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]).real(),\
          (psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]).imag(),points1[n].comp1,points1[n].comp2,points1[n].comp3);
         }
          //fprintf(file.out,"\n");
        count_k++;
       } // close loop on k1
      } // close loop on k
       //ResetComplexMatrix(bse_eigenvector);
      //(bse_eigenvector->a[0][10]).real() = 1.0;
      //(bse_eigenvector->a[1][21]).real() = 1.0;
  //fprintf(file.out,"hole grids\n");
  count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) {
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; }}

  //for (k = 0; k < 1; k++) {
  ////for (k = 0; k < fermi->knet->unique; k++) {
  for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
    //ibz = fermi->knet->ibz[k];
    MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
    MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
      k_bz = fermi->knet->bz[count_k];
      rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_1->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
      //for (int jj=0;jj<nbfn;jj++) fprintf(file.out,"%3d %3d  %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",\
      k_bz,kq_bz,jj,\
      (scf_eigvec_1->a[10][jj]).real(),(scf_eigvec_1->a[10][jj]).imag(),\
      (scf_eigvec_1->a[6][jj]).real(),(scf_eigvec_1->a[6][jj]).imag(),\
      (scf_eigvec_1->a[7][jj]).real(),(scf_eigvec_1->a[7][jj]).imag(),\
      (scf_eigvec_1->a[8][jj]).real(),(scf_eigvec_1->a[8][jj]).imag());
      //printf("particle grid num_points %3d\n",num_points);
      vector_offset = k_bz * nocc * nvir;
      for (p = 0; p < grid_size; p++) {
        //fprintf(file.out,"p = %6d\n",p);
        wavefunction_gridpoint_crystal(&psiamp1[p * fermi->nktot * nocc + k_bz * nocc],&grid[p],number_Rvec,Rvec,k_bz,fermi->knet,\
        &scf_eigvec_1->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
       }
      count_k++;
     }
    printf("c %3d %3d\n",job->taskid,k);
    }

  ////for (p = 0; p < grid_size; p++) {
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) {
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; }}

    ////for (k = 0; k < fermi->knet->unique; k++) {
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      //printf("1 %3d %3d\n",job->taskid,k);
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        //for (i = 0; i < nocc; i++) \
        fprintf(file.out,"k %3d k1 %3d k_bz %3d p %5d o %3d  %14.8f %14.8f %10.4f%10.4f%10.4f\n",\
        k,k1,k_bz,p,i,(psiamp1[i]).real(),(psiamp1[i]).imag(),grid[p].comp1,grid[p].comp2,grid[p].comp3);
        vector_offset = k_bz * nocc * nvir;
        for (m = 0; m < num_bse_states; m++) {
          for (n = 0; n < num_points; n++) {
  for (p = 0; p < grid_size; p++) {
            count = 0;
            for (i = 0; i < nocc; i++) {
              for (j = 0; j < nvir; j++) {
                ////Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                bse_eigenvector->a[m][vector_offset + count] * psiamp1[p * fermi->nktot * nocc + k_bz * nocc + i] * \
                conj(psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]);
                Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                bse_eigenvector->a[m][vector_offset + count] * psiamp1[p * fermi->nktot * nocc + k_bz * nocc + i] * \
                conj(psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]);
                //fprintf(file.out,"K %3d k1 %3d k_bz %3d m %5d n %3d p %3d o %3d v %3d %14.8f %14.8f %14.8f %10.4f%10.4f%10.4f\n",\
                k,k1,k_bz,m,n,p,i,j,(psiamp1[p * fermi->nktot * nocc + k_bz * nocc + i]).real(),\
                (psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]).real(),\
                (bse_eigenvector->a[m][vector_offset + count]).real(),grid[p].comp1,grid[p].comp2,grid[p].comp3);
                //points1[n].comp1,points1[n].comp2,points1[n].comp3);
                //fprintf(file.out,"   %3d %3d %3d %3d %12.6lf %14.8lf %14.8lf %16.10lf   %10.4lf %10.4lf %10.4lf\n",\
                p,m,i,j,(bse_eigenvector->a[m][count]).real(),(psiamp1[i]).real(),(psiamp2[j]).real(),\
                Psiamp[m * grid_size + p],grid[p].comp1,grid[p].comp2,grid[p].comp3);
                count++;
               }
              }
             //fprintf(file.out,"k k1 %3d %3d m n p %3d %3d %3d Psi %10.4f\n",k,k1,m,n,p,\
             Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]);
        } // close loop on p
             }
            }
           count_k++;
          } // close loop on k1
         } // close loop on k
        ////} // close loop on p
       } // close loop on s

  MPI_Reduce(Psiamp1_buffer, Psiamp1, 2 * buffer1_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //printf("2 %3d %3d\n",job->taskid,k);

  if (job->taskid == 0)
  for (s = 0; s < job->spin_dim; s++) {
    for (p = 0; p < grid_size; p++) {
      for (m = 0; m < num_bse_states; m++) {
        for (n = 0; n < num_points; n++) {
          //Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] = \
          (psiamp1[p * fermi->nktot * nocc + 5 * nocc]).real();
          Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] = \
         (conj(Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]) * \
               Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]).real();
//fprintf(file.out,"%3d %3d %3d %3d %10.4f %10.4f\n",s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p,p,m,n,Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p],(Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]).real());
       } // close loop on n
      } // close loop on m
     } // close loop on p
    } // close loop on s

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

  if (job->taskid == 0)
  write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,num_bse_states,num_points,crystal,atoms,shells,gaussians,job,file);

  free(Psiamp);
  free(Psiamp1);
  free(Psiamp1_buffer);
  free(psiamp1);
  free(psiamp2);
  //free(psiamp1_buffer);
  //free(psiamp2_buffer);
  DestroyComplexMatrix(&bse_eigenvector,job);
  DestroyComplexMatrix(&scf_eigvec_0,job);
  DestroyComplexMatrix(&scf_eigvec_1,job);
  DestroyComplexMatrix(&scf_eigvec_2,job);
  free(Rvec);
  free(grid);
  //free_PAIR_TRAN(&pair_p,job);
  free_k_points(&knet,job);
  MPI_File_close(&fh);

}

void plot_correlated_electron_hole_crystal(double *grid_sum, int bse_state[2],  int *grid_par, VECTOR_DOUBLE *points, int num_points, VECTOR_DOUBLE *points1, FERMI *fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurface of fixed electron-correlated hole for BSE in crystals    *
  // ******************************************************************************************

  int i, j, k, m, n, p, s;
  int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;

  // ******************************************************************************************
  // * Generate knet                                                                          *
  // ******************************************************************************************

  KPOINT_TRAN knet;
  count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  if (job->kss == 0)      fermi->nkunique = knet.nktot;
  else if (job->kss == 1) fermi->nkunique = knet.unique;
  allocate_k_points(&knet,crystal,job,file);
  generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  allocate_fermi(fermi,atoms,job,file);
  fermi->knet = &knet;
  fermi->nktot = knet.nktot;
  fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  int grid_size;
  int radmx = 2;
  int plot_dim = 0;
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

  //for(i = 0; i < number_Rvec; i++) \
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);

  // ******************************************************************************************
  // * Open scf eigenvector file and allocate memory for eigenvectors                         *
  // ******************************************************************************************

  MPI_File fh;
  int dimk;
  int vector_size = nbands * nbfn;
  int block_size = vector_size * sizeof(Complex);
  char buf2[110], filename[15] = "/scf_evec";
  strcpy(buf2,file.bse_eigvec);
  strcat(buf2,filename);
  //printf("scf eigenvectors directory: %s \n",buf2);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  ComplexMatrix *scf_eigvec_0, *scf_eigvec_1, *scf_eigvec_2;
  dimk = job->spin_dim * nbands;

  AllocateComplexMatrix(&scf_eigvec_0,&nbands,&nbfn,job);
  AllocateComplexMatrix(&scf_eigvec_1,&nbands,&nbfn,job);
  AllocateComplexMatrix(&scf_eigvec_2,&nbands,&nbfn,job);

  // ******************************************************************************************
  // * Open bse eigenvector file and allocate memory for eigenvectors                         *
  // ******************************************************************************************

  int ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1) * fermi->nktot;
  int num_bse_states = bse_state[1] - bse_state[0] + 1;
  int spin_grid_size = job->spin_dim * num_bse_states * num_points * grid_size;

  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in BSE_PLOT. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }

  char xc[22] = "/bse_eigenvectors";
  char bufbse[120];
  FILE *bse_evectors;
  long memsize = num_bse_states * ntransitions;
  ComplexMatrix *bse_eigenvector;

  AllocateComplexMatrix(&bse_eigenvector,&num_bse_states,&ntransitions,job);
  ResetComplexMatrix(bse_eigenvector);
  strcpy(bufbse,file.bse_eigvec);
  strcat(bufbse,xc);
  //printf("bse eigenvectors directory: %s \n",bufbse);

  MPI_File gh;
  MPI_File_open(MPI_COMM_WORLD,bufbse,MPI_MODE_RDONLY,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, (bse_state[0] - 1) * ntransitions * sizeof(Complex), MPI_SEEK_SET) ;
  MPI_File_read(gh, &bse_eigenvector->a[0][0], 2 * memsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&gh);

  //ResetComplexMatrix(bse_eigenvector);
  //(bse_eigenvector->a[0][255]).real() = 1.0;
  //(bse_eigenvector->a[0][17*256+255]).real() = 1.0;

  for (j = 0; j < bse_eigenvector->iRows; j++) { \
    for (k = 0; k < bse_eigenvector->iCols; k++) { \
      fprintf(file.out,"%3d %6d kpt %3d occ vir %3d %3d %10.4f %10.4f\n",\
      j,k,k/(nocc*nvir),(k-(k/(nocc*nvir))*(nocc*nvir))/nvir,k%nvir,\
      (bse_eigenvector->a[j][k]).real(),(bse_eigenvector->a[j][k]).imag());
      //bse_eigenvector->a[j][k] /= two * sqrt(bse_eigenvalues[j + bse_state[0] - 1]); 
     } \
    fprintf(file.out,"\n"); \
    }

  int begin_k[job->numtasks], end_k[job->numtasks];
  int total_tasks = fermi->knet->unique;
  mpi_begin_end(begin_k,end_k,total_tasks,job->numtasks,job,file);
  int local_kpts;
  local_kpts = 0;
  for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++)
  local_kpts += (end_k[job->taskid] - begin_k[job->taskid]) * fermi->knet->num[k];
  if (local_kpts == 0) local_kpts = 1;
  //printf("%3d %3d %3d %3d\n",job->taskid,begin_k[job->taskid],end_k[job->taskid],local_kpts);

/*
  int nkunique, ksize, t, tot, q, q1;
  int symmetry_number_of_operators;
  int *ibz, *bz1, *bz2, *opr, *num;
  knet_size(&ksize,fermi->is,crystal);
  if (job->kss == 0)      symmetry_number_of_operators = 1;
  else if (job->kss == 1) symmetry_number_of_operators = symmetry->number_of_operators;
  ibz = (int *) malloc(ksize * symmetry->number_of_operators * sizeof(int));
  if (ibz == NULL) { fprintf(stderr, "ERROR: not enough memory for int ibz\n"); exit(1); }
  bz1 = (int *) malloc(ksize * symmetry->number_of_operators * sizeof(int));
  if (bz1 == NULL) { fprintf(stderr, "ERROR: not enough memory for int bz1\n"); exit(1); }
  bz2 = (int *) malloc(ksize * symmetry->number_of_operators * sizeof(int));
  if (bz2 == NULL) { fprintf(stderr, "ERROR: not enough memory for int bz2\n"); exit(1); }
  opr = (int *) malloc(ksize * symmetry->number_of_operators * sizeof(int));
  if (opr == NULL) { fprintf(stderr, "ERROR: not enough memory for int opr\n"); exit(1); }
  num = (int *) malloc(ksize * symmetry->number_of_operators * sizeof(int));
  if (num == NULL) { fprintf(stderr, "ERROR: not enough memory for int num\n"); exit(1); }

q1 = 0;
int q3, k_count, k_bz1;

  //for (q1 = begin_q[job->taskid]; q1 < end_q[job->taskid]; q1++) {
 k_count = 0;;
  for (q1 = 0; q1 < fermi->nkunique; q1++) {
    for (q3 = 0; q3 < fermi->knet->num[q1]; q3++) {
  q = fermi->knet->ibz[q1];
  KPOINT_TRAN knet_little_q_group;
  SYMMETRY symmetry_little_q_group;
  count_little_k_group_operators(q1,&symmetry_little_q_group,symmetry,crystal,fermi->knet,fermi,job,file);
  allocate_SYMMETRY(&symmetry_little_q_group,job,file);
  generate_little_k_group(q1, &symmetry_little_q_group, fermi, fermi->knet, symmetry, crystal, job, file);
  count_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  allocate_k_points(&knet_little_q_group,crystal,job,file);
  generate_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  //generate_q_lattice(&fermi->knet->oblique[q], &q_G, fermi, G, crystal, job, file);
  rotate_q_crystal1(q1,&nkunique,ibz,bz1,bz2,opr,num,crystal,symmetry,fermi,job,file);
  print_knet(&knet_little_q_group, fermi->is, crystal, job, file);
        if (job->kss == 0) k_bz1 = k_count;                  // use for scf_evec_sym
        if (job->kss == 1) k_bz1 = fermi->knet->bz[k_count]; // use for scf_evec_no_sym
if (job->taskid == 0) printf("%3d %3d %3d  %3d\n",q1,q3,k_count,k_bz1);
 k_count++;
 }
 }

MPI_Finalize();
exit(0);
*/

  // ******************************************************************************************
  // Calculate transition density on grid                                                     *
  // ******************************************************************************************
 
  int count, offset;
  //int q, q1, q3;
  int count_bz, k_bz, kq_bz, fbz1, fbz2;
  Complex *psiamp1_buffer, *psiamp2_buffer;
  Complex *Psiamp1, *Psiamp1_buffer, *psiamp1, *psiamp2;
  double *Psiamp;
  char buf[50] = "bse_electron_hole_density_crystal_isosurface.xsf";
  VECTOR_DOUBLE Rvec1[radmx];

  //printf("%3d %3d %3d %3d %3d %3d %3d\n",num_bse_states,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],ntransitions);
  Rvec1[0].comp1 = 0.0;
  Rvec1[0].comp2 = 0.0;
  Rvec1[0].comp3 = 0.0;

  Psiamp = (double *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(double));
  if (Psiamp == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp \n");
  MPI_Finalize(); exit(1); }

  Psiamp1 = (Complex *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(Complex));
  if (Psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp1 \n");
  MPI_Finalize(); exit(1); }

  Psiamp1_buffer = (Complex *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(Complex));
  if (Psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp1 = (Complex *) calloc(job->spin_dim * grid_size * local_kpts * nocc, sizeof(Complex));
  if (psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp2 = (Complex *) calloc(job->spin_dim * num_points * local_kpts * nvir, sizeof(Complex));
  if (psiamp2 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2 \n");
  MPI_Finalize(); exit(1); }

  int k1, count_k, local_k, vector_offset;
  int buffer1_size = job->spin_dim * num_bse_states * num_points * grid_size;

  ResetComplexArray(Psiamp1,&spin_grid_size);
  ResetComplexArray(Psiamp1_buffer,&spin_grid_size);

  offset = (fermi->bands[0] - 1) * nbfn;
  //printf("offset %3d %3d\n",offset,job->spin_dim);
  for (s = 0; s < job->spin_dim; s++) {
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        //printf("k %3d count_k %3d k_bz %3d   %3d %3d %3d\n",k,count_k,k_bz,\
        fermi->knet->oblique[k_bz].comp1,fermi->knet->oblique[k_bz].comp2,fermi->knet->oblique[k_bz].comp3);
        rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_2->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
        for (n = 0; n < num_points; n++) {
          int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
          if (k_bz == 1) {
          for (int ii = 0; ii < nbands; ii++) { for (int jj = 0; jj < dim1; jj++) {
          fprintf(file.out,"%3d %3d %3d %10.4f %10.4f\n",k,ii,jj,(scf_eigvec_2->a[ii][jj]).real(),(scf_eigvec_2->a[ii][jj]).imag());
          }fprintf(file.out,"\n");}}
          //wavefunction_gridpoint_crystal(&psiamp2[n * fermi->nktot * nvir + k_bz * nvir],&points1[n],1,Rvec1,k_bz,\
          fermi->knet,&scf_eigvec_2->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
          ////wavefunction_gridpoint_crystal(&psiamp2[n * local_kpts * nvir + local_k * nvir],&points1[n],1,Rvec1,k_bz,\
          fermi->knet,&scf_eigvec_2->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
          // particles
          wavefunction_gridpoint_crystal(&psiamp2[n * local_kpts * nvir + local_k * nvir],&points1[n],number_Rvec,Rvec,k_bz,\
          fermi->knet,&scf_eigvec_2->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
          //for (j = 0; j < nvir; j++) \
          fprintf(file.out,"k %3d k1 %3d k_bz %3d n %3d vir %3d %14.8lf %14.8lf %10.4lf %10.4lf %10.4lf\n",\
          k,k1,k_bz,n,j,(psiamp2[n * local_kpts * nvir + local_k * nvir + j]).real(),\
          (psiamp2[n * local_kpts * nvir + local_k * nvir + j]).imag(),points1[n].comp1,points1[n].comp2,points1[n].comp3);
         }
        fprintf(file.out,"\n");
        local_k++;
        count_k++;
       } // close loop on k1
      } // close loop on k
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        //printf("k %3d count_k %3d k_bz %3d   %3d %3d %3d\n",k,count_k,k_bz,\
        fermi->knet->oblique[k_bz].comp1,fermi->knet->oblique[k_bz].comp2,fermi->knet->oblique[k_bz].comp3);
        rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_1->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
        //vector_offset = k_bz * nocc * nvir;
        for (p = 0; p < grid_size; p++) {
          //wavefunction_gridpoint_crystal(&psiamp1[p*fermi->nktot * nocc + k_bz * nocc],&grid[p],number_Rvec,Rvec,k_bz,fermi->knet,\
          &scf_eigvec_1->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
          // holes
          wavefunction_gridpoint_crystal(&psiamp1[p * local_kpts * nocc + local_k * nocc],&grid[p],number_Rvec,Rvec,k_bz,\
          fermi->knet,&scf_eigvec_1->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
         }
        local_k++;
        count_k++;
       }
      printf("c %3d %3d\n",job->taskid,k);
      }
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        vector_offset = k_bz * nocc * nvir;
        for (m = 0; m < num_bse_states; m++) {
          for (n = 0; n < num_points; n++) {
            for (p = 0; p < grid_size; p++) {
              count = 0;
              for (i = 0; i < nocc; i++) {
                for (j = 0; j < nvir; j++) {
                  //Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                  bse_eigenvector->a[m][vector_offset + count] * psiamp1[p * fermi->nktot * nocc + k_bz * nocc + i] * \
                  conj(psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]);
                  Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                  bse_eigenvector->a[m][vector_offset + count] * conj(psiamp1[p * local_kpts * nocc + local_k * nocc + i]) * \
                  psiamp2[n * local_kpts * nvir + local_k * nvir + j];
                  //bse_eigenvector->a[m][vector_offset + count] * psiamp1[p * local_kpts * nocc + local_k * nocc + i] * \
                  conj(psiamp2[n * local_kpts * nvir + local_k * nvir + j]);
                  //Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                  bse_eigenvector->a[m][vector_offset + count] * \
                  conj(psiamp2[n * local_kpts * nvir + local_k * nvir + j]);
                  //fprintf(file.out,"K %3d k1 %3d k_bz %3d m %5d n %3d p %3d o %3d v %3d%14.8f %14.8f %14.8f %10.4f%10.4f%10.4f\n",\
                  k,k1,k_bz,m,n,p,i,j,(psiamp1[p * local_kpts * nocc + local_k * nocc + i]).real(),\
                  (psiamp2[n * local_kpts * nvir + local * nvir + j]).real(),\
                  (bse_eigenvector->a[m][vector_offset + count]).real(),grid[p].comp1,grid[p].comp2,grid[p].comp3);
                  count++;
                 }
                }
               } // close loop on p
              }
             }
            local_k++;
            count_k++;
           } // close loop on k1
          } // close loop on k
         } // close loop on s
    MPI_Reduce(Psiamp1_buffer, Psiamp1, 2 * buffer1_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (job->taskid == 0)
  for (s = 0; s < job->spin_dim; s++) {
    for (p = 0; p < grid_size; p++) {
      for (m = 0; m < num_bse_states; m++) {
        for (n = 0; n < num_points; n++) {
          Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] = \
         (conj(Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]) * \
               Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]).real();
               //(Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]).real();
         } // close loop on n
        } // close loop on m
       } // close loop on p
      } // close loop on s

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

  if (job->taskid == 0)
  write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,num_bse_states,num_points,crystal,atoms,shells,gaussians,job,file);

  free(Psiamp);
  free(Psiamp1);
  free(Psiamp1_buffer);
  free(psiamp1);
  free(psiamp2);
  DestroyComplexMatrix(&bse_eigenvector,job);
  DestroyComplexMatrix(&scf_eigvec_0,job);
  DestroyComplexMatrix(&scf_eigvec_1,job);
  DestroyComplexMatrix(&scf_eigvec_2,job);
  free(Rvec);
  free(grid);
  free_k_points(&knet,job);
  MPI_File_close(&fh);

}

void plot_electron_hole_transition_density_crystal(double *grid_sum, int bse_state[2],  int *grid_par, VECTOR_DOUBLE *points, int num_points, VECTOR_DOUBLE *points1, FERMI *fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates isosurface of fixed electron-correlated hole for BSE in crystals    *
  // ******************************************************************************************

  int i, j, k, m, n, p, s;
  int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;

  // ******************************************************************************************
  // * Generate knet                                                                          *
  // ******************************************************************************************

  KPOINT_TRAN knet;
  count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  if (job->kss == 0)      fermi->nkunique = knet.nktot;
  else if (job->kss == 1) fermi->nkunique = knet.unique;
  allocate_k_points(&knet,crystal,job,file);
  generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  allocate_fermi(fermi,atoms,job,file);
  fermi->knet = &knet;
  fermi->nktot = knet.nktot;
  fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

  // ******************************************************************************************
  // * Generate real space grid                                                               *
  // ******************************************************************************************

  int grid_size;
  int radmx = 3;
  int plot_dim = 0;
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

  //for(i = 0; i < number_Rvec; i++) \
  fprintf(file.out,"%d  %10.4lf %10.4lf %10.4lf\n",i,Rvec[i].comp1, Rvec[i].comp2, Rvec[i].comp3);

  // ******************************************************************************************
  // * Open scf eigenvector file and allocate memory for eigenvectors                         *
  // ******************************************************************************************

  MPI_File fh;
  int dimk;
  int vector_size = nbands * nbfn;
  int block_size = vector_size * sizeof(Complex);
  char buf4[110], filename[15] = "/scf_evec";
  strcpy(buf4,file.bse_eigvec);
  strcat(buf4,filename);
  //printf("scf eigenvectors directory: %s \n",buf4);
  MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  ComplexMatrix *scf_eigvec_0, *scf_eigvec_1, *scf_eigvec_2;
  dimk = job->spin_dim * nbands;

  AllocateComplexMatrix(&scf_eigvec_0,&nbands,&nbfn,job);
  AllocateComplexMatrix(&scf_eigvec_1,&nbands,&nbfn,job);
  AllocateComplexMatrix(&scf_eigvec_2,&nbands,&nbfn,job);

  // ******************************************************************************************
  // * Open bse eigenvector file and allocate memory for eigenvectors                         *
  // ******************************************************************************************

  int ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1) * fermi->nktot;
  int num_bse_states = bse_state[1] - bse_state[0] + 1;
  int spin_grid_size = job->spin_dim * num_bse_states * num_points * grid_size;

  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in BSE_PLOT. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }

  char xc[22] = "/bse_eigenvectors";
  char bufbse[120];
  FILE *bse_evectors;
  long memsize = num_bse_states * ntransitions;
  ComplexMatrix *bse_eigenvector;

  AllocateComplexMatrix(&bse_eigenvector,&num_bse_states,&ntransitions,job);
  ResetComplexMatrix(bse_eigenvector);
  strcpy(bufbse,file.bse_eigvec);
  strcat(bufbse,xc);

  MPI_File gh;
  MPI_File_open(MPI_COMM_WORLD,bufbse,MPI_MODE_RDONLY,MPI_INFO_NULL,&gh) ;
  MPI_File_seek(gh, (bse_state[0] - 1) * ntransitions * sizeof(Complex), MPI_SEEK_SET) ;
  MPI_File_read(gh, &bse_eigenvector->a[0][0], 2 * memsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&gh);

  //for (j = 0; j < bse_eigenvector->iRows; j++) { \
    for (k = 0; k < bse_eigenvector->iCols; k++) { \
      fprintf(file.out,"%3d %3d %10.4f %10.4f\n",j,k,(bse_eigenvector->a[j][k]).real(),(bse_eigenvector->a[j][k]).imag()); \
      bse_eigenvector->a[j][k] /= two * sqrt(bse_eigenvalues[j + bse_state[0] - 1]);  \
     } \
    fprintf(file.out,"\n"); \
    }

  int begin_k[job->numtasks], end_k[job->numtasks];
  int total_tasks = fermi->knet->unique;
  mpi_begin_end(begin_k,end_k,total_tasks,job->numtasks,job,file);
  int local_kpts;
  local_kpts = 0;
  for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++)
  local_kpts += (end_k[job->taskid] - begin_k[job->taskid]) * fermi->knet->num[k];
  if (local_kpts == 0) local_kpts = 1;
local_kpts = 1;
  //printf("%3d %3d %3d %3d\n",job->taskid,begin_k[job->taskid],end_k[job->taskid],local_kpts);

  // ******************************************************************************************
  // Calculate transition density on grid                                                     *
  // ******************************************************************************************
 
  char yy[24] = "crystal_isosurface.xsf";
  char buf1[29] = "hole_";
  char buf2[33] = "electron_";
  char buf3[43] = "transition_density_";
  strcat(buf1,yy);
  strcat(buf2,yy);
  strcat(buf3,yy);

  int count, offset;
  int q, q1, q3;
  int count_bz, k_bz, kq_bz, fbz1, fbz2;
  Complex *psiamp1_buffer, *psiamp2_buffer;
  Complex *Psiamp1, *Psiamp1_buffer, *psiamp1, *psiamp2;
  double *Psiamp;
  //char buf[50] = "bse_electron_hole_transition_density_crystal.xsf";
  VECTOR_DOUBLE Rvec1[radmx];

  //printf("%3d %3d %3d %3d %3d %3d %3d\n",num_bse_states,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],ntransitions);
  Rvec1[0].comp1 = 0.0;
  Rvec1[0].comp2 = 0.0;
  Rvec1[0].comp3 = 0.0;

  Psiamp = (double *) calloc(job->spin_dim * grid_size * local_kpts * nvir, sizeof(double));
  //Psiamp = (double *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(double));
  if (Psiamp == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp \n");
  MPI_Finalize(); exit(1); }

  Psiamp1 = (Complex *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(Complex));
  if (Psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp1 \n");
  MPI_Finalize(); exit(1); }

  Psiamp1_buffer = (Complex *) calloc(job->spin_dim * num_bse_states * num_points * grid_size, sizeof(Complex));
  if (Psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for Psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp1 = (Complex *) calloc(job->spin_dim * grid_size * local_kpts * nocc, sizeof(Complex));
  if (psiamp1 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp1 \n");
  MPI_Finalize(); exit(1); }

  psiamp2 = (Complex *) calloc(job->spin_dim * grid_size * local_kpts * nvir, sizeof(Complex));
  if (psiamp2 == NULL) { fprintf(file.out, "ERROR: There is not enough memory for psiamp2 \n");
  MPI_Finalize(); exit(1); }

  int ibz, k1, count_k, local_k, vector_offset;
  int buffer1_size = job->spin_dim * num_bse_states * num_points * grid_size;

  ResetComplexArray(Psiamp1,&spin_grid_size);
  ResetComplexArray(Psiamp1_buffer,&spin_grid_size);

  offset = (fermi->bands[0] - 1) * nbfn;
  //printf("offset %3d %3d\n",offset,job->spin_dim);
  for (s = 0; s < job->spin_dim; s++) {
    count_k = 0;
    //for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    //for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
    k = 3;
    count_k = 13;
    MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
    MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    k_bz = fermi->knet->bz[count_k];
    printf("k %3d count_k %3d k_bz %3d   %3d %3d %3d\n",k,count_k,k_bz,\
    fermi->knet->oblique[k_bz].comp1,fermi->knet->oblique[k_bz].comp2,fermi->knet->oblique[k_bz].comp3);
    rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_1->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
    //ResetComplexMatrix(scf_eigvec_1);
    //(scf_eigvec_1->a[0][11]).imag() = 1.0;
    for (p = 0; p < grid_size; p++) {
      wavefunction_gridpoint_crystal(&psiamp1[p * local_kpts * nocc + local_k * nocc],&grid[p],number_Rvec,Rvec,k_bz,\
      fermi->knet,&scf_eigvec_1->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
     }
    local_k++;
    count_k++;
 //} // close loop on k
  } // close loop on s

  local_k = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < nocc; i++) {
      for (p = 0; p < grid_size; p++) {
        Psiamp[s * nocc * grid_size + i * grid_size + p] = (psiamp1[p * local_kpts * nocc + local_k * nocc + i]).real();
      } // close loop on p
     } // close loop on i
    } // close loop on s


  if (job->taskid == 0)
  write_isosurface_xsf(buf1,Psiamp,grid_size,grid_par,points,nocc,1,crystal,atoms,shells,gaussians,job,file);

  offset = (fermi->bands[0] - 1) * nbfn;
  //printf("offset %3d %3d\n",offset,job->spin_dim);
  for (s = 0; s < job->spin_dim; s++) {
    count_k = 0;
    //for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    //for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      k = 3;
      count_k = 13;
      MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      k_bz = fermi->knet->bz[count_k];
      printf("k %3d count_k %3d k_bz %3d   %3d %3d %3d\n",k,count_k,k_bz,\
      fermi->knet->oblique[k_bz].comp1,fermi->knet->oblique[k_bz].comp2,fermi->knet->oblique[k_bz].comp3);
      rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_2->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
      //ResetComplexMatrix(scf_eigvec_2);
      //(scf_eigvec_2->a[nocc][11]).real() = 0.5;
      for (p = 0; p < grid_size; p++) {
        wavefunction_gridpoint_crystal(&psiamp2[p * local_kpts * nvir + local_k * nvir],&grid[p],number_Rvec,Rvec,k_bz,\
        fermi->knet,&scf_eigvec_2->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
       }
      local_k++;
      count_k++;
   //} // close loop on k
    } // close loop on s

  local_k = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (j = 0; j < nvir; j++) {
      for (p = 0; p < grid_size; p++) {
        Psiamp[s * nvir * grid_size + j * grid_size + p] = (psiamp2[p * local_kpts * nvir + local_k * nvir + j]).imag();
      } // close loop on m
     } // close loop on p
    } // close loop on s
  if (job->taskid == 0)
  write_isosurface_xsf(buf2,Psiamp,grid_size,grid_par,points,nvir,1,crystal,atoms,shells,gaussians,job,file);

/*
  for (s = 0; s < job->spin_dim; s++) {
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      MPI_File_seek(fh,((s * fermi->knet->unique + k) * nbfn * nbfn + offset) * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_eigvec_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        rotate_psi(&scf_eigvec_0->a[0][0],&scf_eigvec_1->a[0][0],nbands,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
        vector_offset = k_bz * nocc * nvir;
        for (p = 0; p < grid_size; p++) {
          //wavefunction_gridpoint_crystal(&psiamp1[p*fermi->nktot * nocc + k_bz * nocc],&grid[p],number_Rvec,Rvec,k_bz,fermi->knet,\
          &scf_eigvec_1->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
          wavefunction_gridpoint_crystal(&psiamp1[p * local_kpts * nocc + local_k * nocc],&grid[p],number_Rvec,Rvec,k_bz,\
          fermi->knet,&scf_eigvec_1->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
         }
        local_k++;
        count_k++;
       }
      printf("c %3d %3d\n",job->taskid,k);
      }
         } // close loop on s

  for (s = 0; s < job->spin_dim; s++) {
    count_k = 0;
    for (k = 0; k < begin_k[job->taskid]; k++) { for (k1 = 0; k1 < fermi->knet->num[k]; k1++) { count_k++; } }
    local_k = 0;
    for (k = begin_k[job->taskid]; k < end_k[job->taskid]; k++) {
      for (k1 = 0; k1 < fermi->knet->num[k]; k1++) {
        k_bz = fermi->knet->bz[count_k];
        vector_offset = k_bz * nocc * nvir;
        for (m = 0; m < num_bse_states; m++) {
          for (n = 0; n < num_points; n++) {
            for (p = 0; p < grid_size; p++) {
              count = 0;
              for (i = 0; i < nocc; i++) {
                for (j = 0; j < nvir; j++) {
                  //Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                  bse_eigenvector->a[m][vector_offset + count] * psiamp1[p * fermi->nktot * nocc + k_bz * nocc + i] * \
                  conj(psiamp2[n * fermi->nktot * nvir + k_bz * nvir + j]);
                  Psiamp1_buffer[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] += \
                  bse_eigenvector->a[m][vector_offset + count] * psiamp1[p * local_kpts * nocc + local_k * nocc + i] * \
                  conj(psiamp2[n * local_kpts * nvir + local_k * nvir + j]);
                  //fprintf(file.out,"K %3d k1 %3d k_bz %3d m %5d n %3d p %3d o %3d v %3d%14.8f %14.8f %14.8f %10.4f%10.4f%10.4f\n",\
                  k,k1,k_bz,m,n,p,i,j,(psiamp1[p * local_kpts * nocc + local_k * nocc + i]).real(),\
                  (psiamp2[n * local_kpts * nvir + local * nvir + j]).real(),\
                  (bse_eigenvector->a[m][vector_offset + count]).real(),grid[p].comp1,grid[p].comp2,grid[p].comp3);
                  count++;
                 }
                }
               } // close loop on p
              }
             }
            local_k++;
            count_k++;
           } // close loop on k1
          } // close loop on k
         } // close loop on s
    MPI_Reduce(Psiamp1_buffer, Psiamp1, 2 * buffer1_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (job->taskid == 0)
  for (s = 0; s < job->spin_dim; s++) {
    for (p = 0; p < grid_size; p++) {
      for (m = 0; m < num_bse_states; m++) {
        for (n = 0; n < num_points; n++) {
          Psiamp[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p] = \
         (conj(Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]) * \
               Psiamp1[s * num_bse_states * num_points * grid_size + m * num_points * grid_size + n * grid_size + p]).real();
         } // close loop on n
        } // close loop on m
       } // close loop on p
      } // close loop on s
*/

  // ******************************************************************************************
  // * Generate plot                                                                          *
  // ******************************************************************************************

  //if (job->taskid == 0)
  //write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,nvir,1,crystal,atoms,shells,gaussians,job,file);
  //write_isosurface_xsf(buf3, Psiamp, grid_size, grid_par, points, nocc, nvir, crystal, atoms, shells, gaussians, job,file);
  //write_isosurface_xsf(buf,Psiamp,grid_size,grid_par,points,num_bse_states,num_points,crystal,atoms,shells,gaussians,job,file);

  free(Psiamp);
  free(Psiamp1);
  free(Psiamp1_buffer);
  free(psiamp1);
  free(psiamp2);
  DestroyComplexMatrix(&bse_eigenvector,job);
  DestroyComplexMatrix(&scf_eigvec_0,job);
  DestroyComplexMatrix(&scf_eigvec_1,job);
  DestroyComplexMatrix(&scf_eigvec_2,job);
  free(Rvec);
  free(grid);
  free_k_points(&knet,job);
  MPI_File_close(&fh);

}

void plot_electron_hole_transition_density(int grid_par[3], VECTOR_DOUBLE *points, FERMI* fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  int i, j, m, n, p, s;
  int dimk = 1;
  int nbands = fermi->bands[1] - fermi->bands[0] + 1;
  int nvir = fermi->bands[1] - fermi->homo[0];
  int nocc = fermi->homo[0] - fermi->bands[0] + 1;
  int ntransitions;
  ntransitions = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  PAIR_TRAN pair_p;

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity >= 1)
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
  // * Read and broadcast scf eigenvectors                                                    *
  // ******************************************************************************************

  int vector_size;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim2 = dim1 * dim1;
  char zz2[24] = "scf_evectors";
  size_t result;
  FILE *scf_evectors;
  DoubleMatrix *eigvec;
  ComplexMatrix *scf_eigenvectors;
  dimk = job->spin_dim * nbands;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&scf_eigenvectors,&dimk,&dim1,job);
  AllocateDoubleMatrix(&eigvec,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);

  if (job->taskid == 0) {
  scf_evectors = fopen(zz2, "rb");
  fseek(scf_evectors, dim1 * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
  result = fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
  fclose(scf_evectors);

  for (i = 0; i < nbands; i++) {
    for (j = 0; j < dim1; j++) {
      eigvec->a[i][j] = (scf_eigenvectors->a[i][j]).real();
      //fprintf(file.out,"eigvec %3d %3d %10.4lf\n",i,j,eigvec->a[i][j]);
     }
    }
   }
  MPI_Bcast(&eigvec->a[0][0],nbands * dim1,MPI_DOUBLE,0,MPI_COMM_WORLD);

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

  printf("%3d %3d %3d %3d %3d %3d %3d\n",grid_size,nocc,nvir,fermi->bands[0],fermi->bands[1],fermi->homo[0],ntransitions);
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

  fprintf(file.out,"hole grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(&psiamp1[p * nocc],&grid[p],Rvec1,&eigvec->a[0][0],nocc,nbands,atoms,shells,gaussians,job,file);
   }

  fprintf(file.out,"particle grids\n");
  for (p = 0; p < grid_size; p++) {
    wavefunction_gridpoint(&psiamp2[p * nvir],&grid[p],Rvec2,&eigvec->a[nocc][0],nvir,nbands,atoms,shells,gaussians,job,file);
   }

  for (p = 0; p < grid_size; p++) {
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < nocc; i++) {
        Psiamp[s * nocc * grid_size + i * grid_size + p] = psiamp1[s * grid_size * nocc + p * nocc + i];
       }
      }
     }

  write_isosurface_xsf(buf1, Psiamp, grid_size, grid_par, points, nocc, 1, crystal, atoms, shells, gaussians, job,file);

  for (p = 0; p < grid_size; p++) {
    for (s = 0; s < job->spin_dim; s++) {
      for (i = 0; i < nvir; i++) {
        Psiamp[s * nvir * grid_size + i * grid_size + p] = psiamp2[s * grid_size * nvir + p * nvir + i];
       }
      }
     }

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

  write_isosurface_xsf(buf3, Psiamp, grid_size, grid_par, points, nocc, nvir, crystal, atoms, shells, gaussians, job,file);

  free(grid);
  free(Rvec);
  free(Psiamp);
  free(psiamp1);
  free(psiamp2);
  free(psiamp1_buffer);
  free(psiamp2_buffer);
  free_PAIR_TRAN(&pair_p,job);
  DestroyDoubleMatrix(&eigvec,job);
  DestroyComplexMatrix(&scf_eigenvectors,job);
//MPI_Finalize();
//exit(0);
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

void wavefunction_gridpoint_crystal(Complex *psiamp, VECTOR_DOUBLE *rpoint, int number_Rvec, VECTOR_DOUBLE *Rvec, int k_bz, KPOINT_TRAN *knet, Complex *eigvec, int nstates, int nbands, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  // get psi for nbands and spins on rpoint

  int m, s;
  int index_i, sheli, shelposi, bfposi, gausposi;
  int i1, i2, j2, i3, j3, i4, j4;
  int ip, jp;
  int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  double psiexp, rsqrd, kdotr;
  double prefac[7];
  double r[3];
  Complex temp;
  VECTOR_DOUBLE *p_Rvec;

  //for (i1 = 0; i1 < job->spin_dim * nstates; i1++)
  //psiamp[i1] = Complex(k_zero, k_zero);

  p_Rvec = Rvec;
  for (jp = 0; jp < number_Rvec; jp++) {
    kdotr = double_vec_dot(&knet->cart[k_bz],p_Rvec);
    temp = Complex(cos(kdotr), sin(kdotr));
    //if (jp > 200 && jp < 220)  \
    fprintf(file.out,"k_bz %3d k %10.4f %10.4f %10.4f R %10.4f %10.4f %10.4f  %3d kdotr %10.4f  %10.4f %10.4f\n",\
    k_bz,knet->cart[k_bz].comp1,knet->cart[k_bz].comp2,knet->cart[k_bz].comp3,p_Rvec->comp1,p_Rvec->comp2,p_Rvec->comp3,\
    jp,kdotr,cos(kdotr),sin(kdotr));
    //if (jp == 13) fprintf(file.out,"k_bz %3d R %10.4f %10.4f %10.4f  %3d kdotr %10.4f  %10.4f %10.4f\n",\
    k_bz,p_Rvec->comp1,p_Rvec->comp2,p_Rvec->comp3,jp,kdotr,cos(kdotr),sin(kdotr));
    for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
    //r[0] = atoms->cell_vector[ip].comp1 + Rvec->comp1;
    //r[1] = atoms->cell_vector[ip].comp2 + Rvec->comp2;
    //r[2] = atoms->cell_vector[ip].comp3 + Rvec->comp3;
    r[0] = atoms->cell_vector[ip].comp1 + p_Rvec->comp1;
    r[1] = atoms->cell_vector[ip].comp2 + p_Rvec->comp2;
    r[2] = atoms->cell_vector[ip].comp3 + p_Rvec->comp3;
    rsqrd = (r[0] - rpoint->comp1) * (r[0] - rpoint->comp1) + \
            (r[1] - rpoint->comp2) * (r[1] - rpoint->comp2) + \
            (r[2] - rpoint->comp3) * (r[2] - rpoint->comp3);
    //fprintf(file.out,"ip %3d jp %3d rsq %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n", \
    ip,jp,rsqrd,r[0],r[1],r[2],rpoint->comp1,rpoint->comp2,rpoint->comp3);
    if (rsqrd > 200.0) continue;
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
    //for (s = 0; s < job->spin_dim; s++) {
      //for (i = 0; i < atoms->number_of_sh_bfns_in_unit_cell; i++) {
      for (m = 0; m < nstates; m++) {
        for (i3 = 0; i3 < sheli; i3++) {
         psiamp[m] += eigvec[m * dim1 + bfposi + i3] * psiexp * temp * prefac[i3];
         ////psiamp[s * nstates + m] += eigvec[s * nbands * dim1 + m * dim1 + bfposi + i3] * psiexp * prefac[i3];
         //fprintf(file.out,"psiamp %3d %3d %3d %3d %3d %3d %10.4lf %10.4lf  %10.4lf %10.4lf\n",ip,bfposi,index_i,s,m,i3, \
         (eigvec[s * nbands * dim1 + m * dim1 + bfposi + i3]).real(), psiexp, prefac[i3],(psiamp[s * nstates + m]).real());
        } // close loop over i3
       } // close loop on m
      //} // close loop on s
     gausposi += shells->ng_sh[index_i];
     bfposi += shells->type_sh[index_i];
    } // close loop on index_i
   } // close loop on ip
  p_Rvec++;
 } // close loop on jp

    if (job->taskid == 0 && job->verbosity > 1) {
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < nstates; m++) {
          fprintf(file.out,"%3d %3d %10.4lf %10.4lf\n",s,m,(psiamp[s * nstates + m]).real(),(psiamp[s * nstates + m]).imag());
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
  //file_prt = fopen("isosurface.xsf", "w");
  if (file_prt == NULL) { fprintf(file.out, "CANNOT OPEN FILE FOR PLOTTING IN write_xsf_isosurface_file\n"); exit(1); }
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

void calc_grid(int *grid_par, VECTOR_DOUBLE *grid, int grid_size, int *number_of_grid_points, VECTOR_DOUBLE *points, JOB_PARAM *job, FILES file)

{

  int i, j, k;
  VECTOR_DOUBLE vec[3];
  double incrementX[3], incrementY[3], incrementZ[3];
  unsigned char plot_dim = 3; //DEFAULT is a 3D plot

  // define three vectors for ploting
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
  for (i = 0; i < grid_par[0]; i++) {
    for (j = 0; j < grid_par[1]; j++) {
      for (k = 0; k < grid_par[2]; k++) {
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


