#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include <xc.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "myconstants.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "LIMITS.h"
#include "SYMMETRY.h"
#include "TOOLS.h"
#include "PARALLEL.h"
#include "MATRIX_UTIL.h"
#include "RECURSION.h"
#include "E_COEFFICIENTS.h"
//#include "INTEGRALS.h"
#include "LEBEDEV_LAIKOV.h"
#include "SCF_ATOM.h"
#include "DFT_ATOM.h"
#include "DFT.h"

using namespace std;

void count_dft_savin_grid(int *ngrid_points, PAIR_TRAN *pair_p, ATOM *atoms, CRYSTAL *crystal, REAL_LATTICE *R, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, m, p, q, l;
int g1, g2, p1, q1, q2;
int L = job->xc_lmx + 1, nradial = job->xc_rad, nsphere;
int L_range[10];
double r_range[10];
double r,rm, x3, cell_vol, w_tot;
double r1, r2, s1, s2;
double expnt1, expnt2;
double weights[(L * 2 + 3) * (L * 2 + 3) / 3];
double theta[(L * 2 + 3) * (L * 2 + 3) / 3];
double phi[(L * 2 + 3) * (L * 2 + 3) / 3];
double grid_X[(nradial + 1) * ((L * 2 + 3) * (L * 2 + 3) / 3)];
double grid_Y[(nradial + 1) * ((L * 2 + 3) * (L * 2 + 3) / 3)];
double grid_Z[(nradial + 1) * ((L * 2 + 3) * (L * 2 + 3) / 3)];
double y[nradial], radial_weight;
double w_grid;
VECTOR_DOUBLE Rvec_tmp[3], grid;


 *ngrid_points = 0;
  cell_vol = k_zero;
  for (m = 0; m < atoms->number_of_atoms_in_unit_cell; m++) {
    l = 0;
    dft_grid_parameters(L_range,r_range,job,file);
    Gauss_Chebyshev_weights_abscissa(nradial,y,&radial_weight,job,file);
    Lebedev_grid(theta, phi, grid_X, grid_Y, grid_Z, weights, &nsphere, L_range[l], file);
    rm = k_one;
    for (i = 0; i < nradial; i++) {
      r = rm * (k_one + y[i]) / (k_one - y[i]);
      x3 = four * pi * two * r * r * sqrt(r * rm) / (k_one - y[i]) * radial_weight;
      if (r > r_range[l]) {
      Lebedev_grid(theta, phi, grid_X, grid_Y, grid_Z, weights, &nsphere, L_range[l], file);
      l++;
     }
      for (j = 0; j < nsphere; j++) {
        grid.comp1 = r * grid_X[j];
        grid.comp2 = r * grid_Y[j];
        grid.comp3 = r * grid_Z[j];
        w_grid = k_zero;
        w_tot = k_zero;
        for (p = 0; p < pair_p->nump; p++) {
          p1  = pair_p->posn[p];
          for (q = 0; q < pair_p->numb[p]; q++) {
            q1 = pair_p->cell1[p1 + q];
            q2 = pair_p->cell2[p1 + q];
            g1 = pair_p->latt2[p1 + q];
            if (pair_p->cell1[p1 + q] != m) continue;
            if (q1 == q2 && g1 == 0) {
              r1 = sqrt(double_vec_dot(&grid,&grid));   
              //expnt1 = 1.0 * (double)atoms->atomic_number[q1];
              expnt1 = 5.0 ;
              s1 = exp(-expnt1 * r1);
             }
            else {
              Rvec_tmp[1].comp1 = atoms->cell_vector[q1].comp1 + grid.comp1 - atoms->cell_vector[q2].comp1 - R->vec_ai[g1].comp1;
              Rvec_tmp[1].comp2 = atoms->cell_vector[q1].comp2 + grid.comp2 - atoms->cell_vector[q2].comp2 - R->vec_ai[g1].comp2;
              Rvec_tmp[1].comp3 = atoms->cell_vector[q1].comp3 + grid.comp3 - atoms->cell_vector[q2].comp3 - R->vec_ai[g1].comp3;
              r2 = sqrt(double_vec_dot(&Rvec_tmp[1],&Rvec_tmp[1]));   
              //expnt2 = 1.0 * (double)atoms->atomic_number[q2];
              expnt2 = 5.0 ;
              w_tot += exp(-expnt2 * r2);
             }
            } // loop over q
           } // loop over p
            if (fabs(s1) < 1e-20 && fabs(w_tot) < 1e-20) w_grid = k_zero;
            else w_grid = s1 / (s1 + w_tot);
            cell_vol += x3 * w_grid * weights[j];
            if (w_grid > k_tol / 100.0) {
            (*ngrid_points)++;
           }
          } // loop over grid points
         } // loop over grid points
        } // loop over m

       if (job->taskid == 0 && crystal->type[0] != 'M') {
       if (fabs(k_one - cell_vol/crystal->primitive_cell_volume) > 0.0000001)
       fprintf(file.out,"WARNING:");
       if (fabs(k_one - cell_vol/crystal->primitive_cell_volume) > 0.0000001 || job->verbosity > 1)
       fprintf(file.out,"Primitive Cell Volume Integration in DFT Grid  Radial Points %5d Lebedev Grid Points %5d"
       "\nPrimitive Cell Volume by integration %10.7e Actual Volume %10.7e Ratio %10.7e\n",nradial,nsphere,cell_vol,\
       crystal->primitive_cell_volume,cell_vol/crystal->primitive_cell_volume);
      }

}

void generate_dft_savin_grid(DFT_GRID *dft_grid, PAIR_TRAN *pair_p, ATOM *atoms, CRYSTAL *crystal, REAL_LATTICE *R, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, m, p, q, l;
int g1, g2, p1, q1, q2;
int npoints;
int L = job->xc_lmx + 1, nradial = job->xc_rad, nsphere;
int L_range[10];
double r_range[10];
double r,rm, x3, cell_vol, w_tot;
double r1, r2, s1, s2;
double expnt1, expnt2;
double weights[(L * 2 + 3) * (L * 2 + 3) / 3];
double theta[(L * 2 + 3) * (L * 2 + 3) / 3];
double phi[(L * 2 + 3) * (L * 2 + 3) / 3];
double grid_X[(nradial + 1) * ((L * 2 + 3) * (L * 2 + 3) / 3)];
double grid_Y[(nradial + 1) * ((L * 2 + 3) * (L * 2 + 3) / 3)];
double grid_Z[(nradial + 1) * ((L * 2 + 3) * (L * 2 + 3) / 3)];
double y[nradial], radial_weight;
double w_grid;
VECTOR_DOUBLE Rvec_tmp[3], grid;

  npoints = 0;
  cell_vol = k_zero;
  for (m = 0; m < atoms->number_of_atoms_in_unit_cell; m++) {
    l = 0;
    dft_grid_parameters(L_range,r_range,job,file);
    Gauss_Chebyshev_weights_abscissa(nradial,y,&radial_weight,job,file);
    Lebedev_grid(theta, phi, grid_X, grid_Y, grid_Z, weights, &nsphere, L_range[l], file);
    rm = k_one;
    dft_grid->numb[m] = 0;
    dft_grid->posn[m] = npoints;
    for (i = 0; i < nradial; i++) {
      r = rm * (k_one + y[i]) / (k_one - y[i]);
      x3 = four * pi * two * r * r * sqrt(r * rm) / (k_one - y[i]) * radial_weight;
      if (r > r_range[l]) {
      Lebedev_grid(theta, phi, grid_X, grid_Y, grid_Z, weights, &nsphere, L_range[l], file);
      l++;
     }
      for (j = 0; j < nsphere; j++) {
        grid.comp1 = r * grid_X[j];
        grid.comp2 = r * grid_Y[j];
        grid.comp3 = r * grid_Z[j];
        w_grid = k_zero;
        w_tot = k_zero;
        for (p = 0; p < pair_p->nump; p++) {
          p1  = pair_p->posn[p];
          for (q = 0; q < pair_p->numb[p]; q++) {
            q1 = pair_p->cell1[p1 + q];
            q2 = pair_p->cell2[p1 + q];
            g1 = pair_p->latt2[p1 + q];
            if (pair_p->cell1[p1 + q] != m) continue;
            if (q1 == q2 && g1 == 0) {
              r1 = sqrt(double_vec_dot(&grid,&grid));   
              //expnt1 = 1.0 * (double)atoms->atomic_number[q1];
              expnt1 = 5.0 ;
              s1 = exp(-expnt1 * r1);
             }
            else {
              Rvec_tmp[1].comp1 = atoms->cell_vector[q1].comp1 + grid.comp1 - atoms->cell_vector[q2].comp1 - R->vec_ai[g1].comp1;
              Rvec_tmp[1].comp2 = atoms->cell_vector[q1].comp2 + grid.comp2 - atoms->cell_vector[q2].comp2 - R->vec_ai[g1].comp2;
              Rvec_tmp[1].comp3 = atoms->cell_vector[q1].comp3 + grid.comp3 - atoms->cell_vector[q2].comp3 - R->vec_ai[g1].comp3;
              r2 = sqrt(double_vec_dot(&Rvec_tmp[1],&Rvec_tmp[1]));   
              //expnt2 = 1.0 * (double)atoms->atomic_number[q2];
              expnt2 = 5.0 ;
              w_tot += exp(-expnt2 * r2);
             }
            } // loop over q
           } // loop over p
            if (fabs(s1) < 1e-20 && fabs(w_tot) < 1e-20) w_grid = k_zero;
            else w_grid = s1 / (s1 + w_tot);
            if (w_grid > k_tol / 100.0) {
            dft_grid->x[npoints].comp1 = grid.comp1;
            dft_grid->x[npoints].comp2 = grid.comp2;
            dft_grid->x[npoints].comp3 = grid.comp3;
            dft_grid->r[npoints]       = r;
            dft_grid->y[npoints]       = y[i];
            dft_grid->weight[npoints]  = x3 * w_grid * weights[j];
           (dft_grid->numb[m])++;
            //fprintf(file.out,"%4d grid %10.4lf %10.4lf %10.4lf r %10.4lf y %10.4lf weight %10.3e %10.3e\n", \
            npoints,grid.comp1,grid.comp2,grid.comp3,r,y[i],w_grid,dft_grid->weight[npoints] );
            cell_vol += x3 * w_grid * weights[j];
            npoints++;
           }
          } // loop over grid points
         } // loop over grid points
        } // loop over m
       dft_grid->npoints = npoints;

       if (job->taskid == 0 && crystal->type[0] != 'M') {
       if (fabs(k_one - cell_vol/crystal->primitive_cell_volume) > 0.0000001)
       fprintf(file.out,"WARNING:");
       if (fabs(k_one - cell_vol/crystal->primitive_cell_volume) > 0.0000001 || job->verbosity > 1)
       fprintf(file.out,"Primitive Cell Volume Integration in DFT Grid  Radial Points %5d Lebedev Grid Points %5d"
       "\nPrimitive Cell Volume by integration %10.7e Actual Volume %10.7e Ratio %10.7e\n",nradial,nsphere,cell_vol,\
       crystal->primitive_cell_volume,cell_vol/crystal->primitive_cell_volume);
      }

}

void generate_orbital_grid(DFT_GRID *dft_grid, double *grid_basis, double *F, TRIPLE_TRAN *Triple, PAIR_TRAN *pair_p, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int gl, i, j, k, j4, n, p, p1, q, q0, q1, q2, g1, g2, r, t, t1;
int bfposi, index_j, shelj, shelposj, gausposj;
int F_off, F_offset, count, count1;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
double xsqrd, ysqrd, zsqrd, rsqrd, alpha_rsqrd;
double time1, time2;
VECTOR_DOUBLE R_vec_tmp[2];

  count = 0;
  for (p = 0; p < Triple->tot; p++) {
    n  = dft_grid->posn[Triple->cell1[p]];
    q0 = Triple->cell1[p];
    q1 = Triple->cell2[p];
    q2 = Triple->cell3[p];
    g1 = Triple->latt2[p];
    g2 = Triple->latt3[p];
    gl = 0;
    //fprintf(file.out,"T %3d n %6d q0 %3d q1 %3d q2 %3d g1 %3d g2 %3d\n",p,n,q0,q1,q2,g1,g2);
    R_vec_tmp[0].comp1 = atoms->cell_vector[q0].comp1 - atoms->cell_vector[q2].comp1 - R->vec_ai[g2].comp1;
    R_vec_tmp[0].comp2 = atoms->cell_vector[q0].comp2 - atoms->cell_vector[q2].comp2 - R->vec_ai[g2].comp2;
    R_vec_tmp[0].comp3 = atoms->cell_vector[q0].comp3 - atoms->cell_vector[q2].comp3 - R->vec_ai[g2].comp3;
    F_off = pair_p->off[pair_p->ptr[gl * dim2 + q1 * dim1 + q2]];
    for (i = 0; i < dft_grid->numb[q0]; i++) {
      R_vec_tmp[1].comp1 = R_vec_tmp[0].comp1 + dft_grid->x[n + i].comp1;
      R_vec_tmp[1].comp2 = R_vec_tmp[0].comp2 + dft_grid->x[n + i].comp2;
      R_vec_tmp[1].comp3 = R_vec_tmp[0].comp3 + dft_grid->x[n + i].comp3;
      rsqrd = double_vec_dot(&R_vec_tmp[1],&R_vec_tmp[1]);
      shelposj = atoms->shelposn_sh[q2];
      gausposj = atoms->gausposn_sh[q2];
      count1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[q2]; index_j++) {
        shelj = shells->type_sh[index_j];
        alpha_rsqrd = k_zero;
        for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
          alpha_rsqrd += gaussians->c_sh[gausposj + j4] * exp(-gaussians->expo_sh[gausposj + j4] * rsqrd) ;
         }
          for (k = 0; k < atoms->bfnnumb_sh[q1]; k++) {
            //if (alpha_rsqrd > k_tol) {
            ////count = t * dft_grid->numb[q0] * atoms->bfnnumb_sh[q1] + i * atoms->bfnnumb_sh[q1] + k;
            F_offset = F_off +  k * atoms->bfnnumb_sh[q2] + count1;
            //fprintf(file.out,"count %3d  %3d %3d %3d i  %3d count %3d count + k %3d j %3d count1 %3d  F_O %3d\n", \
            t,q0,q1,q2,i,count,count+k,index_j,count1,F_offset);
            switch (shelj) {
              case 1:
                grid_basis[count + k] += alpha_rsqrd * F[F_offset];
                //fprintf(file.out,"i %3d count %3d r^2 %10.4lf ar^2 %10.4lf F_offset %3d %10.4lf\n",i,count,rsqrd,alpha_rsqrd,F_offset,F[F_offset]);
                break;
              case 3:
                grid_basis[count + k] += R_vec_tmp[1].comp1 * alpha_rsqrd * F[F_offset];
                grid_basis[count + k] += R_vec_tmp[1].comp2 * alpha_rsqrd * F[F_offset + 1];
                grid_basis[count + k] += R_vec_tmp[1].comp3 * alpha_rsqrd * F[F_offset + 2];
                //grid_basis[count] += R_vec_tmp[1].comp1 * alpha_rsqrd * F[F_offset];
                //grid_basis[count] += R_vec_tmp[1].comp2 * alpha_rsqrd * F[F_offset + 1];
                //grid_basis[count] += R_vec_tmp[1].comp3 * alpha_rsqrd * F[F_offset + 2];
                break;
              case 5:
                xsqrd = R_vec_tmp[1].comp1 * R_vec_tmp[1].comp1;
                ysqrd = R_vec_tmp[1].comp2 * R_vec_tmp[1].comp2;
                zsqrd = R_vec_tmp[1].comp3 * R_vec_tmp[1].comp3;
                grid_basis[count + k] += (zsqrd - (xsqrd + ysqrd) / two) /sqrt(three) * alpha_rsqrd * F[F_offset];
                grid_basis[count + k] +=  R_vec_tmp[1].comp1 * R_vec_tmp[1].comp3     * alpha_rsqrd * F[F_offset + 1];
                grid_basis[count + k] +=  R_vec_tmp[1].comp2 * R_vec_tmp[1].comp3     * alpha_rsqrd * F[F_offset + 2];
                grid_basis[count + k] += (xsqrd - ysqrd) / two                        * alpha_rsqrd * F[F_offset + 3];
                grid_basis[count + k] +=  R_vec_tmp[1].comp1 * R_vec_tmp[1].comp2     * alpha_rsqrd * F[F_offset + 4];
                break;
              case 7:
                fprintf(file.out,"f coefficients not entered yet\n");
                exit(1);
                break;
               } // close switch
            //}
             }
            count1 += shelj;
            gausposj += shells->ng_sh[index_j];
           }
          count += atoms->bfnnumb_sh[q1];
          } // close loop on i
         } // close loop on p

}

void generate_dft_rho_grid(double *rho_grid, double *grid_basis, DFT_GRID *dft_grid, TRIPLE_TRAN *Triple, PAIR_TRAN *pair_p, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)
{

int g1, i, i4, j, n, k, p, q, q0, q1, r, t, t1, t2;
int index_i, shelposi, gausposi, sheli;
int count, count1, count2;
double xsqrd, ysqrd, zsqrd, rsqrd, alpha_rsqrd;
double time1, time2;
VECTOR_DOUBLE R_vec_tmp[2];

  time1 = MPI_Wtime();

  count = 0;
  for (p = 0; p < Triple->tot; p++) {
    n  = dft_grid->posn[Triple->cell1[p]];
    q0 = Triple->cell1[p];
    q1 = Triple->cell2[p];
    g1 = Triple->latt2[p];
    //fprintf(file.out,"p %d n %d q0 %d q1 %d g1 %d\n",p,n,q0,q1,g1);
    shelposi = atoms->shelposn_sh[q1];
    for (i = 0; i < dft_grid->numb[q0]; i++) {
      R_vec_tmp[0].comp1 = atoms->cell_vector[q0].comp1 + dft_grid->x[n + i].comp1 - atoms->cell_vector[q1].comp1 - R->vec_ai[g1].comp1;
      R_vec_tmp[0].comp2 = atoms->cell_vector[q0].comp2 + dft_grid->x[n + i].comp2 - atoms->cell_vector[q1].comp2 - R->vec_ai[g1].comp2;
      R_vec_tmp[0].comp3 = atoms->cell_vector[q0].comp3 + dft_grid->x[n + i].comp3 - atoms->cell_vector[q1].comp3 - R->vec_ai[g1].comp3;
      rsqrd = double_vec_dot(&R_vec_tmp[0],&R_vec_tmp[0]);
      //fprintf(file.out,"i %d rsqrd %10.4lf \n",i,rsqrd);
      shelposi = atoms->shelposn_sh[q1];
      gausposi = atoms->gausposn_sh[q1];
      count1 = 0;
      for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[q1]; index_i++) {
        sheli = shells->type_sh[index_i];
        alpha_rsqrd = k_zero;
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          alpha_rsqrd += gaussians->c_sh[gausposi + i4] * exp(-gaussians->expo_sh[gausposi + i4] * rsqrd);
         }
          //if (alpha_rsqrd > k_tol) {
            count2 = count + count1;
            switch (sheli) {
               case 1:
               //rho_grid[n + i] += grid_basis[t * dft_grid->numb[*atm0] * atoms->bfnnumb_sh[q1] +  i * atoms->bfnnumb_sh[q1] + count] * \
               alpha_rsqrd;
               rho_grid[n + i] += grid_basis[count2] * alpha_rsqrd;
               break;
               case 3:
               rho_grid[n + i] += R_vec_tmp[0].comp1 * grid_basis[count2    ] * alpha_rsqrd;
               rho_grid[n + i] += R_vec_tmp[0].comp2 * grid_basis[count2 + 1] * alpha_rsqrd;
               rho_grid[n + i] += R_vec_tmp[0].comp3 * grid_basis[count2 + 2] * alpha_rsqrd;
               break;
               case 5: // THESE NEED TO BE CHECKED !!
               xsqrd = R_vec_tmp[0].comp1 * R_vec_tmp[0].comp1;
               ysqrd = R_vec_tmp[0].comp2 * R_vec_tmp[0].comp2;
               zsqrd = R_vec_tmp[0].comp3 * R_vec_tmp[0].comp3;
               rho_grid[n + i] += (zsqrd - (xsqrd + ysqrd) / two) / rtthree * grid_basis[count2    ] * alpha_rsqrd;
               rho_grid[n + i] +=   R_vec_tmp[0].comp1 * R_vec_tmp[0].comp3 * grid_basis[count2 + 1] * alpha_rsqrd;
               rho_grid[n + i] +=   R_vec_tmp[0].comp2 * R_vec_tmp[0].comp3 * grid_basis[count2 + 2] * alpha_rsqrd;
               rho_grid[n + i] += (xsqrd - ysqrd) / rttwo                   * grid_basis[count2 + 3] * alpha_rsqrd;
               rho_grid[n + i] +=   R_vec_tmp[0].comp1 * R_vec_tmp[0].comp2 * grid_basis[count2 + 4] * alpha_rsqrd;
               break;
               case 7:
               fprintf(file.out,"f coefficients not entered yet\n");
               exit(1);
               break;
              } // close switch
           //}
            count1 += sheli;
            gausposi += shells->ng_sh[index_i];
          } // close loop on index_i
          count += atoms->bfnnumb_sh[q1];
         } // close loop on i over grid points
        } // close loop on p

}

void vxc_grid_dft(double *xc_grid, DFT_GRID *dft_grid, double *rho_grid, JOB_PARAM *job, FILES file)

{

int i, m;
double w;
xc_func_type func[job->xc_num];

  for (i = 0; i < job->xc_num; i++) {
    w = k_zero;
    xc_func_init(&func[i], job->xc_typ[i], job->spin_dim);
    if (func[i].info->family == 32)
    if (job->taskid == 0) fprintf(file.out,"Fix call to xc_func at line 1190 in SCF.cpp\n");
    //w = xc_hyb_gga_exx_coef(func[i].gga);
    //fprintf(file.out,"weight %lf\n",w);
   }

  for (i = 0; i < job->xc_num * dft_grid->npoints; i++) {
    xc_grid[i] = k_zero;
   }

/*
  for (i = 0; i < job->xc_num * dft_grid->npoints; i++) {
    double R_rho = pow(fabs(rho_grid[i]),0.333333333333333);
    xc_grid[i] = -0.738558766 * R_rho;
   }
*/
//rho_grid[0] = two;

  for (m = 0; m <job->xc_num; m++) {
   if (job->xc_typ[m] > 0) {
   switch(func[m].info->family) {
     case XC_FAMILY_LDA:
     //fprintf(file.out,"XC_FAMILY_LDA %3d %3d %3d\n",m,func[m].info->family,job->xc_typ[m]);
     //xc_lda_exc_vxc(&func[m], npoints, rho_grid, &exc_grid[m * dft_grid->npoints], &xc_grid[m * dft_grid->npoints]);
     xc_lda_exc(&func[m], dft_grid->npoints, rho_grid, &xc_grid[m * dft_grid->npoints]);
     //printf("grid %20.15e\n",xc_grid[0]);
     break;
     case XC_FAMILY_GGA:
     case XC_FAMILY_HYB_GGA:
     //xc_gga_exc_vxc(&func[m], dft_grid->npoints, rho_grid, &sigma_grid[0], &exc_grid[m * dft_grid->npoints], vrho, vsigma);
     //xc_gga_exc(&func, dft_grid->npoints, rho, sigma, zk);
     break;
    }
   }
  }

}

void exchange_matrix_dft(double *Kohn_2e, double *xc_grid, DFT_GRID *dft_grid, double *rho_grid, double *P, TRIPLE_TRAN *Triple, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)
{

int i, i1, ip, jp, gi, gj, i4, j4, j, k, p, q, q0, q1, q2, g2, r, s, t1;
int m, n, nd1, nd2, nd12, nd3, nd4, nd34, mm, mm0, mm1, mm2, mm3;
int index_i, shelposi, gausposi;
int index_j, shelposj, gausposj;
int dim1;
int bfposi, bfposi1, bfposj1, bfposj;
int sheli, sheli1, shelj, shelj1;
int n1, n2, n3, n4, n5, n6;
int t, u, v;
int tmax, umax, vmax;
int imax, jmax;
int off1, off2, off3;
int Kohn_count;
double R_AB_1esqrd, rsqrd;
double time1, time2;
double E1_max;
double fac, fac1;
double f[13][13][13][13];
//CHANGES2015double f[9][9][9][9];
double f1[1][55];
double *fgtuv, *p_fgtuv;
VECTOR_DOUBLE R_vec_tmp[3];
VECTOR_DOUBLE R_AB_1e;

  //for (p = 0; p < Triple->tot; p++) {
    p = 0;
    q0 = Triple->cell1[p];
    ip = Triple->cell2[p];
    gi = 0;
    jp = Triple->cell3[p];
    gj = Triple->latt3[p];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double Kohn_2e_cart[nd12];
    double Kohn_2e_sh[nd34];
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli = shells->type1_sh[index_i];
      sheli1 = shells->type_sh[index_i];
      imax  = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj = shells->type1_sh[index_j];
        shelj1 = shells->type_sh[index_j];
        jmax  = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);
        mm = imax + jmax;
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        off2 = (imax + jmax + 1) * off3;
        off1 = (jmax + 1) * off2;
        ////for (t1 = 0; t1 < pair_t.numt; t1++) {
        for (t1 = 0; t1 < 1; t1++) {
          for (i = 0; i < nd12; i++) Kohn_2e_cart[i] = k_zero;
          for (i = 0; i < nd34; i++) Kohn_2e_sh[i] = k_zero;
          fgtuv = (double *) calloc(mm0, sizeof(double));
          if (fgtuv == NULL) {
          fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
          exit(1);
         }
          m = q0;
          //m = pair_t.cell1[pair_t.posn[t1]];
          n = dft_grid->posn[m];
          r = 0;
          //r = pair_t.latt1[pair_t.posn[t1]];
          R_vec_tmp[0].comp1 = atoms->cell_vector[m].comp1 + R->vec_ai[r].comp1;
          R_vec_tmp[0].comp2 = atoms->cell_vector[m].comp2 + R->vec_ai[r].comp2;
          R_vec_tmp[0].comp3 = atoms->cell_vector[m].comp3 + R->vec_ai[r].comp3;
              for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
                R_vec_tmp[1].comp1 = R_vec_tmp[0].comp1 - R_AB[i4].comp1;
                R_vec_tmp[1].comp2 = R_vec_tmp[0].comp2 - R_AB[i4].comp2;
                R_vec_tmp[1].comp3 = R_vec_tmp[0].comp3 - R_AB[i4].comp3;
                for (i1 = 0; i1 < dft_grid->numb[m];i1++) {
                  fac = sab[i4] / pab_inv[i4] / sqrt(pab_inv[i4]) / pi32 * dft_grid->weight[n + i1];
                  fac1 = k_zero;
                  // this loop causing problems in Ne atom; only needed for > 1 DFT functional
                  //for (m = 0; m <job->xc_num; m++) {
                  //fac1 += xc_grid[m * dft_grid->npoints + n + i1];
                 //}
                  fac1 = xc_grid[n + i1];
                  //fac1 = k_one; // uncomment for overlap
                  fac *= fac1;
                  R_vec_tmp[2].comp1 = R_vec_tmp[1].comp1 + dft_grid->x[n + i1].comp1;
                  R_vec_tmp[2].comp2 = R_vec_tmp[1].comp2 + dft_grid->x[n + i1].comp2;
                  R_vec_tmp[2].comp3 = R_vec_tmp[1].comp3 + dft_grid->x[n + i1].comp3;
                  rsqrd = double_vec_dot(&R_vec_tmp[2],&R_vec_tmp[2]);
                  f1[0][0] = exp(-rsqrd / pab_inv[i4]);
                  f1[0][1] = -two / pab_inv[i4] * f1[0][0];
                  f1[0][2] = -two / pab_inv[i4] * f1[0][1];
                  f1[0][3] = -two / pab_inv[i4] * f1[0][2];
                  f1[0][4] = -two / pab_inv[i4] * f1[0][3];
                  f1[0][5] = -two / pab_inv[i4] * f1[0][4];
                  f1[0][6] = -two / pab_inv[i4] * f1[0][5];
                  f1[0][7] = -two / pab_inv[i4] * f1[0][6];
                  f1[0][8] = -two / pab_inv[i4] * f1[0][7];
                  f1[0][9] = -two / pab_inv[i4] * f1[0][8];
                  f1[0][10] = -two / pab_inv[i4] * f1[0][9];
                  f1[0][11] = -two / pab_inv[i4] * f1[0][10];
                  f1[0][12] = -two / pab_inv[i4] * f1[0][11];
                  non_recursive_ftuvn(mm, 0, f, f1, &R_vec_tmp[2]);
                  for (t = 0; t <= mm; t++) {
                    for (u = 0; u <= mm; u++) {
                      for (v = 0; v <= mm; v++) {
                        if (t + u + v > mm) break;
                        p_fgtuv = fgtuv + t * mm1 + u * mm2 + v * mm3  + i4;
                       *p_fgtuv += f[t][u][v][0] * fac;
                       }
                      }
                     }
                    } // close loop over i1
                   } // close loop over i4
                   for (i = 0; i < sheli; i++) {
                     for (j = 0; j < shelj; j++) {
                       Kohn_count = (bfposi + i) * nd2 + bfposj + j;
                       tmax = shells->tuv[sheli][i][0] + shells->tuv[shelj][j][0];
                       umax = shells->tuv[sheli][i][1] + shells->tuv[shelj][j][1];
                       vmax = shells->tuv[sheli][i][2] + shells->tuv[shelj][j][2];
                       n1 = shells->tuv[sheli][i][0];
                       n2 = shells->tuv[sheli][i][1];
                       n3 = shells->tuv[sheli][i][2];
                       n4 = shells->tuv[shelj][j][0];
                       n5 = shells->tuv[shelj][j][1];
                       n6 = shells->tuv[shelj][j][2];
                       for (t = 0; t <= tmax; t++) {
                         for (u = 0; u <= umax; u++) {
                           for (v = 0; v <= vmax; v++) {
                             for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
                               p_fgtuv   = fgtuv + t * mm1 + u * mm2 + v * mm3 + i4;
                               //fprintf(file.out,"%4d %4d %4d %4d",nd12,Kohn_count,t * mm1 + u * mm2 + v * mm3 + i4,n1 * off1 + n4 * off2 + t * off3 + i4);
                               Kohn_2e_cart[Kohn_count] += *p_fgtuv * E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                                                                      E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                                                                      E1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                             } // End i4 loop
                            }
                           }
                          }
                         } // close loop over j
                        } // close loop over i
                       //fprintf(file.out,"Kohn_2e_cart\n"); int count = 0;
                       //for(j=0;j<nd1;j++) { for(k=0;k<nd2;k++) { fprintf(file.out," %11.4e",Kohn_2e_cart[count]); count++; }
                       //fprintf(file.out,"\n"); } fprintf(file.out,"\n");
                       ////two_center_cartesian_to_sh_shell(Kohn_2e_cart,Kohn_2e_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,\
                       shells,job,file);
                       //fprintf(file.out,"Kohn_2e_rot\n"); count = 0;
                       //for(j=0;j<nd3;j++) { for(k=0;k<nd4;k++) { fprintf(file.out," %11.4e",Kohn_2e_sh[count]); count++; }
                       //fprintf(file.out,"\n"); } fprintf(file.out,"\n");
                       //////rotate_triple_shell(Kohn_2e_sh,&Kohn_2e[dim1],index_i,index_j,bfposi,bfposj,nd4,&pair_t,t1,shells,symmetry,job,file);
                       //fprintf(file.out,"Kohn_2e\n"); int count = 0;
                       //for(j=0;j<nd3;j++) { for(k=0;k<nd4;k++) { fprintf(file.out," %11.4e",Kohn_2e[count]); count++; }
                       //fprintf(file.out,"\n"); } fprintf(file.out,"\n");
                       free(fgtuv);
                      } // close loop over t1
                     bfposj   += shelj;
                     bfposj1  += shelj1;
                     gausposj += shells->ng_sh[index_j];
                    } // close loop over index_j
                   bfposi   += sheli;
                   bfposi1  += sheli1;
                   gausposi += shells->ng_sh[index_i];
                  } // close loop over index_i
                 //dim1 += nd34;

                 //} // close loop over p

}

void dft_grid_parameters(int *L_range, double *r_range, JOB_PARAM *job, FILES file)

{

  switch(job->xc_grd) {

     case 0:  // STANDARD GRID
      r_range[0] = 0.4;
      r_range[1] = 0.6;
      r_range[2] = 0.8;
      r_range[3] = 0.9;
      r_range[4] = 1.1;
      r_range[5] = 2.3;
      r_range[6] = 2.4;
      r_range[7] = 2.6;
      r_range[8] = 2.8;
      r_range[9] = 9999.0;
      L_range[0] = 1;
      L_range[1] = 2;
      L_range[2] = 5;
      L_range[3] = 8;
      L_range[4] = 11;
      L_range[5] = 13;
      L_range[6] = 11;
      L_range[7] = 8;
      L_range[8] = 5;
      L_range[9] = 1;
      break;

     case 1:  // LARGE GRID
      r_range[0] = 0.17;
      r_range[1] = 0.5;
      r_range[2] = 0.9;
      r_range[3] = 3.0;
      r_range[4] = 9999.0;
      L_range[0] = 2;
      L_range[1] = 6;
      L_range[2] = 8;
      L_range[3] = 13;
      L_range[4] = 8;
      L_range[5] = 8;
      break;

     case 2:  // XLARGE GRID
      r_range[0] = 0.17;
      r_range[1] = 0.5;
      r_range[2] = 0.9;
      r_range[3] = 3.5;
      r_range[4] = 9999.0;
      L_range[0] = 4;
      L_range[1] = 8;
      L_range[2] = 12;
      L_range[3] = 13;  // 16
      L_range[4] = 12;
      L_range[5] = 12;
      break;

     case 3:  // XXLARGE GRID
      r_range[0] = 0.17;
      r_range[1] = 0.5;
      r_range[2] = 0.9;
      r_range[3] = 3.5;
      r_range[4] = 9999.0;
      L_range[0] = 6;
      L_range[1] = 10;
      L_range[2] = 14;
      L_range[3] = 13; // 18
      L_range[4] = 14;
      L_range[5] = 14;
      break;

 }

}

