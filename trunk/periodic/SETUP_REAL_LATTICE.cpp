/*
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
*/
#include "conversion_factors.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "LIMITS.h"
#include "MATRIX_UTIL.h"
#include "SETUP_REAL_LATTICE.h"

using namespace std;

  // ******************************************************************************************
  // *  Generate real space lattice vectors and order in shells                               *
  // ******************************************************************************************

void generate_real_lattice(CRYSTAL *crystal, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Generate real space lattice vectors and order in shells                               *
  // ******************************************************************************************

  int count;
  int MXRX, MXRY, MXRZ;
  int i, j, k;
  double tempRmag;
  double smallest = 10000.0;
  double root_gamma_0;
  VECTOR_DOUBLE tempRvec_ai;

  root_gamma_0 = sqrt(k_one / G->gamma_0_inv);

  R->last_ewald_vector = 0;

  count = 0;

  switch (crystal->type[0]) {

    case 'C':

      MXRX = job->mxr; 
      MXRY = job->mxr; 
      MXRZ = job->mxr;
      for (i = -MXRX; i < MXRX + 1; i++) {
        for (j = -MXRY; j < MXRY + 1; j++) {
          for (k = -MXRZ; k < MXRZ + 1; k++) {
            R->vec_ai[count].comp1 = i * crystal->primitive_cell[0].comp1 + j * crystal->primitive_cell[1].comp1 + k
                * crystal->primitive_cell[2].comp1;
            R->vec_ai[count].comp2 = i * crystal->primitive_cell[0].comp2 + j * crystal->primitive_cell[1].comp2 + k
                * crystal->primitive_cell[2].comp2;
            R->vec_ai[count].comp3 = i * crystal->primitive_cell[0].comp3 + j * crystal->primitive_cell[1].comp3 + k
                * crystal->primitive_cell[2].comp3;
            R->mag[count] = sqrt(double_vec_dot(&R->vec_ai[count], &R->vec_ai[count]));
            if (erfc(root_gamma_0 * R->mag[count]) / R->mag[count] > real_space_cutoff) 
            R->last_ewald_vector++;
            if (i > MXRX / 2 && j > MXRY / 2 && k > MXRZ / 2 && R->mag[count] < smallest) smallest = R->mag[count];
            count++;
          }
        }
      }
      //fprintf(file.out,"MXRX %3d %3d %3d smallest radius %10.4lf\n",MXRX,MXRY,MXRZ,smallest);
      //for (i = 0; i < count; i++)
      //fprintf(file.out,"Rmag %3d %10.4lf %10.4lf %10.4lf %10.4lf \n",i,\
      R->vec_ai[i].comp1,R->vec_ai[i].comp2,R->vec_ai[i].comp3,R->mag[i]);

      break;

    case 'S':

      MXRX = job->mxr; 
      MXRY = job->mxr;
      for (i = -MXRX; i < MXRX + 1; i++) {
        for (j = -MXRY; j < MXRY + 1; j++) {
            R->vec_ai[count].comp1 = i * crystal->primitive_cell[0].comp1 + j * crystal->primitive_cell[1].comp1; 
            R->vec_ai[count].comp2 = i * crystal->primitive_cell[0].comp2 + j * crystal->primitive_cell[1].comp2;
            R->vec_ai[count].comp3 = k_zero;
            R->mag[count] = sqrt(double_vec_dot(&R->vec_ai[count], &R->vec_ai[count]));
            if (erfc(root_gamma_0 * R->mag[count]) / R->mag[count] > real_space_cutoff) 
            R->last_ewald_vector++;
            if (i > MXRX / 2 && j > MXRY / 2 && R->mag[count] < smallest) smallest = R->mag[count];
            count++;
          }
        }

      break;

    case 'P':

      MXRX = job->mxr; 
      MXRY = job->mxr; 
      MXRZ = job->mxr;
      for (i = -MXRZ; i < MXRZ + 1; i++) {
            R->vec_ai[count].comp1 = k_zero;
            R->vec_ai[count].comp2 = k_zero;
            R->vec_ai[count].comp3 = i * crystal->primitive_cell[0].comp3;
            R->mag[count] = sqrt(double_vec_dot(&R->vec_ai[count], &R->vec_ai[count]));
            R->last_ewald_vector++;
            // set condition on last_ewald_vector
            count++;
          }

      break;

    case 'M':

            R->vec_ai[0].comp1 = k_zero;
            R->vec_ai[0].comp2 = k_zero;
            R->vec_ai[0].comp3 = k_zero;
            R->last_ewald_vector = 1;
            count = 1;

      break;

  } // close switch(crystal.type

  for (i = 0; i < count; i++)
    R->vec_ai_ord[i] = i;

  for (i = 1; i < count; i++) {
    for (j = count - 1; j >= i; j--) {

      if (R->mag[j - 1] > R->mag[j]) {

        tempRmag = R->mag[j - 1];
        R->mag[j - 1] = R->mag[j];
        R->mag[j] = tempRmag;

        k = R->vec_ai_ord[j - 1];
        R->vec_ai_ord[j - 1] = R->vec_ai_ord[j];
        R->vec_ai_ord[j] = k;

        tempRvec_ai = R->vec_ai[j - 1];
        R->vec_ai[j - 1] = R->vec_ai[j];
        R->vec_ai[j] = tempRvec_ai;

      }
    }
  }

  if (job->taskid == 0 &&job->verbosity > 1) {
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"|                                 REAL LATTICE VECTORS AND SHELLS                                         |\n");
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"\n");
   }
    //fprintf(file.out,"%5d REAL SPACE VECTORS GENERATED\n",count);

}

void generate_real_lattice_shells(CRYSTAL *crystal, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Generate real lattice vector shell arrays                                             *
  // ******************************************************************************************

  int i, j;
  int *Rvec_mag;
  double *Rvecmag;

    Rvec_mag = (int *) malloc((R->max_vector + 1) * sizeof(int));
    //Changes2021 Rvec_mag = (int *) malloc(R->max_vector * sizeof(int));
    if (Rvec_mag == NULL) {
    fprintf(stderr, "ERROR: not enough memory for int Rvec_mag \n");
    exit(1);
   }

    Rvecmag = (double *) malloc(R->max_vector * sizeof(double));
    if (Rvecmag == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double Rvecmag \n");
    exit(1);
   }

  switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':

  R->number_of_shells = 0;
  for (i = 0; i < R->max_vector; i++) 
  Rvec_mag[i] = -1;
  j = 1;
  Rvec_mag[0] = 0;
  for (i = 1; i < R->max_vector; i++) {
    if (fabs(R->mag[i] - R->mag[i - 1]) > 0.000001) {
    Rvec_mag[j] = i;
    Rvecmag[j] = R->mag[i];
    if (i <= R->last_vector) 
    //CHANGESJULY2017if (i < R->last_vector) 
    R->number_of_shells = j;
    if (i < R->last_ewald_vector) 
    R->number_of_ewald_shells = j;
    j++;
   }
  }

  R->num[0] = 1;  // zero vector
  for (i = 1; i < R->max_vector; i++) 
  R->num[i] = Rvec_mag[i + 1] - Rvec_mag[i];

  if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"                  SHELLS %5d WITHIN R->last_vector       %5d MAGNITUDE %12.6lf\n", \
    R->number_of_shells,R->last_vector,R->mag[R->last_vector]);
    fprintf(file.out,"                  SHELLS %5d WITHIN R->ewald_last_vector %5d MAGNITUDE %12.6lf\n\n", \
    R->number_of_ewald_shells,R->last_ewald_vector,R->mag[R->last_ewald_vector]);
    fprintf(file.out, "                           REAL SPACE SHELLS NUMBER CLOSING MAGNITUDE\n");
    fprintf(file.out, "\n");
    for (i = 0; i < R->number_of_ewald_shells; i++)
    fprintf(file.out, "%36d %5d %5d %12.6lf\n", i, R->num[i], Rvec_mag[i + 1], Rvecmag[i]);
    fprintf(file.out, "\n");
    fprintf(file.out, "                            REAL SPACE VECTORS IN SHELLS %5d\n", R->max_vector);
    fprintf(file.out, "\n");
    for (i = 0; i < R->max_vector; i++)
    fprintf(file.out, "%18d %12.6lf%12.6lf%12.6lf     %12.6lf\n", i, R->vec_ai[i].comp1, R->vec_ai[i].comp2, R->vec_ai[i].comp3, R->mag[i]);
    fprintf(file.out, "\n");
  }
     break;

  } // close switch(crystal.type

  free(Rvec_mag);
  free(Rvecmag);

}

void generate_real_lattice_tables(CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Generate sum, difference and rotation tables for lattice vectors                      *
  // ******************************************************************************************

  int MXRX, MXRY, MXRZ;
  int shr1, shr2;
  int dif_index, sum_index;
  int i, j, k;
  int l1, m1, n1, l2, m2, n2, vec1, vec2;
  int *Rvec_ai_inv;
  double *p_irr;
  VECTOR_DOUBLE Rvec_tmp;

    Rvec_ai_inv = (int *) malloc(R->max_vector * sizeof(int));
    if (Rvec_ai_inv == NULL) {
    fprintf(stderr, "ERROR: not enough memory for int Rvec_ai_inv \n");
    exit(1);
   }

  for (i = 0; i < R->max_vector; i++) {
    Rvec_ai_inv[R->vec_ai_ord[i]] = i;
    R->vec_ai_inv[R->vec_ai_ord[i]] = i;
   }

  for (i = 0; i < R_tables->margin_vector; i++) {
    for (j = 0; j < R_tables->margin_vector; j++) {
      R_tables->sumvec[i * R_tables->margin_vector + j] = -1;
    }
  }

  for (i = 0; i < R_tables->margin_vector; i++) {
    for (j = 0; j < R_tables->margin_vector; j++) {
      R_tables->diffvec[i * R_tables->margin_vector + j] = -1;
    }
  }

  switch (crystal->type[0]) {

    case 'C':

      MXRX = job->mxr; 
      MXRY = job->mxr; 
      MXRZ = job->mxr;
      shr1 =  2 * MXRZ + 1;
      shr2 = (2 * MXRY + 1) * (2 * MXRZ + 1);
      for (i = 0; i < R_tables->margin_vector; i++) {
        vec1 = R->vec_ai_ord[i];
        l1 =  vec1 / shr2 - MXRX;
        m1 = (vec1 - (l1 + MXRX) * shr2) / shr1 - MXRY;
        n1 =  vec1 - (l1 + MXRX) * shr2 - (m1 + MXRY) * shr1 - MXRZ;
        for (j = 0; j < R_tables->margin_vector; j++) {
          vec2 = R->vec_ai_ord[j];
          l2 =  vec2 / shr2 - MXRX;
          m2 = (vec2 - (l2 + MXRX) * shr2) / shr1 - MXRY;
          n2 =  vec2 - (l2 + MXRX) * shr2 - (m2 + MXRY) * shr1 - MXRZ;
          dif_index = (MXRX + l1 - l2) * shr2 + (MXRY + m1 - m2) * shr1 + MXRZ + n1 - n2;
          sum_index = (MXRX + l1 + l2) * shr2 + (MXRY + m1 + m2) * shr1 + MXRZ + n1 + n2;
          R_tables->diffvec[i * R_tables->margin_vector + j] = Rvec_ai_inv[dif_index];
          R_tables->sumvec[i * R_tables->margin_vector + j]  = Rvec_ai_inv[sum_index];
          //fprintf(file.out,"i j %5d %5d vec12 %5d %5d   %5d %5d %5d   %5d %5d %5d  sumvec difvec %5d %5d   %5d %5d\n", \
          i,j,vec1,vec2,l1,m1,n1,l2,m2,n2,sum_index,dif_index,Rvec_ai_inv[sum_index],Rvec_ai_inv[dif_index]);
          if (sum_index > R->max_vector || dif_index > R->max_vector) {
          if (job->taskid == 0)
          fprintf(file.out,"sum_index %d or dif_index %d exceeded R->max_vector %d. Increase job->mxr with SET_MXR\n",\
          sum_index,dif_index,R->max_vector);
          MPI_Finalize();
          exit(1);
         }
        }
       }

     break;

    case 'S':

      MXRX = job->mxr; 
      MXRY = job->mxr; 
      shr1 =  2 * MXRY + 1;
      for (i = 0; i < R_tables->margin_vector; i++) {
        vec1 = R->vec_ai_ord[i];
        l1 = vec1 / shr1 - MXRX;
        m1 = vec1 - (l1 + MXRX) * shr1 - MXRY;
        for (j = 0; j < R_tables->margin_vector; j++) {
          vec2 = R->vec_ai_ord[j];
          l2 = vec2 / shr1 - MXRX;
          m2 = vec2 - (l2 + MXRX) * shr1 - MXRY;
          dif_index = (MXRX + l1 - l2) * shr1 + (MXRY + m1 - m2);
          sum_index = (MXRX + l1 + l2) * shr1 + (MXRY + m1 + m2);
          R_tables->diffvec[i * R_tables->margin_vector + j] = Rvec_ai_inv[dif_index];
          R_tables->sumvec[i * R_tables->margin_vector + j]  = Rvec_ai_inv[sum_index];
          //fprintf(file.out,"i j %5d %5d vec12 %5d %5d   %5d %5d   %5d %5d  sumvec difvec %5d %5d   %5d %5d\n", \
          i,j,vec1,vec2,l1,m1,l2,m2,sum_index,dif_index,Rvec_ai_inv[sum_index],Rvec_ai_inv[dif_index]);
          if (sum_index > R->max_vector || dif_index > R->max_vector) {
          if (job->taskid == 0)
          fprintf(file.out,"sum_index %d or dif_index %d exceeded R->max_vector %d. Increase job->mxr with SET_MXR\n",\
          sum_index,dif_index,R->max_vector);
          MPI_Finalize();
          exit(1);
         }
        }
       }

     break;

    case 'P':

      MXRZ = job->mxr; 
      for (i = 0; i < R_tables->margin_vector; i++) {
        vec1 = R->vec_ai_ord[i];
        l1 = vec1 - MXRZ;
        for (j = 0; j < R_tables->margin_vector; j++) {
          vec2 = R->vec_ai_ord[j];
          l2 = vec2 - MXRZ;
          dif_index = MXRZ + l1 - l2;
          sum_index = MXRZ + l1 + l2;
          R_tables->diffvec[i * R_tables->margin_vector + j] = Rvec_ai_inv[dif_index];
          R_tables->sumvec[i * R_tables->margin_vector + j]  = Rvec_ai_inv[sum_index];
          //fprintf(file.out,"i j %5d %5d vec12 %5d %5d   %5d   %5d  sumvec difvec %5d %5d   %5d %5d\n", \
          i,j,vec1,vec2,l1,l2,sum_index,dif_index,Rvec_ai_inv[0*sum_index],Rvec_ai_inv[0*dif_index]);
          if (sum_index > R->max_vector || dif_index > R->max_vector) {
          if (job->taskid == 0)
          fprintf(file.out,"sum_index %d or dif_index %d exceeded R->max_vector %d. Increase job->mxr with SET_MXR\n",\
          sum_index,dif_index,R->max_vector);
          MPI_Finalize();
          exit(1);
         }
        }
       }

     break;

    case 'M':

        R_tables->diffvec[0] = 0;
        R_tables->sumvec[0] = 0;
        for (k = 0; k < symmetry->number_of_operators; k++) 
        R_tables->lattvec[k] = 0;
        //fprintf(file.out,"R_tables->max_vector %3d R_tables->last_vector %3d\n",R_tables->max_vector,R_tables->last_vector); fflush(file.out);
        //fprintf(file.out,"R_tables->sumvec %3d R_tables->diffvec %3d\n",R_tables->sumvec[0],R_tables->diffvec[0]); fflush(file.out);
        //for (k = 0; k < symmetry->number_of_operators; k++) 
        //if (job->taskid == 0)
        //fprintf(file.out,"R_tables->lattvec[%2d] %3d\n",k,R_tables->lattvec[k]); fflush(file.out);

     break;

  } // close switch(crystal.type

  for (i = 0; i < R_tables->margin_vector; i++) {
    for (j = 0; j < R_tables->margin_vector; j++) {
      if (R_tables->diffvec[i * R_tables->margin_vector + j] == -1) {
        fprintf(file.out,"%d %d %d diffvec array not initialised: R_tables->max_vector = %5d decrease R->cutoff %lf or increase MXR %d\n", \
        i,j,R_tables->diffvec[i * R_tables->margin_vector + j], R_tables->margin_vector, R->cutoff * bohr_to_AA, job->mxr);
        exit(1);
      }
    }
  }

  for (i = 0; i < R_tables->margin_vector; i++) {
    for (j = 0; j < R_tables->margin_vector; j++) {
      if (R_tables->sumvec[i * R_tables->margin_vector + j] == -1) {
        if (job->taskid == 0)
        fprintf(file.out,"%d %d sumvec array not initialised: R_tables->max_vector = %5d decrease R->cutoff %lf or increase MXR %d\n", \
        i,j,R_tables->margin_vector, R->cutoff * bohr_to_AA, job->mxr);
        exit(1);
      }
    }
  }

  switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':

  for (i = 0; i < R_tables->margin_vector; i++) {
    for (k = 0; k < symmetry->number_of_operators; k++) {
      p_irr = symmetry->irr + k * 9;
      rotate_vector3(p_irr, &R->vec_ai[i], &Rvec_tmp);
      //fprintf(file.out,"i k %6d %3d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",i,k,R->vec_ai[i].comp1,R->vec_ai[i].comp2,\
      R->vec_ai[i].comp3,Rvec_tmp.comp1,Rvec_tmp.comp2,Rvec_tmp.comp3) ;
      for (j = 0; j < R_tables->max_vector; j++) {
        if (check_vec(&R->vec_ai[j], &Rvec_tmp) == 1) {
          R_tables->lattvec[i * symmetry->number_of_operators + k] = j;
          //fprintf(file.out,"check_vec %3d %3d   %3d\n",i,k,j);
          break;
        }
       } // close loop on j
       //fprintf(file.out,"i k lattvec[i].op[k]  %3d %3d   %5d\n",i,k,R_tables->lattvec[i * symmetry->number_of_operators + k]) ;
      } // close loop on k
     } // close loop on i

  for (i = 0; i < R_tables->margin_vector; i++) {
    for (k = 0; k < symmetry->number_of_operators; k++) {
      if (R_tables->lattvec[i * symmetry->number_of_operators + k] == -1) {
        if (job->taskid == 0) {
          fprintf(file.out,"%d %d lattvec array not initialised: R_tables->margin_vector = %5d reduce R->cutoff %lf\n",\
          i,j,R_tables->margin_vector, R->cutoff);
          MPI_Finalize();
          exit(1);
         }
        }
       }
      }

     break;

  } // close switch(crystal.type

  free(Rvec_ai_inv);

}
