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
#include "conversion_factors.h"
*/
#include "myconstants.h"
#include "USER_DATA.h"
#include "LIMITS.h"
#include "MATRIX_UTIL.h"
#include "SETUP_RECIPROCAL_LATTICE.h"

using namespace std;

  // ******************************************************************************************
  // * Count and generate reciprocal lattice vectors                                          *
  // ******************************************************************************************

void count_reciprocal_lattice_vectors(CRYSTAL crystal, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Count reciprocal lattice vectors                                                       *
  // ******************************************************************************************

  int MXGX, MXGY, MXGZ;
  int i, j, k;
  double tempGmag;
  double expfac, tempexpnt1;
  VECTOR_DOUBLE tempgvec_bi;

  G->max_vector = 0;
  G->last_vector = 0;
  expfac = eight * pi / crystal.primitive_cell_volume;

  switch (crystal.type[0]) {

    case 'C':

      G->gamma_0_inv = 1.0 * pow(crystal.primitive_cell_volume, two_thirds) / four ;

      MXGX = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[0], &crystal.reciprocal_cell[0])));
      MXGY = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[1], &crystal.reciprocal_cell[1])));
      MXGZ = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[2], &crystal.reciprocal_cell[2])));
      for (i = -MXGX; i < MXGX + 1; i++) {
        for (j = -MXGY; j < MXGY + 1; j++) {
          for (k = -MXGZ; k < MXGZ + 1; k++) {
            tempgvec_bi.comp1 = i * crystal.reciprocal_cell[0].comp1 + j * crystal.reciprocal_cell[1].comp1 + k
                * crystal.reciprocal_cell[2].comp1;
            tempgvec_bi.comp2 = i * crystal.reciprocal_cell[0].comp2 + j * crystal.reciprocal_cell[1].comp2 + k
                * crystal.reciprocal_cell[2].comp2;
            tempgvec_bi.comp3 = i * crystal.reciprocal_cell[0].comp3 + j * crystal.reciprocal_cell[1].comp3 + k
                * crystal.reciprocal_cell[2].comp3;
            tempGmag = sqrt(double_vec_dot(&tempgvec_bi, &tempgvec_bi));
            tempexpnt1 = expfac * exp(-tempGmag*tempGmag / four * G->gamma_0_inv) / tempGmag / tempGmag ;
            if (tempGmag < 0.000000001) tempexpnt1 = k_one;
            if (tempexpnt1 > reciprocal_space_cutoff) {
            if (i > 0 || (i == 0 && j > 0) || (i == 0 && j == 0 && k >= 0)) {
              (G->last_vector)++;
              }
              (G->max_vector)++;
            }
          }
        }
      }
      if (job->taskid == 0) printf("G last vector %4d ",G->last_vector);

      break;

    case 'S':

      G->gamma_0_inv = crystal.primitive_cell_volume / 5.0 ;

      MXGX = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[0], &crystal.reciprocal_cell[0])));
      MXGY = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[1], &crystal.reciprocal_cell[1])));
      for (i = -MXGX; i < MXGX + 1; i++) {
        for (j = -MXGY; j < MXGY + 1; j++) {
          tempgvec_bi.comp1 = i * crystal.reciprocal_cell[0].comp1 + j * crystal.reciprocal_cell[1].comp1 + k
              * crystal.reciprocal_cell[2].comp1;
          tempgvec_bi.comp2 = i * crystal.reciprocal_cell[0].comp2 + j * crystal.reciprocal_cell[1].comp2 + k
              * crystal.reciprocal_cell[2].comp2;
          tempgvec_bi.comp3 = i * crystal.reciprocal_cell[0].comp3 + j * crystal.reciprocal_cell[1].comp3 + k
              * crystal.reciprocal_cell[2].comp3;
          tempGmag = sqrt(double_vec_dot(&tempgvec_bi, &tempgvec_bi));
          tempexpnt1 = expfac * exp(-tempGmag*tempGmag / four * G->gamma_0_inv) / tempGmag / tempGmag ;
          if (tempGmag < 0.000000001) tempexpnt1 = k_one;
          if (tempexpnt1 > reciprocal_space_cutoff) {
          if (i > 0 || (i == 0 && j >= 0)) {
            (G->last_vector)++;
            }
            (G->max_vector)++;
          }
        }
      }

      break;

    case 'P':

      G->gamma_0_inv = crystal.primitive_cell_volume / four ; // needs to be set

      MXGZ = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[0], &crystal.reciprocal_cell[0])));
      for (i = -MXGZ; i < MXGZ + 1; i++) {
        tempgvec_bi.comp1 = k_zero;
        tempgvec_bi.comp2 = k_zero;
        tempgvec_bi.comp3 = i * crystal.reciprocal_cell[0].comp3;
        tempGmag = sqrt(double_vec_dot(&tempgvec_bi, &tempgvec_bi));
        tempexpnt1 = expfac * exp(-tempGmag * tempGmag / four * G->gamma_0_inv) / tempGmag / tempGmag ;
        if (tempGmag < 0.000000001) tempexpnt1 = k_one;
        if (tempexpnt1 > reciprocal_space_cutoff) {
        if (i > 0) {
            (G->last_vector)++;
            }
            (G->max_vector)++;
           }
         }

      break;

    case 'M':
    
     G->max_vector = 1;
     G->last_vector = 1;

      break;

    } // close switch(crystal.type

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"|                              RECIPROCAL LATTICE VECTORS AND SHELLS                                      |\n");
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"\n");
   }

}

void generate_reciprocal_lattice(CRYSTAL crystal, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Generate reciprocal lattice vectors and populate G structure                          *
  // ******************************************************************************************

  int MXGX, MXGY, MXGZ;
  int count;
  int count_b2;
  int i, j, k;
  int m, n, p, max = 11;
  int index_G, G_shell;
  int gvec_ord[G->last_vector];
  int *Gvec_mag;
  double tempGmag;
  double Gmag_bi[G->max_vector];
  double Gsq0, Gsq1, expnt1;
  double expfac, tempexpnt1;
  VECTOR_DOUBLE tempgvec_bi;

  count  = 0;
  count_b2 = 0;
  expfac = eight * pi / crystal.primitive_cell_volume;

  switch (crystal.type[0]) {

    case 'C':

      MXGX = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[0], &crystal.reciprocal_cell[0])));
      MXGY = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[1], &crystal.reciprocal_cell[1])));
      MXGZ = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[2], &crystal.reciprocal_cell[2])));

      for (i = -MXGX; i < MXGX + 1; i++) {
        for (j = -MXGY; j < MXGY + 1; j++) {
          for (k = -MXGZ; k < MXGZ + 1; k++) {
            tempgvec_bi.comp1 = i * crystal.reciprocal_cell[0].comp1 + j * crystal.reciprocal_cell[1].comp1 + k
                * crystal.reciprocal_cell[2].comp1;
            tempgvec_bi.comp2 = i * crystal.reciprocal_cell[0].comp2 + j * crystal.reciprocal_cell[1].comp2 + k
                * crystal.reciprocal_cell[2].comp2;
            tempgvec_bi.comp3 = i * crystal.reciprocal_cell[0].comp3 + j * crystal.reciprocal_cell[1].comp3 + k
                * crystal.reciprocal_cell[2].comp3;
            tempGmag = sqrt(double_vec_dot(&tempgvec_bi, &tempgvec_bi));
            tempexpnt1 = expfac * exp(-tempGmag * tempGmag / four * G->gamma_0_inv) / tempGmag / tempGmag ;
            if (tempGmag < 0.000000001) tempexpnt1 = k_one;
            if (tempexpnt1 > reciprocal_space_cutoff) {
            G->vec_bi[count].comp1 = i * crystal.reciprocal_cell[0].comp1 + j * crystal.reciprocal_cell[1].comp1 + k
                * crystal.reciprocal_cell[2].comp1;
            G->vec_bi[count].comp2 = i * crystal.reciprocal_cell[0].comp2 + j * crystal.reciprocal_cell[1].comp2 + k
                * crystal.reciprocal_cell[2].comp2;
            G->vec_bi[count].comp3 = i * crystal.reciprocal_cell[0].comp3 + j * crystal.reciprocal_cell[1].comp3 + k
                * crystal.reciprocal_cell[2].comp3;
            if (i > 0 || (i == 0 && j > 0) || (i == 0 && j == 0 && k >= 0)) {
               G->vec_b2[count_b2] = G->vec_bi[count];
	       G->mag[count_b2] = tempGmag;
               //fprintf(file.out,"tempexp %3d %3d %10.4e   %10.4lf %10.4lf %10.4lf  %10.4lf\n",count,count_b2,tempexpnt1,\
               G->vec_b2[count_b2].comp1,G->vec_b2[count_b2].comp2, G->vec_b2[count_b2].comp3,G->mag[count_b2]); fflush(file.out);
               count_b2++;
               }
               Gmag_bi[count] = tempGmag;
               count++;
              }
             }
            }
           }
            fprintf(file.out,"reset real and reciprocal space cutoffs to previous values!\n");
           // fprintf(file.out,"%3d %3d %3d %3d %3d %14.8lf gm \n",i,j,k, count,count_b2,G->mag[0]);
           //fprintf(file.out,"count %3d %3d %3d %3d\n",count,count_b2,G->last_vector,G->max_vector);

      break;

    case 'S':

      MXGX = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[0], &crystal.reciprocal_cell[0])));
      MXGY = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[1], &crystal.reciprocal_cell[1])));
      for (i = -MXGX; i < MXGX + 1; i++) {
        for (j = -MXGY; j < MXGY + 1; j++) {
          tempgvec_bi.comp1 = i * crystal.reciprocal_cell[0].comp1 + j * crystal.reciprocal_cell[1].comp1 + k
              * crystal.reciprocal_cell[2].comp1;
          tempgvec_bi.comp2 = i * crystal.reciprocal_cell[0].comp2 + j * crystal.reciprocal_cell[1].comp2 + k
              * crystal.reciprocal_cell[2].comp2;
          tempgvec_bi.comp3 = k_zero;
          tempGmag = sqrt(double_vec_dot(&tempgvec_bi, &tempgvec_bi));
          tempexpnt1 = expfac / four * exp(-tempGmag*tempGmag / four * G->gamma_0_inv) / tempGmag ;
          if (tempGmag < 0.000000001) tempexpnt1 = k_one;
          if (tempexpnt1 > reciprocal_space_cutoff) {
          G->vec_bi[count].comp1 = i * crystal.reciprocal_cell[0].comp1 + j * crystal.reciprocal_cell[1].comp1 + k
              * crystal.reciprocal_cell[2].comp1;
          G->vec_bi[count].comp2 = i * crystal.reciprocal_cell[0].comp2 + j * crystal.reciprocal_cell[1].comp2 + k
              * crystal.reciprocal_cell[2].comp2;
          G->vec_bi[count].comp3 = k_zero;
          if (i > 0 || (i == 0 && j >= 0)) {
          //fprintf(file.out,"tempexp %lf %lf %lf %lf\n",tempexpnt1,G->vec_bi[count].comp1,G->vec_bi[count].comp2,G->vec_bi[count].comp3);
             G->vec_b2[count_b2] = G->vec_bi[count];
	     G->mag[count_b2] = tempGmag;
             count_b2++;
             }
             Gmag_bi[count] = tempGmag;
             count++;
            }
          }
        }
         //fprintf(file.out,"count %d %d %d\n",G->num[0],G->last_vector,G->last_vector);

         break;

    case 'P':

      MXGZ = (int) ceil(G->cutoff / sqrt(double_vec_dot(&crystal.reciprocal_cell[0], &crystal.reciprocal_cell[0])));

      for (i = -MXGZ; i < MXGZ + 1; i++) {
        tempgvec_bi.comp1 = k_zero;
        tempgvec_bi.comp2 = k_zero;
        tempgvec_bi.comp3 = i * crystal.reciprocal_cell[0].comp3;
        tempGmag = sqrt(double_vec_dot(&tempgvec_bi, &tempgvec_bi));
        tempexpnt1 = expfac * exp(-tempGmag*tempGmag / four * G->gamma_0_inv) / tempGmag / tempGmag ;
        if (tempGmag < 0.000000001) tempexpnt1 = k_one;
        if (tempexpnt1 > reciprocal_space_cutoff) {
        G->vec_bi[count].comp1 = k_zero;
        G->vec_bi[count].comp2 = k_zero;
        G->vec_bi[count].comp3 = i * crystal.reciprocal_cell[0].comp3;
        if (i > 0) {
           G->vec_b2[count_b2] = G->vec_bi[count];
	   G->mag[count_b2] = tempGmag;
           count_b2++;
          }
           Gmag_bi[count] = tempGmag;
           count++;
          }
        }
         //fprintf(file.out,"count %d %d %d\n",G->num[0],G->last_vector,G->last_vector);

         break;

    case 'M':

         break;

  } // close switch(crystal.type

    Gvec_mag = (int *) malloc(G->max_vector * sizeof(double));
    if (Gvec_mag == NULL) {
    fprintf(stderr, "ERROR: not enough memory for int Gvec_mag \n");
    exit(1);
   }

  for (i = 0; i < count_b2; i++)
  gvec_ord[i] = i;
  for (i = 1; i < count_b2; ++i) {
    for (j = count_b2 - 1; j >= i; --j) {
      if ((G->mag[j - 1]) > (G->mag[j])) {

        tempGmag = G->mag[j - 1];
        G->mag[j - 1] = G->mag[j];
        G->mag[j] = tempGmag;

        k = gvec_ord[j - 1];
        gvec_ord[j - 1] = gvec_ord[j];
        gvec_ord[j] = k;

        tempgvec_bi = G->vec_b2[j - 1];
        G->vec_b2[j - 1] = G->vec_b2[j];
        G->vec_b2[j] = tempgvec_bi;

      }
    }
  }

  for (i = 1; i < count; ++i) {
    for (j = count - 1; j >= i; --j) {
      if ((Gmag_bi[j - 1]) > (Gmag_bi[j])) {

        tempGmag = Gmag_bi[j - 1];
        Gmag_bi[j - 1] = Gmag_bi[j];
        Gmag_bi[j] = tempGmag;

        tempgvec_bi = G->vec_bi[j - 1];
        G->vec_bi[j - 1] = G->vec_bi[j];
        G->vec_bi[j] = tempgvec_bi;

      }
    }
  }

      Gvec_mag[0] = 0;
      G->sqr[0] = k_zero;
      G->EXPFAC[0] = k_one;
      G->num[0] = 1;
      G->x[0] = k_one;
      G->y[0] = k_one;
      G->z[0] = k_one;
      G_shell = 0;
      for (index_G = 1; index_G < G->last_vector; index_G++) {
        for (i = 0; i < 9; i++) {
          G->x[index_G * 9 + i] = pow(G->vec_b2[index_G].comp1,double(i)) ;
          G->y[index_G * 9 + i] = pow(G->vec_b2[index_G].comp2,double(i)) ;
          G->z[index_G * 9 + i] = pow(G->vec_b2[index_G].comp3,double(i)) ;
         }
          G->num[index_G] = 1;
          Gsq0 = double_vec_dot(&G->vec_b2[index_G - 1], &G->vec_b2[index_G - 1]);
          Gsq1 = double_vec_dot(&G->vec_b2[index_G], &G->vec_b2[index_G]);
          //fprintf(file.out,"%3d %3d  %10.4lf %10.4lf %10.4lf\n",index_G,G->last_vector,G->vec_b2[index_G].comp1,G->vec_b2[index_G].comp2,\
          G->vec_b2[index_G].comp3);
          if (fabs(Gsq1 - Gsq0) < 0.0000001) {
          G->num[G_shell] += 1;
        } else {
          G_shell++;
          G->sqr[G_shell] = Gsq1;
          Gvec_mag[G_shell] = index_G;
          for (i = 0; i < 9; i++) 
          G->shell_mag[G_shell * 9 + i] = pow(sqrt(Gsq1),double(i));
          G->invsqr[G_shell] = k_one / Gsq1;

          switch (crystal.type[0]) {

            case 'C':
            case 'P':

            expnt1 = expfac * exp(-G->sqr[G_shell] / four * G->gamma_0_inv) * G->invsqr[G_shell];
            G->EXPFAC[G_shell] = expnt1;
            //fprintf(file.out,"expnt %3d %10.4e %10.4e %10.4e\n", G_shell,expnt1,G->EXPFAC[G_shell],G->invsqr[G_shell]);

            break;

            case 'S':

            expnt1 = expfac / four * exp(-G->sqr[G_shell] / four * G->gamma_0_inv) / sqrt(G->sqr[G_shell]);
            G->EXPFAC[G_shell] = expnt1;

            // set up binomial coefficients
 
            for (m = 0; m < max; m++) {
            for (n = 0; n < max; n++) {
            G->A->a[m][n] = k_zero;
           }
           }
            for (m = 0; m < max; m++) {
            G->A->a[m][0] = k_one;
            G->A->a[m][m] = k_one;
            for (n = 0; n < m; n++) {
            G->A->a[m + 1][n + 1] = G->A->a[m][n] + G->A->a[m][n + 1];
           }
           }
            //fprintf(file.out,"A matrix \n\n");
            //for (m = 0; m < max; m++) {
            //for (n = 0; n <= m; n++) {
            //fprintf(file.out,"%8.0lf ", G->A->a[m][n]);
           //}
            //fprintf(file.out,"\n");
           //}

            // set up matrix for derivatives of erf and erfc

            for (m = 0; m < max; m++) {
            for (n = 0; n < max; n++) {
            G->B->a[m][n] = k_zero;
           }
           }
            G->B->a[1][0] = +k_one;
            for (m = 2; m < max; m++) {
            for (n = 0; n < max - 1; n++) {
            G->B->a[m][n] = (n + 1) * G->B->a[m - 1][n + 1];
            if (n > 0) G->B->a[m][n] -= two * G->B->a[m - 1][n - 1];
           } 
           } 
            //fprintf(file.out,"B matrix \n\n");
            //for (m = 0; m < max; m++) {
            //for (n = 0; n <= m; n++) {
            //fprintf(file.out,"%8.0lf ", G->B->a[m][n]);
           //}
            //fprintf(file.out,"\n");
           //}

            break;

            case 'M':

            break;

           }
           if (G->sqr[G_shell] * expnt1 < reciprocal_space_cutoff) break;
          } // close else
         } // close loop over index_G

          G->number_of_shells = G_shell;

  if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out, "\n");
    fprintf(file.out, "                               HALF-SHELLS %5d WITHIN G->last_vector       %5d magnitude %10.4f\n\n", \
    G->number_of_shells,G->last_vector,sqrt(G->sqr[G->last_vector]));
    fprintf(file.out, "                    RECIPROCAL SPACE HALF-SHELL NUMBER CLOSING MAGNITUDE EXPFAC\n");
    fprintf(file.out, "\n");
    for (i = 0; i < G->number_of_shells; i++)
    fprintf(file.out, "%28d %5d %5d %12.6lf %12.5e\n", i, G->num[i], Gvec_mag[i + 1], G->mag[Gvec_mag[i]], G->EXPFAC[i]);
    fprintf(file.out, "\n");
    fprintf(file.out, "                        RECIPROCAL SPACE VECTORS IN HALF-SHELLS %5d\n", G->last_vector);
    fprintf(file.out, "\n");
    for (i = 0; i < G->last_vector - 1; i++)
    fprintf(file.out, "%18d %12.6lf%12.6lf%12.6lf     %12.6lf\n", i, G->vec_b2[i].comp1, G->vec_b2[i].comp2, G->vec_b2[i].comp3, G->mag[i]);
    fprintf(file.out, "\n");
    fprintf(file.out, "\n");
    fprintf(file.out, "                        RECIPROCAL SPACE VECTORS IN FULL-SHELLS %5d\n", G->max_vector);
    fprintf(file.out, "\n");
    for (i = 0; i < G->max_vector; i++)
    fprintf(file.out, "%18d %12.6lf%12.6lf%12.6lf     %12.6lf\n", i, G->vec_bi[i].comp1, G->vec_bi[i].comp2, G->vec_bi[i].comp3, Gmag_bi[i]);
    fprintf(file.out, "\n");
   }

  free(Gvec_mag);

}

void generate_q_lattice(VECTOR_INT *qvec, Q_LATTICE *q_G, FERMI *fermi, RECIPROCAL_LATTICE *G, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, index_G; 
int size = 4 * job->l_max + 1;
VECTOR_DOUBLE q_vec;

  switch (crystal->type[0]) {

    case 'C':

    q_G->gamma_0_inv = G->gamma_0_inv;
    q_G->last_vector = G->last_vector;
    q_G->max_vector  = G->max_vector;
    getcart(qvec, &q_vec, fermi->is, crystal);
    //fprintf(file.out,"q_G %10.4f %10.4f %10.4f\n",q_vec.comp1,q_vec.comp2,q_vec.comp3);
    for (index_G = 0; index_G < q_G->last_vector; index_G++) {
      q_G->vec[index_G].comp1 = q_vec.comp1 + G->vec_bi[index_G].comp1;
      q_G->vec[index_G].comp2 = q_vec.comp2 + G->vec_bi[index_G].comp2;
      q_G->vec[index_G].comp3 = q_vec.comp3 + G->vec_bi[index_G].comp3;
      //fprintf(file.out,"q_G %10.4f %10.4f %10.4f\n",q_G->vec[index_G].comp1,q_G->vec[index_G].comp2,q_G->vec[index_G].comp3);
      for (i = 0; i < size; i++) {
        q_G->x[index_G * size + i] = pow(q_G->vec[index_G].comp1,double(i)) ;
        q_G->y[index_G * size + i] = pow(q_G->vec[index_G].comp2,double(i)) ;
        q_G->z[index_G * size + i] = pow(q_G->vec[index_G].comp3,double(i)) ;
       }
        q_G->sqr[index_G] = q_G->vec[index_G].comp1 * q_G->vec[index_G].comp1 + q_G->vec[index_G].comp2 * q_G->vec[index_G].comp2 + \
        q_G->vec[index_G].comp3 * q_G->vec[index_G].comp3;
       }
    q_G->x[0] = k_one;
    q_G->y[0] = k_one;
    q_G->z[0] = k_one;

    break;

    case 'S':

    q_G->gamma_0_inv = G->gamma_0_inv;
    q_G->last_vector = G->last_vector;
    q_G->max_vector  = G->max_vector;

    break;


    case 'P':

    break;


    case 'M':

    break;

   } // close switch

  if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out, "\n");
    fprintf(file.out, "\n");
    fprintf(file.out, "                        RECIPROCAL SPACE Q VECTORS IN FULL-SHELLS %5d\n", q_G->last_vector);
    fprintf(file.out, "\n");
    for (i = 0; i < q_G->last_vector; i++)
    fprintf(file.out, "%18d %12.6lf%12.6lf%12.6lf     %12.6lf\n", \
    i, q_G->vec[i].comp1, q_G->vec[i].comp2, q_G->vec[i].comp3, sqrt(q_G->sqr[i]));
    fprintf(file.out, "\n");
   }

}
