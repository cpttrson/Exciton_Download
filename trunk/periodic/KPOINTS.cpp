/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include <mpi.h>
*/
#include <cstdlib>
#include "myconstants.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "LIMITS.h"
#include "ALLOCATE_MEMORY.h"
#include "MATRIX_UTIL.h"
#include "KPOINTS.h"

using namespace std;

void knet_size(int *ksize, int is[3], CRYSTAL *crystal)

{

  switch (crystal->type[0]) {

    case 'C':
      *ksize = is[0] * is[1] * is[2] ;
      break;

    case 'S':
      *ksize = is[0] * is[1];
      break;

    case 'P':
      *ksize = is[2];
      break;

    case 'M':
      *ksize = 1;
      break;

  } // close switch (crystal->type

}

void count_k_points(KPOINT_TRAN *knet, int is[3], CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  //int i, j, k, l, p, ksize, rot_index, num_poly;
  int i, j, k, l, p, ksize, nkunique, rot_index = 0;
  int *p_inr, *taken;
  VECTOR_INT kvec[3];

  // apply symmops to recip lattice vectors, if they transform into each other and is[0] != is[1], stop

  //knet_size(ksize,is,crystal);
  knet_size(&ksize,is,crystal);
  knet->nktot = ksize;

  //taken = (int *) malloc(*ksize * sizeof(int));
  taken = (int *) malloc(ksize * sizeof(int));
  if (taken == NULL) { fprintf(stderr, "ERROR: not enough memory for int taken\n"); exit(1); }
  for (i = 0; i < ksize; i++) 
  //for (i = 0; i < *ksize; i++) 
  taken[i] = -1;

  nkunique = 0;
  //*nkunique = 0;

  switch (crystal->type[0]) {

  case 'C':

      for (i = 0; i < is[0]; i++) {
        for (j = 0; j < is[1]; j++) {
          for (k = 0; k < is[2]; k++) {
            kvec[0].comp1 = i;
            kvec[0].comp2 = j;
            kvec[0].comp3 = k;
              p = i * is[1] * is[2] + j * is[2] + k;
                if (taken[p] == -1) {
                   //(*nkunique)++;
                    nkunique++;
                    }
                     for (l = 0; l < symmetry->number_of_operators; l++) {
                     //for (l = 0; l < num_symm; l++) {
                       p_inr = symmetry->inr + symmetry->inverse[l] * 9;
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       //p_inr = symmetry->inr + l * 9;
                       //rotate_vector_int_transpose(p_inr, &kvec[0], &kvec[1]);
                       //fprintf(file.out,"kvec0 %3d %4d %4d %4d  %4d  %4d  %4d %4d %4d    %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                       //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                       //kvec[1].comp3 -= kvec[1].comp3 / is[2];
                       kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                       kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                       kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
                       //fprintf(file.out,"kvec1 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                       if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                       if (kvec[1].comp3 < 0) kvec[1].comp3 += is[2];
                       //fprintf(file.out,"kvec3 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                       if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                       if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                       rot_index = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
                       //fprintf(file.out,"kvec4 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       taken[rot_index] = 1;
                       //TIME_REVERSAL
                       if (job->trs == 1) {
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       //rotate_vector_int_transpose(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp1 *= -1;
                       kvec[1].comp2 *= -1;
                       kvec[1].comp3 *= -1;
                       //fprintf(file.out,"kvec5 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                       kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                       kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
                       //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                       //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                       //kvec[1].comp3 -= kvec[1].comp3 / is[2];
                       //fprintf(file.out,"kvec6 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                       if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                       if (kvec[1].comp3 < 0) kvec[1].comp3 += is[2];
                       //fprintf(file.out,"kvec7 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l, \
                       rot_index,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                       if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                       if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                       rot_index = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
                       //fprintf(file.out,"kvec8 %3d %4d %4d %4d  %4d  %4d   %4d %4d %4d\n",nkunique,i,j,k,l,rot_index, \
                       kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                       //fflush(file.out);
                       taken[rot_index] = 1;
                      }
                       //TIME_REVERSAL
                      }
                     }
                    }
                   }
                     knet->unique = nkunique;

      break;

  case 'S':

      for (i = 0; i < is[0]; i++) {
        for (j = 0; j < is[1]; j++) {
            kvec[0].comp1 = i;
            kvec[0].comp2 = j;
            kvec[0].comp3 = 0;
              p = i * is[1]  + j;
                if (taken[p] == -1) 
                   //(*nkunique)++;
                     nkunique++;
                     for (l = 0; l < symmetry->number_of_operators; l++) {
                     //for (l = 0; l < num_symm; l++) {
                       p_inr = symmetry->inr + symmetry->inverse[l] * 9;
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                       kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                       //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                       //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                       if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                       if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                       if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                       if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                       rot_index = kvec[1].comp1 * is[1] + kvec[1].comp2;
                       taken[rot_index] = 1;
                       //TIME_REVERSAL
                       if (job->trs == 1) {
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp1 *= -1;
                       kvec[1].comp2 *= -1;
                       kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                       kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                       //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                       //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                       if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                       if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                       if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                       if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                       rot_index = kvec[1].comp1 * is[1] + kvec[1].comp2;
                       taken[rot_index] = 1;
                      }
                       //TIME_REVERSAL
                      }
                     }
                    }
                     knet->unique = nkunique;

      break;

  case 'P':

      for (i = 0; i < is[2]; i++) {
            kvec[0].comp1 = 0;
            kvec[0].comp2 = 0;
            kvec[0].comp3 = i;
                if (taken[i] == -1) 
                   //(*nkunique)++;
                     nkunique++;
                     for (l = 0; l < symmetry->number_of_operators; l++) {
                     //for (l = 0; l < num_symm; l++) {
                       p_inr = symmetry->inr + symmetry->inverse[l] * 9;
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp3 -= (kvec[1].comp1 / is[2]) * is[2];
                       if (kvec[1].comp3 <  0) kvec[1].comp3 += is[2];
                       if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                       rot_index = kvec[1].comp3;
                       taken[rot_index] = 1;
                       //TIME_REVERSAL
                       if (job->trs == 1) {
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp3 *= -1;
                       kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
                       if (kvec[1].comp3 <  0) kvec[1].comp3 += is[2];
                       if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                       rot_index = kvec[1].comp3;
                       taken[rot_index] = 1;
                      }
                       //TIME_REVERSAL
                       //fprintf(file.out,"kvec nkunique %3d i %3d l %4d rot_index %4d   %4d %4d %4d\n",*nkunique,i,l,rot_index, \
                       kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                      } // close loop over l
                     } // close loop over i
                      knet->unique = nkunique;

      break;

  case 'M':

            //*nkunique = 1;
            knet->unique = 1;

      break;

  } // close switch (crystal->type[0]) {

   if (job->kss == 0 && crystal->type[0] != 'M') {  // switch off reciprocal space symmetry in k space only
    knet->unique = ksize;
   }

   if (job->verbosity > 1) 
    fprintf(file.out,"count_k_points nktot %d nkunique %d\n",ksize,nkunique);

    free(taken);

}

void generate_k_points(KPOINT_TRAN *knet, int is[3], CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)
//void generate_k_points(KPOINT_TRAN *knet, int num_symm, int is[3], int nkunique, CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, p, ksize, nkunique, rot_index;
  int *p_inr, *taken;
  int count;
  VECTOR_INT kvec[3];
  VECTOR_INT kvectmp;

  // apply symmops to recip lattice vectors, if they transform into each other and is[0] != is[1], stop

  knet_size(&ksize,is,crystal);

  taken = (int *) malloc(ksize * sizeof(int));
  if (taken == NULL) { fprintf(stderr, "ERROR: not enough memory for int taken\n"); exit(1); }
  for (i = 0; i < ksize; i++) 
  taken[i] = -1;

  // generate unique k points, skip if reading from Crystal09

  if (job->C09 != 1) {

  switch (crystal->type[0]) {

  case 'C':

      nkunique = 0;
      for (i = 0; i < is[0]; i++) {
        for (j = 0; j < is[1]; j++) {
          for (k = 0; k < is[2]; k++) {
            kvec[0].comp1 = i;
            kvec[0].comp2 = j;
            kvec[0].comp3 = k;
              p = i * is[1] * is[2] + j * is[2] + k;
                if (taken[p] == -1) {
                  knet->ibz[nkunique] = p;
                   nkunique++;
                    }
                     //for (l = 0; l < num_symm; l++) {
                     for (l = 0; l < symmetry->number_of_operators; l++) {
                          p_inr = symmetry->inr + symmetry->inverse[l] * 9;
                          rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                          //rotate_vector_int_transpose(p_inr, &kvec[0], &kvec[1]);
                          //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                          //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                          //kvec[1].comp3 -= kvec[1].comp3 / is[2];
                          kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                          kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                          kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
                          if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                          if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                          if (kvec[1].comp3 < 0) kvec[1].comp3 += is[2];
                          if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                          if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                          if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                          rot_index = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
                          //fprintf(file.out,"trs %3d op %3d kvec[0] %3d %3d %3d  kvec[1] %3d %3d %3d\n",\
                          0,l,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
                          taken[rot_index] = 1;
                          //TIME_REVERSAL
                          if (job->trs == 1) {
                          rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                          //rotate_vector_int_transpose(p_inr, &kvec[0], &kvec[1]);
                          kvec[1].comp1 *= -1;
                          kvec[1].comp2 *= -1;
                          kvec[1].comp3 *= -1;
                          kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                          kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                          kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
                          //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                          //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                          //kvec[1].comp3 -= kvec[1].comp3 / is[2];
                          if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                          if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                          if (kvec[1].comp3 < 0) kvec[1].comp3 += is[2];
                          if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                          if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                          if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                          rot_index = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
                          taken[rot_index] = 1;
                         }
                          //TIME_REVERSAL
                         }
                        }
                       }
                      }

   if (job->verbosity > 1) {
    fprintf(file.out,"nktot %d nkunique %d\n",ksize,nkunique);
     for (i = 0; i < nkunique; i++) {
        kvec[0].comp1 =  knet->ibz[i] / (is[2] * is[1]);
        kvec[0].comp2 = (knet->ibz[i] - kvec[0].comp1 * is[2] * is[1]) / is[2];
        kvec[0].comp3 =  knet->ibz[i] - kvec[0].comp1 * is[2] * is[1] - kvec[0].comp2 * is[2];
      fprintf(file.out, " point %d, knet %d, %d, %d multiplicity %5d\n", i, kvec[0].comp1, kvec[0].comp2, kvec[0].comp3,knet->num[i]); fflush(file.out);
     }
    }

      break;

  case 'S':

      nkunique = 0;
      for (i = 0; i < is[0]; i++) {
        for (j = 0; j < is[1]; j++) {
            kvec[0].comp1 = i;
            kvec[0].comp2 = j;
            kvec[0].comp3 = 0;
              p = i * is[1]  + j;
                if (taken[p] == -1) {
                  knet->ibz[nkunique] = p;
                   nkunique++;
                  }
                     //for (l = 0; l < num_symm; l++) {
                     for (l = 0; l < symmetry->number_of_operators; l++) {
                       p_inr = symmetry->inr + symmetry->inverse[l] * 9;
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                       kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                       //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                       //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                       if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                       if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                       if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                       if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                       rot_index = kvec[1].comp1 * is[1] + kvec[1].comp2;
                       taken[rot_index] = 1;
                       //TIME_REVERSAL
                       if (job->trs == 1) {
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp1 *= -1;
                       kvec[1].comp2 *= -1;
                       kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
                       kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
                       //kvec[1].comp1 -= kvec[1].comp1 / is[0];
                       //kvec[1].comp2 -= kvec[1].comp2 / is[1];
                       if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
                       if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
                       if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
                       if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
                       rot_index = kvec[1].comp1 * is[1] + kvec[1].comp2;
                       taken[rot_index] = 1;
                      }
                       //TIME_REVERSAL
                      }
                     }
                    }

   if (job->verbosity > 1) {
     for (i = 0; i < nkunique; i++) {
        kvec[0].comp1 = knet->ibz[i] / is[1];
        kvec[0].comp2 = knet->ibz[i] - kvec[0].comp1 * is[1];
        kvec[0].comp3 = 0;
      fprintf(file.out, " point %d, knet %d, %d, %d \n", i, kvec[0].comp1, kvec[0].comp2, kvec[0].comp3);
     }
    }

      break;

  case 'P':

      nkunique = 0;
      for (i = 0; i < is[2]; i++) {
            kvec[0].comp1 = 0;
            kvec[0].comp2 = 0;
            kvec[0].comp3 = i;
              p = i;
                if (taken[p] == -1) {
                  knet->ibz[nkunique] = p;
                   nkunique++;
                  }
                //if (taken[i] == -1) 
                   //(*nkunique)++;
                     //for (l = 0; l < num_symm; l++) {
                     for (l = 0; l < symmetry->number_of_operators; l++) {
                       p_inr = symmetry->inr + symmetry->inverse[l] * 9;
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp3 -= (kvec[1].comp1 / is[2]) * is[2];
                       if (kvec[1].comp3 <  0) kvec[1].comp3 += is[2];
                       if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                       rot_index = kvec[1].comp3;
                       taken[rot_index] = 1;
                       //TIME_REVERSAL
                       if (job->trs == 1) {
                       rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
                       kvec[1].comp3 *= -1;
                       kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
                       if (kvec[1].comp3 <  0) kvec[1].comp3 += is[2];
                       if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
                       rot_index = kvec[1].comp3;
                       taken[rot_index] = 1;
                      } 
                       //TIME_REVERSAL
                      } // close loop over l
                     } // close loop over i

          if (job->taskid == 0 && job->verbosity > 1) {
            for (i = 0; i < nkunique; i++) {
               kvec[0].comp1 = 0;
               kvec[0].comp2 = 0;
               kvec[0].comp3 = knet->ibz[i];
               fprintf(file.out, " point %d, knet %d, %d, %d \n", i, kvec[0].comp1, kvec[0].comp2, kvec[0].comp3);
              }
             }

      break;

  case 'M':

            nkunique = 1;
            knet->ibz[0] = 0;
            kvec[0].comp1 = 0;
            kvec[0].comp2 = 0;
            kvec[0].comp3 = 0;

      break;

  } // close switch (crystal->type[0]) {

           if (knet->unique != nkunique && job->kss == 1) {
           if (job->taskid == 0)
           fprintf(file.out,"Error in calculating number of unique k points in generate_k_points\n");
           MPI_Finalize();
           exit(1);
          }

  } // close if job->C09

  // generate equivalent k points from unique points

  switch (crystal->type[0]) {

  case 'C':

      for (i = 0; i < ksize; i++) 
      taken[i] = -1;

      count = 0;

      for (i = 0; i < knet->unique; i++) {
      //for (i = 0; i < nkunique; i++) {
        kvec[0].comp1 =  knet->ibz[i] / (is[2] * is[1]);
        kvec[0].comp2 = (knet->ibz[i] - kvec[0].comp1 * is[2] * is[1]) / is[2];
        kvec[0].comp3 =  knet->ibz[i] - kvec[0].comp1 * is[2] * is[1] - kvec[0].comp2 * is[2];
          //for (l = 0; l < num_symm; l++) {
          for (l = 0; l < symmetry->number_of_operators; l++) {
            p_inr = symmetry->inr + symmetry->inverse[l] * 9;
            rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
            //p_inr = symmetry->inr + l * 9;
            //fprintf(file.out,"%3d %3d %3d     %3d %3d %3d     %3d %3d %3d   %3d %3d %3d    %3d %3d %3d\n",\
            *(p_inr+0),*(p_inr+1),*(p_inr+2),*(p_inr+3),*(p_inr+4),*(p_inr+5),*(p_inr+6),*(p_inr+7),*(p_inr+8),\
            kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
            //rotate_vector_int_transpose(p_inr, &kvec[0], &kvec[1]);
            kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
            kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
            kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
            if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
            if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
            if (kvec[1].comp3 < 0) kvec[1].comp3 += is[2];
            if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
            if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
            if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
            p = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
            if (taken[p] == -1) {
            getcart(&kvec[1], &knet->cart[p], is, crystal);
            knet->oblique[p] = kvec[1];
            knet->bz[count] = p;
            knet->fbz[p] = i;
//May2013
            ////knet->opr[p] = l;
            // change 4 knet->opr[p] = symmetry->inverse[l];
            knet->opr[p] = l;
//May2013
            knet->trs[p] = 0;
            knet->num[i]++;
            knet->weight[i] = (double) knet->num[i] / (double) ksize;
            //fprintf(file.out,"kvec1 i %4d  l %4d  inv[l] %4d p %4d count %4d     %4d %4d %4d   %4d %4d %4d  num %5d ibz %5d\n",\
            i,l,symmetry->inverse[l],p,count,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3, \
            kvec[1].comp1,kvec[1].comp2,kvec[1].comp3,knet->num[i],knet->ibz[i]); fflush(file.out);
            count++;
            taken[p]++;
           }
           //TIME_REVERSAL
            if (job->trs == 1) {
            rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
            //rotate_vector_int_transpose(p_inr, &kvec[0], &kvec[1]);
            kvec[1].comp1 *= -1;
            kvec[1].comp2 *= -1;
            kvec[1].comp3 *= -1;
            kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
            kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
            kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
            if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
            if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
            if (kvec[1].comp3 < 0) kvec[1].comp3 += is[2];
            if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
            if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
            if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
            p = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
            if (taken[p] == -1) {
            getcart(&kvec[1], &knet->cart[p], is, crystal);
            knet->oblique[p] = kvec[1];
            knet->bz[count] = p;
            knet->fbz[p] = i;
//May2013
            ////knet->opr[p] = l;
            knet->opr[p] = l;
            // change 4 knet->opr[p] = symmetry->inverse[l];
//May2013
            knet->trs[p] = 1;
            knet->num[i]++;
            knet->weight[i] = (double) knet->num[i] / (double) ksize;
            count++;
            //fprintf(file.out,"kvecq %4d  %4d %4d %4d  %4d %4d %4d   %4d %4d %4d  num  %5d ibz%5d\n",\
            l,i,p,count,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,\
            kvec[1].comp1,kvec[1].comp2,kvec[1].comp3,knet->num[i],knet->ibz[i]); fflush(file.out);
            taken[p]++;
           }
           }
          //TIME_REVERSAL
          }
         }

      for (i = 0; i < ksize; i++) {
        if (taken[i] == -1) {
        if (job->taskid == 0) fprintf(file.out,"k point %4d was not initialised in generate_k_points. ksize = %5d\n", i, ksize);
        MPI_Finalize();
        exit(0);
       }
      }

      break;

  case 'S':

      for (i = 0; i < ksize; i++) 
      taken[i] = -1;

      count = 0;

      for (i = 0; i < knet->unique; i++) {
      //for (i = 0; i < nkunique; i++) {
        kvec[0].comp1 = knet->ibz[i] / is[1];
        kvec[0].comp2 = knet->ibz[i] - kvec[0].comp1 * is[1];
        kvec[0].comp3 = 0;
          //for (l = 0; l < num_symm; l++) {
          for (l = 0; l < symmetry->number_of_operators; l++) {
            p_inr = symmetry->inr + symmetry->inverse[l] * 9;
            rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
            kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
            kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
            if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
            if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
            if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
            if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
            p = kvec[1].comp1 * is[1] + kvec[1].comp2;
            if (taken[p] == -1) {
            getcart(&kvec[1], &knet->cart[p], is, crystal);
            knet->oblique[p] = kvec[1];
            knet->bz[count] = p;
            knet->fbz[p] = i;
            knet->opr[p] = l;
            knet->trs[p] = 0;
            knet->num[i]++;
            knet->weight[i] = (double) knet->num[i] / (double) ksize;
            count++;
            //fprintf(file.out,"p %3d %3d %3d %3d  %3d %3d %3d %10.4lf %10.4lf %10.4lf\n",i+1,knet->ibz[i],l,p,knet->oblique[p].comp1,knet->oblique[p].comp2,\
            knet->oblique[p].comp3,knet->cart[p].comp1,knet->cart[p].comp2,knet->cart[p].comp3);
            taken[p]++;
           }
            //TIME_REVERSAL
            if (job->trs == 1) {
            rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
            kvec[1].comp1 *= -1;
            kvec[1].comp2 *= -1;
            kvec[1].comp1 -= (kvec[1].comp1 / is[0]) * is[0];
            kvec[1].comp2 -= (kvec[1].comp2 / is[1]) * is[1];
            if (kvec[1].comp1 < 0) kvec[1].comp1 += is[0];
            if (kvec[1].comp2 < 0) kvec[1].comp2 += is[1];
            if (kvec[1].comp1 == is[0]) kvec[1].comp1 = 0;
            if (kvec[1].comp2 == is[1]) kvec[1].comp2 = 0;
            p = kvec[1].comp1 * is[1] + kvec[1].comp2;
            if (taken[p] == -1) {
            getcart(&kvec[1], &knet->cart[p], is, crystal);
            knet->oblique[p] = kvec[1];
            knet->bz[count] = p;
            knet->fbz[p] = i;
            knet->opr[p] = l;
            knet->trs[p] = 1;
            knet->num[i]++;
            knet->weight[i] = (double) knet->num[i] / (double) ksize;
            count++;
            //fprintf(file.out,"q %3d %3d %3d %3d  %3d %3d %3d %10.4lf %10.4lf %10.4lf\n",i+1,knet->ibz[i],l,p,knet->oblique[p].comp1,knet->oblique[p].comp2,\
            knet->oblique[p].comp3,knet->cart[p].comp1,knet->cart[p].comp2,knet->cart[p].comp3); fflush(file.out);
            taken[p]++;
           }
            //TIME_REVERSAL
           }
          }
         }

      break;

  case 'P':

      for (i = 0; i < ksize; i++) 
      taken[i] = -1;

      count = 0;

      for (i = 0; i < knet->unique; i++) {
      //for (i = 0; i < nkunique; i++) {
        kvec[0].comp1 = 0;
        kvec[0].comp2 = 0;
        kvec[0].comp3 = knet->ibz[i];
        //for (l = 0; l < num_symm; l++) {
        for (l = 0; l < symmetry->number_of_operators; l++) {
          p_inr = symmetry->inr + symmetry->inverse[l] * 9;
          rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
          kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
          //kvec[1].comp3 -= (kvec[1].comp1 / is[2]) * is[2];
          if (kvec[1].comp3 <  0) kvec[1].comp3 += is[2];
          if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
            p = kvec[1].comp3;
            if (taken[p] == -1) {
            getcart(&kvec[1], &knet->cart[p], is, crystal);
            knet->oblique[p] = kvec[1];
            knet->bz[count] = p;
            knet->fbz[p] = i;
            knet->opr[p] = l;
            knet->trs[p] = 0;
            knet->num[i]++;
            knet->weight[i] = (double) knet->num[i] / (double) ksize;
            //fprintf(file.out,"kvec1 %4d  %4d %4d %4d  %4d %4d %4d   %4d %4d %4d   %5d\n",l,i,p,count,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3, \
            kvec[1].comp1,kvec[1].comp2,kvec[1].comp3,knet->num[i]); fflush(file.out);
            count++;
            taken[p]++;
           }
            //TIME_REVERSAL
            if (job->trs == 1) {
            rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
            kvec[1].comp3 *= -1;
            kvec[1].comp3 -= (kvec[1].comp3 / is[2]) * is[2];
            if (kvec[1].comp3 <  0) kvec[1].comp3 += is[2];
            if (kvec[1].comp3 == is[2]) kvec[1].comp3 = 0;
            p = kvec[1].comp3;
            if (taken[p] == -1) {
            getcart(&kvec[1], &knet->cart[p], is, crystal);
            knet->oblique[p] = kvec[1];
            knet->bz[count] = p;
            knet->fbz[p] = i;
            knet->opr[p] = l;
            knet->trs[p] = 1;
            knet->num[i]++;
            knet->weight[i] = (double) knet->num[i] / (double) ksize;
            count++;
            //fprintf(file.out,"kvecq %4d  %4d %4d %4d  %4d %4d %4d   %4d %4d %4d   %5d\n",l,i,p,count,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3, \
            kvec[1].comp1,kvec[1].comp2,kvec[1].comp3,knet->num[i]); fflush(file.out);
            taken[p]++;
           }
           }
            //TIME_REVERSAL
          } // close loop over l
         } // close loop over i

      break;

  case 'M':

            nkunique = 1;
            knet->ibz[0] = 0;
            //CHANGES2017
            knet->opr[0] = 0;
            //CHANGES2017
            //CHANGES2014
            knet->num[0] = 1;
            //CHANGES2014
            kvec[0].comp1 = 0;
            kvec[0].comp2 = 0;
            kvec[0].comp3 = 0;
            //getcart(&kvec[0], &knet->cart[0], is, crystal);

      break;

  } // close switch (crystal->type[0]) {

        int knet_num[ksize];
        for (i = 0; i < knet->unique; i++)
        knet_num[i] = knet->num[i];

        if (job->kss == 0 && crystal->type[0] != 'M') {  // switch off reciprocal space symmetry in k space only
        knet->unique = ksize;
        for (i = 0; i < ksize; i++) {
        knet->ibz[i] = knet->bz[i];
        knet->weight[i] = k_one / (double) ksize;
        knet->trs[i] = 0;
        knet_num[i] = 1;
        // CHP 16/04/2014
        knet->num[i] = 1;
        knet->opr[i] = 0;
        knet->fbz[i] = i;
        knet->ibz[i] = i;
        knet->bz[i] = i;
        kvec[0].comp1 =  knet->ibz[i] / (is[2] * is[1]);
        kvec[0].comp2 = (knet->ibz[i] - kvec[0].comp1 * is[2] * is[1]) / is[2];
        kvec[0].comp3 =  knet->ibz[i] - kvec[0].comp1 * is[2] * is[1] - kvec[0].comp2 * is[2];
        p = kvec[0].comp1 * is[1] * is[2] + kvec[0].comp2 * is[2] + kvec[0].comp3;
        getcart(&kvec[0], &knet->cart[p], is, crystal);
        knet->oblique[p] = kvec[0];
       }
      }

        if (job->taskid == 0 && job->verbosity > 1 && crystal->type[0] != 'M') {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                            K POINT CALCULATION                                          |\n");
        fprintf(file.out,"===========================================================================================================\n");
        int nrows = knet->unique / 4;
        int remainder = knet->unique - nrows * 4;
        int k1, k2, k3, k4,pointer = 0;

        if (nrows > 100) { nrows = 100; remainder = 0; if (job->taskid == 0) fprintf(file.out,"Printing of large k point list truncated.\n"); }

        for (i = 0; i < nrows; i++) {
        k1 = knet->bz[pointer];
        pointer += knet_num[4 * i];
        k2 = knet->bz[pointer];
        pointer += knet_num[4 * i + 1];
        k3 = knet->bz[pointer];
        pointer += knet_num[4 * i + 2];
        k4 = knet->bz[pointer];
        pointer += knet_num[4 * i + 3];
        //fprintf(file.out,"%4d ",pointer);
        //fprintf(file.out,"%4d %4d %4d %4d\n",k1,k2,k3,k4);
        fprintf(file.out,"|   %4d %4d %4d  %7.5lf | %4d %4d %4d  %7.5lf | %4d %4d %4d  %7.5lf | %4d %4d %4d  %7.5lf |\n", \
        knet->oblique[k1].comp1,knet->oblique[k1].comp2,knet->oblique[k1].comp3,knet->weight[4 * i],\
        knet->oblique[k2].comp1,knet->oblique[k2].comp2,knet->oblique[k2].comp3,knet->weight[4 * i + 1],\
        knet->oblique[k3].comp1,knet->oblique[k3].comp2,knet->oblique[k3].comp3,knet->weight[4 * i + 2],\
        knet->oblique[k4].comp1,knet->oblique[k4].comp2,knet->oblique[k4].comp3,knet->weight[4 * i + 3]);
       }
        if (remainder == 1) {
        k1 = knet->bz[pointer];
        fprintf(file.out,"|   %4d %4d %4d  %7.5lf |                         |                         |                         |\n",\
        knet->oblique[k1].comp1,knet->oblique[k1].comp2,knet->oblique[k1].comp3,knet->weight[4 * nrows]);
       }
        if (remainder == 2) {
        k1 = knet->bz[pointer];
        pointer += knet_num[nrows * 4];
        k2 = knet->bz[pointer];
        fprintf(file.out,"|   %4d %4d %4d  %7.5lf | %4d %4d %4d  %7.5lf |                         |                         |\n",\
        knet->oblique[k1].comp1,knet->oblique[k1].comp2,knet->oblique[k1].comp3,knet->weight[4 * nrows],\
        knet->oblique[k2].comp1,knet->oblique[k2].comp2,knet->oblique[k2].comp3,knet->weight[4 * nrows + 1]);
       }
        if (remainder == 3) {
        k1 = knet->bz[pointer];
        pointer += knet_num[nrows * 4];
        k2 = knet->bz[pointer];
        pointer += knet_num[nrows * 4 + 1];
        k3 = knet->bz[pointer];
        fprintf(file.out,"|   %4d %4d %4d  %7.5lf | %4d %4d %4d  %7.5lf | %4d %4d %4d  %7.5lf |                         |\n",\
        knet->oblique[k1].comp1,knet->oblique[k1].comp2,knet->oblique[k1].comp3,knet->weight[4 * nrows],\
        knet->oblique[k2].comp1,knet->oblique[k2].comp2,knet->oblique[k2].comp3,knet->weight[4 * nrows + 1],\
        knet->oblique[k3].comp1,knet->oblique[k3].comp2,knet->oblique[k3].comp3,knet->weight[4 * nrows + 2]);
       }
        fprintf(file.out,"===========================================================================================================\n");
      }

   if (job->taskid == 0 && job->verbosity > 1) 
    { fprintf(file.out,"generate_k_points nktot %d nkunique %d\n",ksize,knet->unique); fflush(file.out); }

    free(taken);

}

void generate_k_point_pairs(KPOINT_TRAN *knet, int is[3], CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, p, q, t1, t2, ksize, unique_pairs, total_pairs, rot_index_p, rot_index_q;
int *p_inr;
VECTOR_INT kvec[3], qvec[3];
IntMatrix *taken;

  // ******************************************************************************************
  // * Count and generate k point pairs                                                       *
  // * taken array is ksize * ksize which becomes very large for a dense k mesh               *
  // * Use k_point_equivalent_pairs for large k meshes                                        *
  // ******************************************************************************************

  knet_size(&ksize,is,crystal);
  knet->nktot = ksize;

  AllocateIntMatrix(&taken,&ksize,&ksize,job);

  for (i = 0; i < ksize; i++) {
  for (j = 0; j < ksize; j++) {
  taken->a[i][j] = -1;
 }
 }

  unique_pairs = 0;
  total_pairs = 0;
  for (i = 0; i < knet->unique; i++) {
    p = knet->ibz[i];
    for (q = 0; q < knet->nktot; q++) {
      //for (t1 = 0; t1 <= job->trs; t1++) {
        kvec[0] = decompose_k_point(is, p, crystal, job, file);
        //if (t1 == 1) { kvec[0].comp1 *= -1; kvec[0].comp2 *= -1; kvec[0].comp3 *= -1; }
        for (l = 0; l < symmetry->number_of_operators; l++) {
          p_inr = symmetry->inr + symmetry->inverse[l] * 9;
          rotate_vector_int_latt(p_inr, is, t1, &kvec[0], &kvec[1], crystal, job, file);
          if (kvec[1].comp1 != kvec[0].comp1 || kvec[1].comp2 != kvec[0].comp2 || kvec[1].comp3 != kvec[0].comp3) continue;
          qvec[0] = decompose_k_point(is, q, crystal, job, file);
          if (taken->a[p][q] == -1) {
            unique_pairs++;
            //fprintf(file.out,"\n");
            rotate_vector_int_latt(p_inr, is, t1, &qvec[0], &qvec[1], crystal, job, file);
            //for (t2 = 0; t2 <= job->trs; t2++) {
              for (m = 0; m < symmetry->number_of_operators; m++) {
                p_inr = symmetry->inr + symmetry->inverse[m] * 9;
                rotate_vector_int_latt(p_inr, is, t2, &kvec[0], &kvec[2], crystal, job, file);
                rot_index_p = kvec[2].comp1 * is[1] * is[2] + kvec[2].comp2 * is[2] + kvec[2].comp3;
                rotate_vector_int_latt(p_inr, is, t2, &qvec[1], &qvec[2], crystal, job, file);
                rot_index_q = qvec[2].comp1 * is[1] * is[2] + qvec[2].comp2 * is[2] + qvec[2].comp3;
                if (taken->a[rot_index_p][rot_index_q] == -1) {
                  total_pairs++;
                  //fprintf(file.out,"i %3d p %3d q %3d l %3d m %3d un %3d to %3d\n",\
                  i,rot_index_p,rot_index_q,l,m,unique_pairs,total_pairs);
                  taken->a[rot_index_p][rot_index_q] = 1;
                 }
                } // close loop on m
            // } // close loop on t2
              } // close loop on if (taken
             }
         // } // close loop on t1
           }
          }

  fprintf(file.out,"unique pairs %3d total pairs %3d\n",unique_pairs,total_pairs);

  DestroyIntMatrix(&taken,job);

}

void k_point_equivalent_pairs(int i, KPOINT_PAIR_TRAN *kpair, KPOINT_TRAN *knet, int is[3], CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int j, k, l, m, n, p, q, t1, t2, flag, ksize, rot_index_p, rot_index_q, total_pairs;
int *p_inr;
VECTOR_INT kvec[3], qvec[3];

  // ******************************************************************************************
  // * Count and generate k point pairs for given unique k point number i                     *
  // ******************************************************************************************

knet_size(&ksize,is,crystal);
knet->nktot = ksize;

for (j = 0; j < ksize; j++) kpair->num[j]  =  0;
for (j = 0; j < ksize; j++) kpair->ibz2[j] = -1;
for (j = 0; j < ksize * symmetry->number_of_operators; j++) kpair->bz1[j] = -1;
for (j = 0; j < ksize * symmetry->number_of_operators; j++) kpair->bz2[j] = -1;

  total_pairs = 0;
  p = knet->ibz[i];
  kvec[0] = decompose_k_point(is, p, crystal, job, file);
  for (q = 0; q < knet->nktot; q++) {
    //for (t1 = 0; t1 <= job->trs; t1++) { // omit time reversal symmetry for now
    //if (t1 == 1) { kvec[0].comp1 *= -1; kvec[0].comp2 *= -1; kvec[0].comp3 *= -1; }
    for (l = 0; l < symmetry->number_of_operators; l++) {
      p_inr = symmetry->inr + symmetry->inverse[l] * 9;
      rotate_vector_int_latt(p_inr, is, t1, &kvec[0], &kvec[1], crystal, job, file);
      if (kvec[1].comp1 != kvec[0].comp1 || kvec[1].comp2 != kvec[0].comp2 || kvec[1].comp3 != kvec[0].comp3) continue;
      qvec[0] = decompose_k_point(is, q, crystal, job, file);
      rotate_vector_int_latt(p_inr, is, t1, &qvec[0], &qvec[1], crystal, job, file);
      rot_index_q = qvec[1].comp1 * is[1] * is[2] + qvec[1].comp2 * is[2] + qvec[1].comp3;
      for (k = 0; k <= total_pairs; k++) {
        if (rot_index_q == kpair->ibz2[k]) break; 
          if (k == total_pairs && (kpair->ibz2[k] != rot_index_q)) {
            kpair->ibz1[k] = p;
            kpair->ibz2[k] = rot_index_q;
            total_pairs++;
            //fprintf(file.out,"un %3d ibz1 %3d ibz2 %3d\n",total_pairs,kpair->ibz1[k],kpair->ibz2[k]);
            break;
           }
          }
         } // close loop on l
     // } // close loop on t1
       } // close loop on q

  if (total_pairs != ksize) {
  if (job->taskid == 0) fprintf(file.out,"k_point_equivalent_pairs: number of pairs of k points does not sum to ksize. Stopping.\n");
  MPI_Finalize();
  exit(1);
 }

  total_pairs = 0;
  kpair->nktot = 0;
  kpair->unique = 0;
  for (k = 0; k < ksize; k++) {
    flag = 0;
    //for (t2 = 0; t2 <= job->trs; t2++) {
    for (m = 0; m < symmetry->number_of_operators; m++) {
      p_inr = symmetry->inr + symmetry->inverse[m] * 9;
      kvec[0] = decompose_k_point(is, kpair->ibz1[k], crystal, job, file);
      //if (t2 == 1) { kvec[0].comp1 *= -1; kvec[0].comp2 *= -1; kvec[0].comp3 *= -1; }
      rotate_vector_int_latt(p_inr, is, t2, &kvec[0], &kvec[1], crystal, job, file);
      rot_index_p = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
      qvec[0] = decompose_k_point(is, kpair->ibz2[k], crystal, job, file);
      //if (t2 == 1) { qvec[0].comp1 *= -1; qvec[0].comp2 *= -1; qvec[0].comp3 *= -1; }
      rotate_vector_int_latt(p_inr, is, t2, &qvec[0], &qvec[1], crystal, job, file);
      rot_index_q = qvec[1].comp1 * is[1] * is[2] + qvec[1].comp2 * is[2] + qvec[1].comp3;
      for (n = 0; n <= total_pairs; n++) {
        if (rot_index_p == kpair->bz1[n] && rot_index_q == kpair->bz2[n]) break;
          if (n == total_pairs && (kpair->bz1[n] != rot_index_p || kpair->bz2[n] != q)) {
            flag = 1;
            kpair->bz1[n]  = rot_index_p;
            kpair->bz2[n]  = rot_index_q;
            kpair->fbz1[n] = i;
            kpair->fbz2[n] = k;
            kpair->opr[n]  = m;
           (kpair->num[kpair->unique])++;
            total_pairs++;
            fprintf(file.out,"un %3d eq %3d op %3d %3d bz1 %3d bz2 %3d fbz1 %3d fbz2 %3d\n",\
            k,n,kpair->opr[n],kpair->num[kpair->unique],kpair->bz1[n],kpair->bz2[n],kpair->fbz1[n],kpair->fbz2[n]);
            break;
           }
          } // close loop on n
         } // close loop on m
      //} // close loop on t2
        if (flag) {
          //fprintf(file.out,"\n");
          kpair->posn[kpair->unique] = kpair->nktot;
          kpair->nktot = total_pairs;
         (kpair->unique)++;
        }
       } // close loop on k

       //fprintf(file.out,"unique pairs %3d total pairs %3d\n",kpair->unique,kpair->nktot);

}

void q_point_equivalent_pairs(int i, KPOINT_PAIR_TRAN *kpair, KPOINT_TRAN *knet, int is[3], CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int j, k, l, m, n, p, q, t1, t2, flag, ksize, rot_index_k, rot_index_p, rot_index_kq, total_pairs;
int kq;
int *p_inr;
VECTOR_INT kvec[3], qvec[3], kqvec[3];

  // ******************************************************************************************
  // * Count and generate k point pairs for given unique q point number i                     *
  // ******************************************************************************************

knet_size(&ksize,is,crystal);
knet->nktot = ksize;
int taken[ksize];

for (j = 0; j < ksize; j++) taken[j]       =  -1;
for (j = 0; j < ksize; j++) kpair->num[j]  =  0;
for (j = 0; j < ksize; j++) kpair->ibz2[j] = -1;
for (j = 0; j < ksize * symmetry->number_of_operators; j++) kpair->bz1[j] = -1;
for (j = 0; j < ksize * symmetry->number_of_operators; j++) kpair->bz2[j] = -1;

  total_pairs = 0;
  kpair->nktot = 0;
  kpair->unique = 0;
  q = knet->ibz[i];
  qvec[0] = decompose_k_point(is, q, crystal, job, file);
  for (k = 0; k < knet->nktot; k++) {
    flag = 0;
    kvec[0] = decompose_k_point(is, k, crystal, job, file);
    //printf("%3d %3d %3d     %3d %3d %3d   %3d %3d %3d\n",\
    i,k,q,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,qvec[0].comp1,qvec[0].comp2,qvec[0].comp3);
    kq = compose_k_point(is, kvec[0], qvec[0], crystal, job, file);
    //kqvec[0] = decompose_k_point(is, kq, crystal, job, file);
    kqvec[0] = add_k_points(kvec[0], qvec[0]);
    if (kqvec[0].comp1 >= is[0] || kqvec[0].comp2 >= is[1] || kqvec[0].comp3 >= is[2]) continue;
    //for (t1 = 0; t1 <= job->trs; t1++) { // omit time reversal symmetry for now
    //if (t1 == 1) { kvec[0].comp1 *= -1; kvec[0].comp2 *= -1; kvec[0].comp3 *= -1; }
    t1 = 0;
    for (l = 0; l < symmetry->number_of_operators; l++) {
      p_inr = symmetry->inr + symmetry->inverse[l] * 9;
      rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
      rotate_vector_int(p_inr, &kqvec[0], &kqvec[1]);
      qvec[1] = subtract_k_points(kqvec[1], kvec[1]);
      if (qvec[1].comp1 != qvec[0].comp1 || qvec[1].comp2 != qvec[0].comp2 || qvec[1].comp3 != qvec[0].comp3) continue;
      rotate_vector_int_latt(p_inr, is, t1,  &kvec[0],  &kvec[1], crystal, job, file);
      rotate_vector_int_latt(p_inr, is, t1, &kqvec[0], &kqvec[1], crystal, job, file);
      rot_index_k  =  kvec[1].comp1 * is[1] * is[2] +  kvec[1].comp2 * is[2] +  kvec[1].comp3;
      rot_index_kq = kqvec[1].comp1 * is[1] * is[2] + kqvec[1].comp2 * is[2] + kqvec[1].comp3;
      //fprintf(file.out,"%3d    %3d %3d %3d\n",rot_index_k,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
      if (taken[rot_index_k] == -1) {
        kpair->bz1[total_pairs]  = rot_index_k;
        kpair->bz2[total_pairs]  = rot_index_kq;
        kpair->fbz1[total_pairs] = k;
        kpair->fbz2[total_pairs] = kq;
        kpair->opr[total_pairs]  = l;
       (kpair->num[kpair->unique])++;
        total_pairs++;
        taken[rot_index_k] = 1;
        flag = 1;
        //fprintf(file.out,"k %3d l %3d tot %3d  k %3d %3d %3d q  %3d %3d %3d  kq %3d %3d %3d  kp %3d %3d %3d  qp %3d %3d %3d\n",\
        k,l,total_pairs,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,qvec[0].comp1,qvec[0].comp2,qvec[0].comp3,\
        kqvec[0].comp1,kqvec[0].comp2,kqvec[0].comp3,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3,\
        qvec[1].comp1,qvec[1].comp2,qvec[1].comp3);
       }
      } // close loop on l
     // } // close loop on t1
        if (flag) {
          fprintf(file.out,"\n");
          kpair->posn[kpair->unique] = kpair->nktot;
          kpair->nktot = total_pairs;
         (kpair->unique)++;
        }
       } // close loop on k

//for (k=0;k<ksize;k++) {
  //qvec[0] = decompose_k_point(is, k, crystal, job, file);
//fprintf(file.out,"%3d %3d   %3d %3d %3d\n",k,taken[k],qvec[0].comp1,qvec[0].comp2,qvec[0].comp3);
//}
       fprintf(file.out,"unique pairs %3d total pairs %3d\n",kpair->unique,kpair->nktot);


  if (total_pairs != ksize) {
  if (job->taskid == 0) fprintf(file.out,"q_point_equivalent_pairs: number of pairs of k points does not sum to ksize. Stopping.\n");
  //MPI_Finalize();
  //exit(1);
 }

/*
  total_pairs = 0;
  kpair->nktot = 0;
  kpair->unique = 0;
  for (k = 0; k < kpair->unique; k++) {
    flag = 0;
    //for (t2 = 0; t2 <= job->trs; t2++) {
    for (l = 0; l < symmetry->number_of_operators; l++) {
      p_inr = symmetry->inr + symmetry->inverse[l] * 9;
      kvec[0]  = decompose_k_point(is, kpair->bz1[kpair->posn[k]], crystal, job, file);
      kqvec[0] = decompose_k_point(is, kpair->bz2[kpair->posn[k]], crystal, job, file);
      //if (t2 == 1) { kvec[0].comp1 *= -1; kvec[0].comp2 *= -1; kvec[0].comp3 *= -1; }
      rotate_vector_int_latt(p_inr, is, t2, &kvec[0], &kvec[1], crystal, job, file);
      rot_index_p = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
      rotate_vector_int_latt(p_inr, is, t2, &kqvec[0], &kqvec[1], crystal, job, file);
      rot_index_kq = kqvec[1].comp1 * is[1] * is[2] + kqvec[1].comp2 * is[2] + kqvec[1].comp3;

      for (n = 0; n <= total_pairs; n++) {
        if (rot_index_p == kpair->bz1[n] && rot_index_q == kpair->bz2[n]) break;
          if (n == total_pairs && (kpair->bz1[n] != rot_index_p || kpair->bz2[n] != q)) {
            flag = 1;
            kpair->bz1[n]  = rot_index_p;
            kpair->bz2[n]  = rot_index_q;
            kpair->fbz1[n] = i;
            kpair->fbz2[n] = k;
            kpair->opr[n]  = m;
           (kpair->num[kpair->unique])++;
            total_pairs++;
            fprintf(file.out,"un %3d eq %3d op %3d %3d bz1 %3d bz2 %3d fbz1 %3d fbz2 %3d\n",\
            k,n,kpair->opr[n],kpair->num[kpair->unique],kpair->bz1[n],kpair->bz2[n],kpair->fbz1[n],kpair->fbz2[n]);
            break;
           }
          } // close loop on n
         } // close loop on l
      //} // close loop on t2
        if (flag) {
          fprintf(file.out,"\n");
          kpair->posn[kpair->unique] = kpair->nktot;
          kpair->nktot = total_pairs;
         (kpair->unique)++;
        }
       } // close loop on k

*/
       //fprintf(file.out,"total pairs %3d\n",total_pairs);
       //fprintf(file.out,"unique pairs %3d total pairs %3d\n",kpair->unique,kpair->nktot);

}

void q_point_rotate_pairs(int i, int *total_pairs, int *grand_total_pairs, int *bz1, int *bz2, FERMI *fermi, KPOINT_TRAN *knet, KPOINT_TRAN *knet_little_group, int is[3], CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int j, k, l, n, q, k1, t1, kq, ksize, rot_index_k, rot_index_q;
int *p_inr;
int symmetry_number_of_operators;
double weight = k_zero;
VECTOR_INT kvec[3], qvec[3];


  // ******************************************************************************************
  // * Count and generate k point pairs for given unique q point number i                     *
  // ******************************************************************************************

  if (job->kss == 0)      symmetry_number_of_operators = 1;
  else if (job->kss == 1) symmetry_number_of_operators = symmetry->number_of_operators;
  knet_size(&ksize,is,crystal);
  for (j = 0; j < ksize * symmetry_number_of_operators; j++) bz1[j] = -1;
  for (j = 0; j < ksize * symmetry_number_of_operators; j++) bz2[j] = -1;

 *total_pairs = 0;
  q = knet->ibz[i];
  qvec[0] = decompose_k_point(is, q, crystal, job, file);
  for (k = 0; k < knet_little_group->unique; k++) {
    k1 = knet_little_group->ibz[k];
    kvec[0]  = decompose_k_point(is, k1, crystal, job, file);
//fprintf(file.out,"kvec %3d %3d %3d    %3d %3d %3d\n",kvec[0].comp1, kvec[0].comp2, kvec[0].comp3,\
qvec[0].comp1, qvec[0].comp2, qvec[0].comp3);
    for (t1 = 0; t1 <= job->trs; t1++) {
      for (l = 0; l < symmetry_number_of_operators; l++) {
        p_inr = symmetry->inr + symmetry->inverse[l] * 9;
        rotate_vector_int_latt(p_inr, is, t1, &kvec[0], &kvec[1], crystal, job, file);
        rot_index_k = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
        rotate_vector_int(p_inr, &qvec[0], &qvec[2]);
        rotate_vector_int_latt(p_inr, is, t1, &qvec[0], &qvec[1], crystal, job, file);
        rot_index_q  =  qvec[1].comp1 * is[1] * is[2] +  qvec[1].comp2 * is[2] +  qvec[1].comp3;
        kq = compose_k_point(fermi->is, kvec[1], qvec[1], crystal, job, file);
        for (n = 0; n <= *total_pairs; n++) {
          //if (rot_index_q == kpair->bz2[n]) break;
          //if (rot_index_k == kpair->bz1[n]) break;
          if (rot_index_k == bz1[n] && rot_index_q == bz2[n]) break;
          //if (rot_index_k == kpair->bz1[n] && rot_index_q == kpair->bz2[n]) break;
          //if (n == *total_pairs && (kpair->bz2[n] != rot_index_q)) {
          //if (n == *total_pairs && (kpair->bz1[n] != rot_index_k)) {
          if (n == *total_pairs && (bz1[n] != rot_index_k || bz2[n] != rot_index_q)) {
          //if (n == *total_pairs && (kpair->bz1[n] != rot_index_k || kpair->bz2[n] != rot_index_q)) {
            bz1[n]  = rot_index_k;
            bz2[n]  = rot_index_q;
//fprintf(file.out,"qvec %3d %3d %3d    %3d %3d %3d   %3d %3d %3d    k %3d k+q %3d bz1,2  %3d  %3d\n",\
qvec[0].comp1, qvec[0].comp2, qvec[0].comp3,qvec[2].comp1, qvec[2].comp2, qvec[2].comp3,\
qvec[1].comp1-qvec[2].comp1, qvec[1].comp2-qvec[2].comp2, qvec[1].comp3-qvec[2].comp3,k1,kq,bz1[n],bz2[n]);
            //fprintf(file.out,"k %3d k+q %3d bz1,2  %3d  %3d\n",k1,kq,bz1[n],bz2[n]);
            (*total_pairs)++;
            break;
           }
          } // close loop on n
         } // close loop on l
        } // close loop on t1
       weight += knet->weight[i] * knet_little_group->weight[k];
       } // close loop on k
       *grand_total_pairs += *total_pairs;
       //fprintf(file.out,"total pairs %3d grand total pairs %3d weight %10.4lf\n",*total_pairs,*grand_total_pairs,weight);

}

void print_knet(KPOINT_TRAN *knet, int is[3], CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, j, ksize;

        knet_size(&ksize,is,crystal);

        if (job->taskid == 0 && job->verbosity >= 1 && crystal->type[0] != 'M') {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                            PRINT K POINT_TRAN                                           |\n");
        fprintf(file.out,"===========================================================================================================\n");
        int nrows = knet->unique / 4;
        int remainder = knet->unique - nrows * 4;
        int k1, k2, k3, k4,count = 0;

        if (nrows > 100) { nrows = 100; remainder = 0; if (job->taskid == 0) fprintf(file.out,"Printing of large k point list truncated.\n"); }

        fprintf(file.out,"|    ibz     fbz      bz    |    opr    num    trs    |   ksize = %5d  nkunique = %5d  nktot = %5d  |\n", \
        ksize,knet->unique,knet->nktot);
        fprintf(file.out,"===========================================================================================================\n");
        for (i = 0; i < knet->unique; i++) {
          for (j = 0; j < knet->num[i]; j++) {
            k1 = knet->bz[count];
            //p = knet->oblique[count].comp1 * is[1] * is[2] + knet->oblique[count].comp2 * is[2] + knet->oblique[count].comp3;
            fprintf(file.out,"| %6d  %6d  %6d    | %6d %6d %6d    | %4d %4d %4d  %7.5lf |                         |\n", \
            knet->ibz[i],knet->fbz[k1],knet->bz[count],knet->opr[k1],knet->num[i],knet->trs[k1], \
            knet->oblique[k1].comp1,knet->oblique[k1].comp2, knet->oblique[k1].comp3,knet->weight[i]);
            count++;
           }
          }
        fprintf(file.out,"===========================================================================================================\n");
       }

}

void tetrahedron_vertices2(int is[3], int kpt, int k_vertex[6][4], KPOINT_TRAN *knet, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int k, l, m, n, l1, m1, n1;
  int k_vertop[6][4];

      switch (crystal->type[0]) {

        case 'C':

          l = kpt / (is[2] * is[1]);
          l1 = l + 1;
          if (l1 == is[0])
            l1 = 0;
          m = (kpt - l * is[2] * is[1]) / is[2];
          m1 = m + 1;
          if (m1 == is[1])
            m1 = 0;
          n = kpt - l * is[2] * is[1] - m * is[2];
          n1 = n + 1;
          if (n1 == is[2])
            n1 = 0;

          k = l * is[2] * is[1] + m * is[2] + n; // 0 0 0
          //fprintf(file.out,"k l m n %3d   %3d %3d %3d   %3d\n",k,l,m,n,knet->fbz[k]);
          k_vertex[0][0] = knet->fbz[k];
          k_vertop[0][0] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n1; // 1 1 1
          k_vertex[0][1] = knet->fbz[k];
          k_vertop[0][1] = knet->opr[k];
          k = l * is[2] * is[1] + m * is[2] + n1; // 0 0 1
          k_vertex[0][2] = knet->fbz[k];
          k_vertop[0][2] = knet->opr[k];
          k = l1 * is[2] * is[1] + m * is[2] + n1; // 1 0 1
          k_vertex[0][3] = knet->fbz[k];
          k_vertop[0][3] = knet->opr[k];

          k = l * is[2] * is[1] + m * is[2] + n; // 0 0 0
          k_vertex[1][0] = knet->fbz[k];
          k_vertop[1][0] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n1; // 1 1 1
          k_vertex[1][1] = knet->fbz[k];
          k_vertop[1][1] = knet->opr[k];
          k = l * is[2] * is[1] + m * is[2] + n1; // 0 0 1
          k_vertex[1][2] = knet->fbz[k];
          k_vertop[1][2] = knet->opr[k];
          k = l * is[2] * is[1] + m1 * is[2] + n1; // 0 1 1
          k_vertex[1][3] = knet->fbz[k];
          k_vertop[1][3] = knet->opr[k];

          k = l * is[2] * is[1] + m * is[2] + n; // 0 0 0
          k_vertex[2][0] = knet->fbz[k];
          k_vertop[2][0] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n1; // 1 1 1
          k_vertex[2][1] = knet->fbz[k];
          k_vertop[2][1] = knet->opr[k];
          k = l * is[2] * is[1] + m1 * is[2] + n; // 0 1 0
          k_vertex[2][2] = knet->fbz[k];
          k_vertop[2][2] = knet->opr[k];
          k = l * is[2] * is[1] + m1 * is[2] + n1; // 0 1 1
          k_vertex[2][3] = knet->fbz[k];
          k_vertop[2][3] = knet->opr[k];

          k = l * is[2] * is[1] + m * is[2] + n; // 0 0 0
          k_vertex[3][0] = knet->fbz[k];
          k_vertop[3][0] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n1; // 1 1 1
          k_vertex[3][1] = knet->fbz[k];
          k_vertop[3][1] = knet->opr[k];
          k = l * is[2] * is[1] + m1 * is[2] + n; // 0 1 0
          k_vertex[3][2] = knet->fbz[k];
          k_vertop[3][2] = knet->opr[k];
          k = l1 * is[2] * is[1] + m * is[2] + n1; // 1 0 1
          k_vertex[3][3] = knet->fbz[k];
          k_vertop[3][3] = knet->opr[k];

          k = l * is[2] * is[1] + m * is[2] + n; // 0 0 0
          k_vertex[4][0] = knet->fbz[k];
          k_vertop[4][0] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n1; // 1 1 1
          k_vertex[4][1] = knet->fbz[k];
          k_vertop[4][1] = knet->opr[k];
          k = l1 * is[2] * is[1] + m * is[2] + n; // 1 0 0
          k_vertex[4][2] = knet->fbz[k];
          k_vertop[4][2] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n; // 1 1 0
          k_vertex[4][3] = knet->fbz[k];
          k_vertop[4][3] = knet->opr[k];

          k = l * is[2] * is[1] + m * is[2] + n; // 0 0 0
          k_vertex[5][0] = knet->fbz[k];
          k_vertop[5][0] = knet->opr[k];
          k = l1 * is[2] * is[1] + m1 * is[2] + n1; // 1 1 1
          k_vertex[5][1] = knet->fbz[k];
          k_vertop[5][1] = knet->opr[k];
          k = l1 * is[2] * is[1] + m * is[2] + n; // 1 0 0
          k_vertex[5][2] = knet->fbz[k];
          k_vertop[5][2] = knet->opr[k];
          k = l1 * is[2] * is[1] + m * is[2] + n1; // 1 0 1
          k_vertex[5][3] = knet->fbz[k];
          k_vertop[5][3] = knet->opr[k];

          break;

        case 'S':

          l = kpt / is[1];
          l1 = l + 1;
          if (l1 == is[0])
            l1 = 0;
          m = kpt - l * is[1];
          m1 = m + 1;
          if (m1 == is[1])
            m1 = 0;

          k = l * is[1] + m; // 0 0
          // CHECK new k_vertex[0][0] = k;
          k_vertex[0][0] = knet->fbz[k];
          //fprintf(file.out,"\n%3d %3d %3d %3d %3d %3d %3d\n",kpt,k,l,l1,m,m1,knet->fbz[k]);
          k = l1 * is[1] + m; // 1 0
          k_vertex[0][1] = knet->fbz[k];
          //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %3d\n",kpt,k,l,l1,m,m1,knet->fbz[k]);
          k = l * is[1] + m1; // 0 1
          k_vertex[0][2] = knet->fbz[k];
          //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %3d\n",kpt,k,l,l1,m,m1,knet->fbz[k]);

          k = l1 * is[1] + m1; // 1 1
          k_vertex[1][0] = knet->fbz[k];
          //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %3d\n",kpt,k,l,l1,m,m1,knet->fbz[k]);
          k = l1 * is[1] + m; // 1 0
          k_vertex[1][1] = knet->fbz[k];
          //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %3d\n",kpt,k,l,l1,m,m1,knet->fbz[k]);
          k = l * is[1] + m1; // 0 1
          k_vertex[1][2] = knet->fbz[k];
          //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %3d\n",kpt,k,l,l1,m,m1,knet->fbz[k]);

          break;

        case 'P':

          l = 0;
          m = 0;
          n = kpt;
          k = n;
          k_vertex[0][0] = knet->fbz[k];
          n1 = n + 1;
          if (n1 == is[2])
            n1 = 0;
          k = n1;
          k_vertex[1][0] = knet->fbz[k];

          break;

        case 'M':

          fprintf(file.out, "Density of States not calculated for finite systems\n");
          exit(1);

          break;

      } // close switch (crystal_type

}

void polygons_vertices_state_density(int is[3], int *polygons, int *vertices, double *vol_fac, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  switch (crystal->type[0]) {

    case 'C':
      *vertices = 4 ;
      *polygons = 6 ;
      *vol_fac = job->spin_fac / au_to_eV / is[0] / is[1] / is[2] / two;

      break;

    case 'S':
      *vertices = 3 ;
      *polygons = 2 ;
      *vol_fac = job->spin_fac / au_to_eV / is[0] / is[1] / two;

      break;

    case 'P':
      *vertices = 2 ;
      *polygons = 1 ;
      break;

    case 'M':
      fprintf(file.out, "Density of States not calculated for finite systems\n");
      exit(1);
      break;

  } // close switch (crystal->type

}

void polygons_vertices_susceptibility(int is[3], int *polygons, int *vertices, double *vol_fac, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  switch (crystal->type[0]) {

    case 'C':
      *vertices = 4 ;
      *polygons = 6 ;
      *vol_fac = job->spin_fac * hbar * hbar * hbar * hbar / el / au_to_eV / au_to_eV / au_to_eV / epsilon_0 / m0 / m0 / \
      a0 / a0 / a0 / a0 / a0 / crystal->primitive_cell_volume / is[0] / is[1] / is[2] / two;

      break;

    case 'S':
      *vertices = 3 ;
      *polygons = 2 ;
      *vol_fac = job->spin_fac * hbar * hbar * hbar * hbar / el / au_to_eV / au_to_eV / au_to_eV / epsilon_0 / m0 / m0 / \
      a0 / a0 / a0 / a0 / crystal->primitive_cell_volume / is[0] / is[1] / two;

      break;

    case 'P':
      *vertices = 2 ;
      *polygons = 1 ;
      break;

    case 'M':
      fprintf(file.out, "Density of States not calculated for finite systems\n");
      exit(1);
      break;

  } // close switch (crystal->type

}

void polygons_vertices_vol(int is[3], int *polygons, int *vertices, double *vol_fac, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  VECTOR_KNET k_tmp, k_tmp1, k_tmp2, k_tmp3;

  switch (crystal->type[0]) {

    case 'C':
      *vertices = 4 ;
      *polygons = 6 ;
      k_tmp1.oblique.comp1 = 1;
      k_tmp1.oblique.comp2 = 0;
      k_tmp1.oblique.comp3 = 0;
      k_tmp2.oblique.comp1 = 1;
      k_tmp2.oblique.comp2 = 1;
      k_tmp2.oblique.comp3 = 0;
      k_tmp3.oblique.comp1 = 1;
      k_tmp3.oblique.comp2 = 1;
      k_tmp3.oblique.comp3 = 1;
      getcart(&k_tmp1.oblique, &k_tmp1.cart, is, crystal);
      getcart(&k_tmp2.oblique, &k_tmp2.cart, is, crystal);
      getcart(&k_tmp3.oblique, &k_tmp3.cart, is, crystal);
      double_vec_cross(&k_tmp2.cart, &k_tmp3.cart, &k_tmp.cart);
      //*vol_fac = fabs(double_vec_dot(&k_tmp1.cart, &k_tmp.cart)) / two / six;
      *vol_fac = job->spin_fac * hbar * hbar * hbar * hbar / el / au_to_eV / au_to_eV / au_to_eV / epsilon_0 / m0 / m0 / \
      a0 / a0 / a0 / a0 / a0 / crystal->primitive_cell_volume / is[0] / is[1] / is[2] / two;
      //fprintf(file.out,"V_two %f %f %f %e\n",k_tmp1.cart.comp1,k_tmp1.cart.comp2,k_tmp1.cart.comp3,*vol_fac) ;

      break;

    case 'S':
      *vertices = 3 ;
      *polygons = 2 ;
      k_tmp1.oblique.comp1 = 1;
      k_tmp1.oblique.comp2 = 0;
      k_tmp1.oblique.comp3 = 0;
      k_tmp2.oblique.comp1 = 0;
      k_tmp2.oblique.comp2 = 1;
      k_tmp2.oblique.comp3 = 0;
      k_tmp3.oblique.comp1 = 0;
      k_tmp3.oblique.comp2 = 0;
      k_tmp3.oblique.comp3 = 0;
      getcart(&k_tmp1.oblique, &k_tmp1.cart, is, crystal);
      getcart(&k_tmp2.oblique, &k_tmp2.cart, is, crystal);
      getcart(&k_tmp3.oblique, &k_tmp3.cart, is, crystal);
      double_vec_cross(&k_tmp1.cart, &k_tmp2.cart, &k_tmp.cart);
      //*vol_fac = sqrt(double_vec_dot(&k_tmp.cart, &k_tmp.cart));
      *vol_fac = job->spin_fac * hbar * hbar * hbar * hbar / el / au_to_eV / au_to_eV / au_to_eV / epsilon_0 / m0 / m0 / \
      a0 / a0 / a0 / a0 / crystal->primitive_cell_volume / is[0] / is[1] / two;
      //fprintf(file_out,"A_two %f %f %f %e\n",k_tmp.cart.comp1,k_tmp.cart.comp2,k_tmp.cart.comp3,A_two) ;

      break;

    case 'P':
      *vertices = 2 ;
      *polygons = 1 ;
      break;

    case 'M':
      fprintf(file.out, "Density of States not calculated for finite systems\n");
      exit(1);
      break;

  } // close switch (crystal->type

}

VECTOR_INT decompose_k_point(int is[3], int kpt, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int l, m, n;
  VECTOR_INT point;

      switch (crystal->type[0]) {

        case 'C':

          l =  kpt / (is[2] * is[1]);
          m = (kpt - l * is[2] * is[1]) / is[2];
          n =  kpt - l * is[2] * is[1] - m * is[2];

          break;

        case 'S':

          l = kpt / is[1];
          m = kpt - l * is[1];
          n = 0;

          break;

        case 'P':

          l = 0;
          m = 0;
          n = kpt;

          break;

        case 'M':

          fprintf(file.out, "ERROR: decompose_k_point k points not calculated for finite systems\n");
          exit(1);

      } // close switch (crystal_type

          point.comp1 = l;
          point.comp2 = m;
          point.comp3 = n;

    return point;
}

int compose_k_point(int is[3], VECTOR_INT point1, VECTOR_INT point2, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int l, m, n, p;

        l = point1.comp1 + point2.comp1;
        m = point1.comp2 + point2.comp2;
        n = point1.comp3 + point2.comp3;
    
        if (l >= is[0]) l -= is[0];
        if (m >= is[1]) m -= is[1];
        if (n >= is[2]) n -= is[2];

      switch (crystal->type[0]) {

        case 'C':


          p = l * is[1] * is[2] + m * is[2] + n;
if (p > is[0] * is[1] * is[2]) {fprintf(file.out,"p is %3d > is^3 \n",p); exit(1); }

          break;

        case 'S':

          p = l * is[1] + m;

          break;

        case 'P':

          p = n;

          break;

        case 'M':

          fprintf(file.out, "ERROR: compose_k_point k points not calculated for finite systems\n");
          exit(1);

      } // close switch (crystal_type

    return p;
}

VECTOR_INT add_k_points(VECTOR_INT point1, VECTOR_INT point2)

{

VECTOR_INT sum;

        sum.comp1 = point1.comp1 + point2.comp1;
        sum.comp2 = point1.comp2 + point2.comp2;
        sum.comp3 = point1.comp3 + point2.comp3;

    return sum;
}
    
VECTOR_INT subtract_k_points(VECTOR_INT point1, VECTOR_INT point2)

{

VECTOR_INT diff;

        diff.comp1 = point1.comp1 - point2.comp1;
        diff.comp2 = point1.comp2 - point2.comp2;
        diff.comp3 = point1.comp3 - point2.comp3;

    return diff;
}
    
void  occupation_fermi(FERMI *fermi, KPOINT_TRAN *knet, double *eigval, int nkunique, int nbands, JOB_PARAM *job, FILES file)

{

int i, j, k, occupied_states;
double fac, fermi_occupation;

  for (i = 0; i < job->spin_dim;i++) {
    for (j = 0; j < nkunique;j++) {
      fermi->occupied[i * nkunique + j] = 0;
        for (k = 0; k < nbands;k++) {
          fermi->occupation[i * nkunique  * nbands + j * nbands + k] = k_zero;
         }
        }
       }

  fac = k_one;
  occupied_states = 0;
  for (i = 0; i < job->spin_dim;i++) {
    for (j = 0; j < nkunique;j++) {
      for (k = 0; k < nbands;k++) {
        if (*(eigval + i * nkunique * nbands + j * nbands + k) <= job->fermi_energy) {
        //if (*(eigval + i * nkunique * nbands + j * nbands + k) < job->fermi_energy) {
          occupied_states += job->spin_fac * knet->num[j];
          if (occupied_states > fermi->target_states) fac = k_one - (double)(occupied_states - fermi->target_states) / (double) knet->num[j];
fprintf(file.out,"%d %d %d %d %12.5lf %10.4lf %12.6lf\n",i,j,k,occupied_states,job->fermi_energy,fac,*(eigval + i * nkunique * nbands + j * nbands + k));
          fermi->occupied[i * nkunique + j]++;
          fermi->occupation[i * nkunique  * nbands + j * nbands + k] += (double) job->spin_fac * fac;
          fac = k_one;
          //if (occupied_states > fermi->target_states) break;
         }
        }
       }
      }

  fermi_occupation = k_zero;
  for (i = 0; i < job->spin_dim;i++) {
    for (j = 0; j < nkunique;j++) {
      for (k = 0; k < nbands;k++) {
        fermi_occupation += fermi->occupation[i * nkunique  * nbands + j * nbands + k] * knet->num[j];
       }
      }
     }

printf("fermi occupation %10.4lf target states %5d\n",fermi_occupation,fermi->target_states);

}

void calculate_fermi_level(FERMI *fermi, double *eigval, int bands[2], KPOINT_TRAN *knet, int nkunique, int nktot, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, k;
int working_value, occupied_states;
int nbands = bands[1] - bands[0] + 1;

  job->fermi_energy = initial_fermi_level_guess;
  occupied_states = 0;

  working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) / 2 - 1;
  //fprintf(file.out,"%d %lf %d %d\n",working_value ,atoms->number_of_electrons_in_unit_cell, bands[1] ,nktot);
  //printf("%d %lf %d %d\n",working_value ,atoms->number_of_electrons_in_unit_cell, bands[1] ,nktot);

  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      if (job->fermi_energy < *(eigval + i * nkunique * nbands + j * nbands + working_value))
        job->fermi_energy = *(eigval + i * nkunique * nbands + j * nbands + working_value) + k_tol;
        //fprintf(file.out,"%e\n",job->fermi_energy);
        //printf("%e\n",job->fermi_energy);
       }
      }

  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (job->fermi_energy > *(eigval + i * nkunique * nbands + j * nbands + k))
          occupied_states += knet->num[j];
          //fprintf(file.out,"%d\n",occupied_states);
          //printf("%d\n",occupied_states);
         }
        }
       }

     fprintf(file.out,"Fermi energy = %e Occupied States = %d States per k point %e\n",job->fermi_energy,occupied_states,(double)(occupied_states/nktot));
     //printf("Fermi energy = %e Occupied States = %d\n",job->fermi_energy,occupied_states/nktot);
     //fflush(file.out);
}

void  calculate_fermi_level_metal(FERMI *fermi, double *eigval, int bands[2], KPOINT_TRAN *knet, int nkunique, int nktot, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, k, count, tmp;
int working_value, occupied_states, bracketed_states, target_states, fermi_level_offset;
int nbands = bands[1] - bands[0] + 1;
double fermi_energy_temp[job->spin_dim * nkunique];
double fermi_energy_max, fermi_energy_min, fermi_energy_mean;
double temp;

  if (job->spin_dim == 2) {
    fprintf(file.out,"Check that calculate_fermi_level_metal in KPOINTS.cpp is working for spin polarised systems\n");
    MPI_Finalize();
    exit(0);
   }

  working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) ;
  target_states = working_value * nktot;
  //working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) / 2 ;
  //target_states = 2 * working_value * nktot;
  //fprintf(file.out,"%d %lf %d %d\n",working_value, atoms->number_of_electrons_in_unit_cell, bands[1] ,nktot);
  //printf("%d %lf %d %d\n",working_value, atoms->number_of_electrons_in_unit_cell, bands[1] ,nktot);
  //if (working_value == 1)
  //working_value = 0;    // fix for 1 band test case

  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
        //fermi_energy_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + working_value);
        fermi_energy_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + working_value / 2 - 1);
        //fprintf(file.out,"fermi_energy_count %d %d %e\n",i,j,fermi_energy_temp[count]);
        //printf("%e\n",fermi_energy_temp[count]);
        count++;
       }
      }

  fermi_energy_max =  initial_fermi_level_guess;
  fermi_energy_min = -initial_fermi_level_guess;

  //printf("E_f min %12.6e E_f max %12.6e\n",fermi_energy_min,fermi_energy_max);
  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
        if (fermi_energy_temp[count] < fermi_energy_min) {
           fermi_energy_min =  fermi_energy_temp[count];
          }
        if (fermi_energy_temp[count] > fermi_energy_max) {
           fermi_energy_max =  fermi_energy_temp[count];
          }
        count++;
       }
        //fprintf(file.out,"E_f min %12.6e E_f max %12.6e\n",fermi_energy_min,fermi_energy_max);
        printf("E_f min %12.6e E_f max %12.6e\n",fermi_energy_min,fermi_energy_max);
      }

  // Move window containing fermi energy in by bisection

     count = 0;

  do {

    fermi_energy_mean = (fermi_energy_max + fermi_energy_min) / two;

    //fprintf(file.out,"E_f min %12.6e E_f max %12.6e E_f mean %12.6e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
    //printf("E_f min %12.6e E_f max %12.6e E_f mean %12.6e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
    occupied_states   = 0;
    bracketed_states  = 0;

  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_mean >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
          occupied_states += job->spin_fac * knet->num[j];
         }
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            bracketed_states += job->spin_fac * knet->num[j];
            //bracketed_states++;
           }
          }
        //fprintf(file.out,"occupied %7d bracketed %4d\n",occupied_states,bracketed_states);
        //printf("occupied %7d bracketed %4d target_states %4d j %4d knet->num %4d\n",occupied_states,bracketed_states, \
        target_states,j,knet->num[j]);
         }
        }

     if (occupied_states <= target_states) {
        fermi_energy_min = fermi_energy_mean;
       }
     else {
        fermi_energy_max = fermi_energy_mean;
       }

      count++;
    
      //printf("target %7d occupied %7d bracketed %4d %4d\n",target_states,occupied_states,bracketed_states,target_states);
      printf("E_f min %18.12e E_f max %18.12e E_f mean %18.12e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);

      if (count > 100) {
      fprintf(file.out,"Fermi energy not found in calculate_fermi_level_metal\n");
      fprintf(file.out,"target %7d occupied %7d bracketed %4d nktot %4d\n",target_states,occupied_states,bracketed_states,nktot);
      fprintf(file.out,"E_f min %18.12e E_f max %18.12e E_f mean %18.12e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
      MPI_Finalize();
      exit(1);
     }

      }  while (bracketed_states > 2 * nktot);

   // one more bisection to get final number of bracketed states in this window

//May2013
for (j = 0; j < nkunique; j++) 
fermi->occupied[j] = 0;
//May2013

  occupied_states   = 0;
  bracketed_states  = 0;
  count             = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_min >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
          occupied_states += job->spin_fac * knet->num[j];
//May2013
(fermi->occupied[j])++;
fprintf(file.out,"%d %d %d\n",j,k,fermi->occupied[j]); fflush(file.out);
//May2013
         }
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            bracketed_states += job->spin_fac * knet->num[j];
            count++;
           }
        //fprintf(file.out,"occupied %7d bracketed %4d\n",occupied_states,bracketed_states);
          }
         }
        }

//May2013
//for (j = 0; j < nkunique; j++) 
        //printf("occupied %7d bracketed %4d %4d fermi->occupied %4d\n",occupied_states,bracketed_states,target_states,fermi->occupied[j]);
//May2013
        printf("occupied %7d bracketed %4d %4d\n",occupied_states,bracketed_states,target_states);

fprintf(file.out,"%d %d\n",count,fermi->occupied[0]); fflush(file.out);
     double eigenval_temp[count];
     int    eigenval_tmp[count];

     fermi_level_offset = target_states - occupied_states;
     //fprintf(file.out,"occupied %7d bracketed %4d offset %3d\n",occupied_states,bracketed_states,fermi_level_offset);
     //if (fermi_level_offset > bracketed_states || fermi_level_offset < 0) {
     if (fermi_level_offset > bracketed_states) {
        fprintf(file.out,"Fermi level in calc_fermi_level_metal not in range \nTarget %7d Occupied %7d bracketed %4d offset %3d\n",
        target_states,occupied_states,bracketed_states,fermi_level_offset);
        MPI_Finalize();
        exit(0);
       }

  // Get eigenvalues in range and degeneracies

  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            eigenval_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + k);
            eigenval_tmp[count]  = job->spin_fac * knet->num[j];
            count++;
           }
          }
         }
        }

fprintf(file.out,"%d %d\n",count,fermi->occupied[0]); fflush(file.out);
  // Sort eigenvalues and degeneracies in ascending order

  for (i = 0; i < count; i++) {
    for (j = count - 1; j >= 1; j--) {
      if (eigenval_temp[j - 1] > eigenval_temp[j]) {

        temp = eigenval_temp[j - 1];
        eigenval_temp[j - 1] = eigenval_temp[j];
        eigenval_temp[j] = temp;

        tmp = eigenval_tmp[j - 1];
        eigenval_tmp[j - 1] = eigenval_tmp[j];
        eigenval_tmp[j] = tmp;

      } 
    } 
  } 

fprintf(file.out,"%d %d\n",count,fermi->occupied[0]); fflush(file.out);
  //for (i = 0; i < count; i++) {
    //fprintf(file.out,"sorted eigenvalues %12.6e\n",eigenval_temp[i]);
   //}

  for (i = 0; i < count; i++) {
    occupied_states += eigenval_tmp[i];
    job->fermi_energy = eigenval_temp[i];
        //fprintf(file.out,"occupied %7d %e\n",occupied_states,job->fermi_energy);
        printf("target %7d count %d occupied %7d %e\n",target_states,count,occupied_states,job->fermi_energy);
    if (occupied_states >= target_states) {
       break;
      }
    }

fprintf(file.out,"%d %d\n",count,fermi->occupied[0]); fflush(file.out);
//May2013
int fermi_occupied_states = occupied_states;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
fprintf(file.out,"%d %d %d\n",j,k,fermi->occupied[j]); fflush(file.out);
      for (k = fermi->occupied[j]; k < nbands; k++) {
        if (*(eigval + i * nkunique * nbands + j * nbands + k) <= job->fermi_energy) {
          fermi_occupied_states += job->spin_fac * knet->num[j];
fprintf(file.out,"%d %d %d\n",j,k,fermi->occupied[j]); fflush(file.out);
//(fermi->occupied[j])++;
    if (fermi_occupied_states >= target_states) {
       break;
      }
}}}}

MPI_Finalize();
exit(0);

fermi_occupied_states = 0;
for (j = 0; j < nkunique; j++) {
    for (k = 0; k < fermi->occupied[j]; k++) {
        fermi_occupied_states += job->spin_fac * knet->num[j];
        //printf("occupied %7d bracketed %4d %4d fermi->occupied %4d\n",fermi_occupied_states,bracketed_states,target_states,fermi->occupied[j]);
}}
//for (j = 0; j < nkunique; j++) 
        //fprintf(file.out,"occupied %7d bracketed %4d %4d fermi->occupied %4d\n",fermi_occupied_states,bracketed_states,target_states,fermi->occupied[j]);
//May2013

        if (job->taskid == 0) {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                          FERMI LEVEL CALCULATION                                        |\n");
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"| FERMI ENERGY     %8.5lf | OCCUPIED STATES %7d | STATES PER K-POINT %8.3lf                       |\n", \
        job->fermi_energy,occupied_states,(double)occupied_states/(double)nktot);
        fprintf(file.out,"===========================================================================================================\n");
       }

     //fprintf(file.out,"Fermi energy = %e Occupied States = %d States per k point %e\n",
     //job->fermi_energy,occupied_states,(double)occupied_states/(double)nktot);
     //fflush(file.out);
     printf("Fermi energy = %e Occupied States = %d States per k point %e\n",
     job->fermi_energy,occupied_states,(double)occupied_states/(double)nktot);
}

/*
void  calculate_fermi_level_metal_old(FERMI *fermi, double *eigval, int bands[2], KPOINT_TRAN *knet, int nkunique, int nktot, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, k, count, tmp;
int working_value, occupied_states, bracketed_states, target_states, fermi_level_offset;
int nbands = bands[1] - bands[0] + 1;
double fermi_energy_temp[job->spin_dim * nkunique];
double fermi_energy_max, fermi_energy_min, fermi_energy_mean;
double temp;

  if (job->spin_dim == 2) {
    fprintf(file.out,"Check that calculate_fermi_level_metal in KPOINTS.cpp is working for spin polarised systems\n");
    MPI_Finalize();
    exit(0);
   }

  working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) ;
  target_states = working_value * nktot;
  //working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) / 2 ;
  //target_states = 2 * working_value * nktot;
  //fprintf(file.out,"%d %lf %d %d\n",working_value, atoms->number_of_electrons_in_unit_cell, bands[1] ,nktot);
  //printf("%d %lf %d %d\n",working_value, atoms->number_of_electrons_in_unit_cell, bands[1] ,nktot);
  //if (working_value == 1)
  //working_value = 0;    // fix for 1 band test case

  fermi_energy_max = -initial_fermi_level_guess;
  fermi_energy_min =  initial_fermi_level_guess;

/ * November 2013
  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
        //fermi_energy_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + working_value);
        fermi_energy_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + working_value / 2 - 1);
        //fprintf(file.out,"fermi_energy_count %d %d %e\n",i,j,fermi_energy_temp[count]);
        //printf("%e\n",fermi_energy_temp[count]);
        count++;
       }
      }

  fermi_energy_max =  initial_fermi_level_guess;
  fermi_energy_min = -initial_fermi_level_guess;

  //printf("E_f min %12.6e E_f max %12.6e\n",fermi_energy_min,fermi_energy_max);
  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
        if (fermi_energy_temp[count] < fermi_energy_min) {
           fermi_energy_min =  fermi_energy_temp[count];
          }
        if (fermi_energy_temp[count] > fermi_energy_max) {
           fermi_energy_max =  fermi_energy_temp[count];
          }
        count++;
       }
        //fprintf(file.out,"E_f min %12.6e E_f max %12.6e\n",fermi_energy_min,fermi_energy_max);
        //printf("E_f min %12.6e E_f max %12.6e\n",fermi_energy_min,fermi_energy_max);
      }

* /

  // Move window containing fermi energy in by bisection

     count = 0;

  do {

    fermi_energy_mean = (fermi_energy_max + fermi_energy_min) / two;

    //fprintf(file.out,"E_f min %12.6e E_f max %12.6e E_f mean %12.6e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
    //printf("E_f min %12.6e E_f max %12.6e E_f mean %12.6e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
    occupied_states   = 0;
    bracketed_states  = 0;

  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_mean >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
          occupied_states += job->spin_fac * knet->num[j];
         }
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            bracketed_states += job->spin_fac * knet->num[j];
            //bracketed_states++;
           }
          }
        //fprintf(file.out,"occupied %7d bracketed %4d\n",occupied_states,bracketed_states);
        //printf("occupied %7d bracketed %4d target_states %4d j %4d knet->num %4d\n",occupied_states,bracketed_states, \
        target_states,j,knet->num[j]);
         }
        }

     if (occupied_states <= target_states) {
        fermi_energy_min = fermi_energy_mean;
       }
     else {
        fermi_energy_max = fermi_energy_mean;
       }

      count++;
    
      //printf("target %7d occupied %7d bracketed %4d %4d\n",target_states,occupied_states,bracketed_states,target_states);
      //printf("E_f min %18.12e E_f max %18.12e E_f mean %18.12e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);

      if (count > 100) {
      fprintf(file.out,"Fermi energy not found in calculate_fermi_level_metal\n");
      fprintf(file.out,"target %7d occupied %7d bracketed %4d nktot %4d\n",target_states,occupied_states,bracketed_states,nktot);
      fprintf(file.out,"E_f min %18.12e E_f max %18.12e E_f mean %18.12e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
      MPI_Finalize();
      exit(1);
     }

      }  while (bracketed_states > 2 * nktot);

   // one more bisection to get final number of bracketed states in this window

  occupied_states   = 0;
  bracketed_states  = 0;
  count             = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_min >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
          occupied_states += job->spin_fac * knet->num[j];
         }
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            bracketed_states += job->spin_fac * knet->num[j];
            count++;
           }
        //fprintf(file.out,"occupied %7d bracketed %4d\n",occupied_states,bracketed_states);
          }
         }
        }

        //printf("occupied %7d bracketed %4d %4d\n",occupied_states,bracketed_states,target_states);

     double eigenval_temp[count];
     int    eigenval_tmp[count];

     fermi_level_offset = target_states - occupied_states;
     //fprintf(file.out,"occupied %7d bracketed %4d offset %3d\n",occupied_states,bracketed_states,fermi_level_offset);
     //if (fermi_level_offset > bracketed_states || fermi_level_offset < 0) {
     if (fermi_level_offset > bracketed_states) {
        fprintf(file.out,"Fermi level in calc_fermi_level_metal not in range \nTarget %7d Occupied %7d bracketed %4d offset %3d\n",
        target_states,occupied_states,bracketed_states,fermi_level_offset);
        MPI_Finalize();
        exit(0);
       }

  // Get eigenvalues in range and degeneracies

  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            eigenval_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + k);
            eigenval_tmp[count]  = job->spin_fac * knet->num[j];
            count++;
           }
          }
         }
        }

  // Sort eigenvalues and degeneracies in ascending order

  for (i = 0; i < count; i++) {
    for (j = count - 1; j >= 1; j--) {
      if (eigenval_temp[j - 1] > eigenval_temp[j]) {

        temp = eigenval_temp[j - 1];
        eigenval_temp[j - 1] = eigenval_temp[j];
        eigenval_temp[j] = temp;

        tmp = eigenval_tmp[j - 1];
        eigenval_tmp[j - 1] = eigenval_tmp[j];
        eigenval_tmp[j] = tmp;

      } 
    } 
  } 

  //for (i = 0; i < count; i++) {
    //fprintf(file.out,"sorted eigenvalues %12.6e\n",eigenval_temp[i]);
   //}

  for (i = 0; i < count; i++) {
    occupied_states += eigenval_tmp[i];
    job->fermi_energy = eigenval_temp[i];
        //fprintf(file.out,"occupied %7d %e\n",occupied_states,job->fermi_energy);
        //printf("target %7d count %d occupied %7d %e\n",target_states,count,occupied_states,job->fermi_energy);
    if (occupied_states >= target_states) {
       break;
      }
    }

        if (job->taskid == 0) {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                          FERMI LEVEL CALCULATION                                        |\n");
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"| FERMI ENERGY     %8.5lf | OCCUPIED STATES %7d | STATES PER K-POINT %8.3lf                       |\n", \
        job->fermi_energy,occupied_states,(double)occupied_states/(double)nktot);
        fprintf(file.out,"===========================================================================================================\n");
       }

     //fprintf(file.out,"Fermi energy = %e Occupied States = %d States per k point %e\n",
     //job->fermi_energy,occupied_states,(double)occupied_states/(double)nktot);
     //fflush(file.out);
     printf("Fermi energy = %e Occupied States = %d States per k point %e\n",
     job->fermi_energy,occupied_states,(double)occupied_states/(double)nktot);
}
*/

void  calculate_fermi_level_metal_old(FERMI *fermi, double *eigval, int bands[2], KPOINT_TRAN *knet, int nkunique, int nktot, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, k, m, count, tmp, lumo;
int working_value, occupied_states, bracketed_states, target_states, fermi_level_offset;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = bands[1] - bands[0] + 1;
double fermi_energy_temp[job->spin_dim * nkunique];
double fermi_energy_max, fermi_energy_min, fermi_energy_mean;
double temp;

  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < fermi->nkunique; j++) {
    fermi->occupied[i * fermi->nkunique + j] = 0;
      for (k = 0; k < dim1; k++) {
        fermi->occupation[i * fermi->nkunique * dim1 + j * dim1 + k] = k_zero;
       }
      }
     }
   
  working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) ;
  target_states = working_value * nktot;
  fermi->target_states = target_states;
  //working_value = ((int) (atoms->number_of_electrons_in_unit_cell) - 2 * bands[0] + 2) / 2 ;
  //target_states = 2 * working_value * nktot;
  //fprintf(file.out,"%d %lf %d %d\n",working_value, atoms->number_of_electrons_in_unit_cell, bands[0] ,nktot);
  //printf("%d %lf %d %d\n",working_value, atoms->number_of_electrons_in_unit_cell, bands[0] ,nktot);
  //if (working_value == 1)
  //working_value = 0;    // fix for 1 band test case

  fermi_energy_max = -initial_fermi_level_guess;
  fermi_energy_min =  initial_fermi_level_guess;

  // Move window containing fermi energy in by bisection

     count = 0;
     int bracket_target = 4 * nktot;
     //if (bracket_target > nbands) bracket_target = target_states;
     //bracket_target = 12;

  do {

    fermi_energy_mean = (fermi_energy_max + fermi_energy_min) / two;

    //fprintf(file.out,"E_f min %12.6e E_f max %12.6e E_f mean %12.6e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
    //printf("E_f min %12.6e E_f max %12.6e E_f mean %12.6e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
    occupied_states   = 0;
    bracketed_states  = 0;

  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        //CHANGES2014if (fermi_energy_min >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
        if (fermi_energy_mean > *(eigval + i * nkunique * nbands + j * nbands + k)) {
          occupied_states += job->spin_fac * knet->num[j];
         }
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            bracketed_states += job->spin_fac * knet->num[j];
            //bracketed_states++;
           }
          }
        //fprintf(file.out,"occupied %7d bracketed %4d\n",occupied_states,bracketed_states);
        //printf("occupied %7d bracketed %4d target_states %4d j %4d knet->num %4d\n",occupied_states,bracketed_states, \
        target_states,j,knet->num[j]);
         }
        }

     if (occupied_states <= target_states) {
     ////if (occupied_states < target_states) {
        fermi_energy_min = fermi_energy_mean;
       }
     else {
        fermi_energy_max = fermi_energy_mean;
       }

      count++;
    
      //fprintf(file.out,"target %7d occupied %7d bracketed %4d %4d\n",target_states,occupied_states,bracketed_states,target_states);
      //fprintf(file.out,"E_f min %18.12e E_f max %18.12e E_f mean %18.12e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);

      if (count > 100) {
      fprintf(file.out,"Fermi energy not found in calculate_fermi_level_metal\n");
      fprintf(file.out,"target %7d occupied %7d bracketed %4d nktot %4d\n",target_states,occupied_states,bracketed_states,nktot);
      fprintf(file.out,"E_f min %18.12e E_f max %18.12e E_f mean %18.12e\n",fermi_energy_min,fermi_energy_max,fermi_energy_mean);
      MPI_Finalize();
      exit(1);
     }
//fprintf(file.out,"%d %d\n",bracketed_states,bracket_target);

      //}  while (bracketed_states > bracket_target);
      }  while (bracketed_states > 4 * nktot);

   // one more bisection to get final number of bracketed states in this window
  
  occupied_states   = 0;
  bracketed_states  = 0;
  count             = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        //CHANGES2014if (fermi_energy_min >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
        if (fermi_energy_min > *(eigval + i * nkunique * nbands + j * nbands + k)) {
          occupied_states += job->spin_fac * knet->num[j];
          fermi->occupied[i * nkunique + j]++;
          fermi->occupation[i * nkunique * dim1 + j * dim1 + k] += k_one;
          job->fermi_energy = *(eigval + i * nkunique * nbands + j * nbands + k);
          job->lumo_energy  = 9999.0;
          job->fermi_k = knet->bz[j];
          job->lumo_k  = 0;
          //fprintf(file.out,"%d %10.4lf %d %10.4lf  %d %d %d\n",job->fermi_k,job->fermi_energy,job->lumo_energy,job->lumo_k, \
          knet->oblique[job->fermi_k].comp1,knet->oblique[job->fermi_k].comp2,knet->oblique[job->fermi_k].comp3);
          //CHANGES2014fermi->occupation[i * nkunique * nbands + j * nbands + k] += (double) job->spin_fac;
          //fprintf(file.out,"%d %d %d %6d %e %12.4lf %12.4lf\n",i,j,k,occupied_states,*(eigval + i * nkunique * nbands + j * nbands + k), \
          fermi_energy_min,fermi_energy_max);
         }
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            bracketed_states += job->spin_fac * knet->num[j];
            count++;
           }
          //fprintf(file.out,"Occupied %7d bracketed %4d Target %4d Fermi %10.4lf\n",occupied_states,bracketed_states,target_states,job->fermi_energy);
          }
         }
        }

        //printf("occupied %7d bracketed %4d %4d\n",occupied_states,bracketed_states,target_states);

     double eigenval_temp[count];
     int    eigenval_tmp[count];
     int    eigenval_spin[count];
     int    eigenval_kpnt[count];
     int    eigenval_band[count];

     fermi_level_offset = target_states - occupied_states;
     //fprintf(file.out,"occupied %7d bracketed %4d offset %3d\n",occupied_states,bracketed_states,fermi_level_offset);
     //if (fermi_level_offset > bracketed_states || fermi_level_offset < 0) {
     if (fermi_level_offset > bracketed_states) {
        fprintf(file.out,"Fermi level in calc_fermi_level_metal not in range \nTarget %7d Occupied %7d bracketed %4d offset %3d\n",
        target_states,occupied_states,bracketed_states,fermi_level_offset);
        MPI_Finalize();
        exit(0);
       }

  // Get eigenvalues in range and degeneracies

  count = 0;
  for (i = 0; i < job->spin_dim; i++) {
    for (j = 0; j < nkunique; j++) {
      for (k = 0; k < nbands; k++) {
        if (fermi_energy_min <= *(eigval + i * nkunique * nbands + j * nbands + k) && 
            fermi_energy_max >= *(eigval + i * nkunique * nbands + j * nbands + k)) {
            eigenval_temp[count] = *(eigval + i * nkunique * nbands + j * nbands + k);
            eigenval_tmp[count]  = job->spin_fac * knet->num[j];
            eigenval_spin[count]  = i;
            eigenval_kpnt[count]  = j;
            eigenval_band[count]  = k;
            count++;
           }
          }
         }
        }

  // Sort eigenvalues and degeneracies in ascending order

  for (i = 0; i < count; i++) {
    for (j = count - 1; j >= 1; j--) {
      if (eigenval_temp[j - 1] > eigenval_temp[j]) {

        temp = eigenval_temp[j - 1];
        eigenval_temp[j - 1] = eigenval_temp[j];
        eigenval_temp[j] = temp;

        tmp = eigenval_tmp[j - 1];
        eigenval_tmp[j - 1] = eigenval_tmp[j];
        eigenval_tmp[j] = tmp;

        tmp = eigenval_spin[j - 1];
        eigenval_spin[j - 1] = eigenval_spin[j];
        eigenval_spin[j] = tmp;

        tmp = eigenval_kpnt[j - 1];
        eigenval_kpnt[j - 1] = eigenval_kpnt[j];
        eigenval_kpnt[j] = tmp;

        tmp = eigenval_band[j - 1];
        eigenval_band[j - 1] = eigenval_band[j];
        eigenval_band[j] = tmp;

      } 
    } 
  } 

  //int sum = 0;
  //for (i = 0; i < count; i++) {
    //sum += eigenval_tmp[i];
    //fprintf(file.out,"sorted eigenvalues %12.6e %5d %5d   spin %5d kpoint %5d \n",eigenval_temp[i],eigenval_tmp[i],sum,eigenval_spin[i],eigenval_kpnt[i]);
   //}

    //CHANGES2014
    //job->fermi_energy = fermi_energy_min;
    //printf("%10.4lf\n",job->fermi_energy);

  lumo = 0;
  if (fermi_level_offset > 0) { // CHANGES2014
  for (m = 0; m < count; m++) {
    //CHANGES2014job->fermi_energy = eigenval_temp[i];
    //fprintf(file.out,"occupieda %7d %e\n",occupied_states,job->fermi_energy);
    if (occupied_states >= target_states) {
       break;
      }
       //fprintf(file.out,"%d %10.4lf %d %10.4lf  %d %d %d\n",job->fermi_k,job->fermi_energy,job->lumo_energy,job->lumo_k, \
       knet->oblique[job->fermi_k].comp1,knet->oblique[job->fermi_k].comp2,knet->oblique[job->fermi_k].comp3);
       occupied_states += eigenval_tmp[m];
       i = eigenval_spin[m];
       j = eigenval_kpnt[m];
       k = eigenval_band[m];
       double fac = k_one;
       if (occupied_states > fermi->target_states) fac = k_one - (double)(occupied_states - fermi->target_states) / (double) job->spin_fac / \
       (double) knet->num[j];
       //CHANGES2014if (occupied_states > fermi->target_states) fac = k_one - (double)(occupied_states - fermi->target_states) / (double) knet->num[j];
       //if (occupied_states > fermi->target_states) printf("%d %d %d %10.4lf %d\n",i,j,k,fac,knet->num[j]);
       fermi->occupied[i * nkunique + j]++;
       fermi->occupation[i * nkunique * dim1 + j * dim1 + k] += fac;
       //fermi->occupation[i * nkunique * nbands + j * nbands + k] += (double) job->spin_fac * fac;
       job->fermi_energy = eigenval_temp[m];
       lumo = 1;
       job->fermi_k = knet->bz[eigenval_kpnt[m]];
       if (occupied_states == fermi->target_states && m < count) { // fermi level falls on non-degenerate k point
       job->lumo_energy = eigenval_temp[m + 1];
       job->lumo_k  = knet->bz[eigenval_kpnt[m + 1]];
      }
       else {                                                      // fermi level falls on degenerate k point
       job->lumo_energy = eigenval_temp[m];
       job->lumo_k  = knet->bz[eigenval_kpnt[m]];
      }

      }
     } // CHANGES2014

     double fermi_occupation = k_zero;
     for (i = 0; i < job->spin_dim;i++) {
       for (j = 0; j < nkunique;j++) {
         for (k = 0; k < nbands;k++) {
           fermi_occupation += fermi->occupation[i * nkunique * dim1 + j * dim1 + k] * (double)job->spin_fac * knet->num[j];
           //fermi_occupation += fermi->occupation[i * nkunique * nbands + j * nbands + k] * knet->num[j];
           //fprintf(file.out,"%d %d %d  %e\n",i,j,k,fermi->occupation[i * nkunique * dim1 + j * dim1 + k]);
          }
          //fprintf(file.out,"%3d %3d   %3d\n",i,j,fermi->occupied[i * nkunique + j]);
         }
        }

        if (job->taskid == 0) {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                                          FERMI LEVEL CALCULATION                                        |\n");
        fprintf(file.out,"===========================================================================================================\n");
        if (lumo == 1) {
        fprintf(file.out,"| TARGET STATES     %7d | STATES FOUND %10.4lf | FERMI K POINT  %2d %2d %2d | LUMO K POINT   %2d %2d %2d |\n", \
        target_states,fermi_occupation,knet->oblique[job->fermi_k].comp1,knet->oblique[job->fermi_k].comp2,knet->oblique[job->fermi_k].comp3, \
        knet->oblique[job->lumo_k].comp1,knet->oblique[job->lumo_k].comp2,knet->oblique[job->lumo_k].comp3);
        fprintf(file.out,"| FERMI ENERGY (eV) %7.2lf | LUMO ENERGY (eV) %6.2lf | ENERGY GAP (eV) %7.3lf |                         |\n", \
        job->fermi_energy * au_to_eV,job->lumo_energy * au_to_eV,(job->fermi_energy - job->lumo_energy) * au_to_eV);
        //job->fermi_energy ,job->lumo_energy * au_to_eV,(job->fermi_energy - job->lumo_energy) * au_to_eV);
       }
        else {
        fprintf(file.out,"| TARGET STATES     %7d | STATES FOUND %10.4lf | FERMI K POINT  %2d %2d %2d |                         |\n", \
        target_states,fermi_occupation,knet->oblique[job->fermi_k].comp1,knet->oblique[job->fermi_k].comp2,knet->oblique[job->fermi_k].comp3);
        fprintf(file.out,"| FERMI ENERGY (eV) %7.2lf |                         |                         |                         |\n", \
        job->fermi_energy);
       }
        fprintf(file.out,"===========================================================================================================\n");
       }

}
