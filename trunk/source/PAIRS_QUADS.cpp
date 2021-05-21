/*
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <mpi.h>
#include <sys/types.h>
*/
#include "myconstants.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PAIRS_QUADS.h"

using namespace std;

void count_pairs4(PAIR_TRAN *pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Count pairs and find value of R_tables->last_vector                                   *
  // ******************************************************************************************

  int count;
  int i, j, k, l, m;
  double Rsqrd;
  VECTOR_DOUBLE Rvec_tmp;

    count = 0;
    R->last_vector = 0;
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
      for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
        for (k = 0; k < R->max_vector; k++) {
          Rvec_tmp.comp1 = atoms->cell_vector[i].comp1 - atoms->cell_vector[j].comp1 - R->vec_ai[k].comp1;
          Rvec_tmp.comp2 = atoms->cell_vector[i].comp2 - atoms->cell_vector[j].comp2 - R->vec_ai[k].comp2;
          Rvec_tmp.comp3 = atoms->cell_vector[i].comp3 - atoms->cell_vector[j].comp3 - R->vec_ai[k].comp3;
          Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
          if (Rsqrd < R->cutoff * R->cutoff) {
          count++;
          if (k > R->last_vector) R->last_vector = k;
          //fprintf(file.out,"pair  at1 at2 lat2 %3d %3d  %3d count last %4d  %3d cutoff Rsq  %10.4lf %10.4lf %18.12lf %18.12lf\n",\
          i,j,k,count,R->last_vector,R->cutoff,sqrt(Rsqrd),Rsqrd,R->cutoff * R->cutoff);
         }
        }
       }
      }

    pair_p->tot = count;
    R->last_vector++;

    R->margin_vector = 0;
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
      for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          for (l = 0; l < R->last_vector; l++) {
            Rvec_tmp.comp1 = R->vec_ai[atom_p->O[i * symmetry->number_of_operators + k]].comp1 + \
            R->vec_ai[atom_p->O[j * symmetry->number_of_operators + k]].comp1 + R->vec_ai[l].comp1;
            Rvec_tmp.comp2 = R->vec_ai[atom_p->O[i * symmetry->number_of_operators + k]].comp2 + \
            R->vec_ai[atom_p->O[j * symmetry->number_of_operators + k]].comp2 + R->vec_ai[l].comp2;
            Rvec_tmp.comp3 = R->vec_ai[atom_p->O[i * symmetry->number_of_operators + k]].comp3 + \
            R->vec_ai[atom_p->O[j * symmetry->number_of_operators + k]].comp3 + R->vec_ai[l].comp3;
            for (m = 0; m < R->max_vector; m++) { 
            if (check_vec(&Rvec_tmp, &R->vec_ai[m]) == 1) {
            //fprintf(file.out,"m %3d\n",m);
            if (m > R->margin_vector) R->margin_vector = m;
            break;
           }
          if (m == R->max_vector) {
          if (job->taskid == 0) fprintf(file.out,"R->margin vector not found in count_pairs4\n");
          MPI_Finalize();
          exit(1);
         }
        }
       }
      }
     }
    }

    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
      for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          for (l = 0; l < R->last_vector; l++) {
            Rvec_tmp.comp1 = R->vec_ai[atom_p->O[i * symmetry->number_of_operators + k]].comp1 - \
            R->vec_ai[atom_p->O[j * symmetry->number_of_operators + k]].comp1 + R->vec_ai[l].comp1;
            Rvec_tmp.comp2 = R->vec_ai[atom_p->O[i * symmetry->number_of_operators + k]].comp2 - \
            R->vec_ai[atom_p->O[j * symmetry->number_of_operators + k]].comp2 + R->vec_ai[l].comp2;
            Rvec_tmp.comp3 = R->vec_ai[atom_p->O[i * symmetry->number_of_operators + k]].comp3 - \
            R->vec_ai[atom_p->O[j * symmetry->number_of_operators + k]].comp3 + R->vec_ai[l].comp3;
            for (m = 0; m < R->max_vector; m++) { 
            if (check_vec(&Rvec_tmp, &R->vec_ai[m]) == 1) {
            //fprintf(file.out,"m %3d\n",m);
            if (m > R->margin_vector) R->margin_vector = m;
            break;
           }
          if (m == R->max_vector) {
          if (job->taskid == 0) fprintf(file.out,"R->margin vector not found in count_pairs4\n");
          MPI_Finalize();
          exit(1);
         }
        }
       }
      }
     }
    }

    R->margin_vector++;

    //fprintf(file.out,"max %4d margin %4d last %4d ewald %4d cutoff %10.4lf pair_p->tot %4d\n",\
    R->max_vector,R->margin_vector,R->last_vector,R->last_ewald_vector,R->cutoff,pair_p->tot);
    if (job->taskid == 0)
    printf("max %4d margin %4d last %4d ewald %4d cutoff %10.4lf pair_p->tot %4d\n",\
    R->max_vector,R->margin_vector,R->last_vector,R->last_ewald_vector,R->cutoff,pair_p->tot);

}

void select_pairs4(PAIR_TRAN *pair_p, ATOM *atoms, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Count pairs and find value of R_tables->last_vector                                   *
  // ******************************************************************************************

  int count;
  int i, j, k;
  double Rsqrd;
  VECTOR_DOUBLE Rvec_tmp;

    count = 0;
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
      for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
        for (k = 0; k < R->max_vector; k++) {
          Rvec_tmp.comp1 = atoms->cell_vector[i].comp1 - atoms->cell_vector[j].comp1 - R->vec_ai[k].comp1;
          Rvec_tmp.comp2 = atoms->cell_vector[i].comp2 - atoms->cell_vector[j].comp2 - R->vec_ai[k].comp2;
          Rvec_tmp.comp3 = atoms->cell_vector[i].comp3 - atoms->cell_vector[j].comp3 - R->vec_ai[k].comp3;
          Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
          if (Rsqrd < R->cutoff * R->cutoff) {
          pair_p->cell1[count] = i;
          pair_p->cell2[count] = j;
          pair_p->latt1[count] = 0;
          pair_p->latt2[count] = k;
          count++;
         }
        }
       }
      }

    pair_p->tot = count;
    pair_p->nump = count;

}

void count_density_pairs4(PAIR_TRAN *pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm, number_of_permutations;
int atm_n1, atm_n2, lat_n1, lat_n2;
int cell1_temp, cell2_temp, O1_temp, O2_temp;
int latt1_temp, latt2_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int keeper;
int pair_index_0, pair_index_1;
int total_pairs;
PAIR_TMP pairs;

  double Rsqrd;
  VECTOR_DOUBLE Rvec_tmp;

    pair_p->nump = 0;
    pair_p->tot = 0;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.cell1[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.cell2[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.latt1[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.latt2[i] = -1;

   for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
     for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
       for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
         Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
         Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
         Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
         Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);

         if (Rsqrd < R->cutoff * R->cutoff) {

         lat_n1 = 0;
         keeper = 1;
         total_pairs = 0;
         pair_index_0 = lat_n2 * dim2 + atm_n1 * dim1 + atm_n2;
         //pair_index_0 = atm_n1 * dim3 + atm_n2 * dim4 + lat_n2;
         //fprintf(file.out,"input %3d %3d %3d   %5d\n",atm_n1,atm_n2,lat_n2,pair_index_0); fflush(file.out);

         for (pm = 0; pm < number_of_permutations; pm++) {
           for (k = 0; k < symmetry->number_of_operators; k++) {

             cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
             cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

             O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
             O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
             latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector +O1_temp];
             latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector +O2_temp];

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

            pair_index_1 = rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2;
            //pair_index_1 = rotated_cell1 * dim3 + rotated_cell2 * dim4 + rotated_latt2;
            //fprintf(file.out,"%3d %3d output %3d %3d %3d   %5d          %3d %3d\n",\
            latt1_temp,latt2_temp,rotated_cell1,rotated_cell2,rotated_latt2,pair_index_1,pm,k); fflush(file.out);

            if (pair_index_1 <  pair_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }
            //if (pair_index_1 <  pair_index_0) { keeper = 0; pm = 2; k = symmetry->number_of_operators; }

         } // close loop over k
        } // close loop over pm

            if (!keeper) continue;

         for (pm = 0; pm < number_of_permutations; pm++) {
           for (k = 0; k < symmetry->number_of_operators; k++) {

             cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
             cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

             O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
             O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

             latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector +O1_temp];
             latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector +O2_temp];

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d %3d   %3d %3d\n",pm,k,pair_index_0,pair_index_1,\
           rotated_cell1,rotated_cell2,rotated_latt1,rotated_latt2,pair_p->tot,total_pairs);
  
      for (j = 0; j <= total_pairs; j++) {
        if (rotated_cell1 == pairs.cell1[j] && rotated_cell2 == pairs.cell2[j] && rotated_latt1 == pairs.latt1[j] && rotated_latt2 == pairs.latt2[j])
            break;
        if (j == total_pairs && rotated_latt2 < R_tables->last_vector && \
           (rotated_cell1 != pairs.cell1[j] || rotated_cell2 != pairs.cell2[j] || rotated_latt1 != pairs.latt1[j] || rotated_latt2 != pairs.latt2[j])) {
            pairs.cell1[j] = rotated_cell1;
            pairs.cell2[j] = rotated_cell2;
            pairs.latt1[j] = rotated_latt1;
            pairs.latt2[j] = rotated_latt2;
            total_pairs++;
            //fprintf(file.out,"pair_p %3d %3d %3d %3d %3d   %3d %3d  %3d %3d\n",j,pair_p->nump,total_pairs,pair_p->tot,rotated_cell1,rotated_cell2,\
            rotated_latt1,rotated_latt2); fflush(file.out);
            break;
           }
          } // close loop over j
         } // close loop over k
        } // close loop over pm
       pair_p->tot += total_pairs;
      (pair_p->nump)++;

      } // close if Rsqrd
      } // close loop over lat_n2
      } // close loop over atm_n2
      } // close loop over atm_n1

     //fprintf(file.out,"unique_pairs %3d total pairs %3d\n",pair_p->nump,pair_p->tot);
     //printf("unique_pairs %3d total pairs %3d\n",pair_p->nump,pair_p->tot);
       
}

void count_range_selected_pairs(PAIR_TRAN *pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm, number_of_permutations;
int atm_n1, atm_n2, lat_n1, lat_n2;
int cell1_temp, cell2_temp, O1_temp, O2_temp;
int latt1_temp, latt2_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int keeper;
int pair_index_0, pair_index_1;
int total_pairs;
double Rsqrd;
VECTOR_DOUBLE Rvec_tmp;
PAIR_TMP pairs;

if (pair_p->cutoff > R->cutoff) {
if (job->taskid == 0) fprintf(file.out,"Pairs with cutoff %7.1f greater than R->cutoff %7.1f requested. Increase R->cutoff\n",\
pair_p->cutoff, R->cutoff);
MPI_Finalize();
exit(0);
}

pair_p->nump = 0;
pair_p->tot = 0;

if (job->pms == 0) number_of_permutations = 1;
else               number_of_permutations = 2;

for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.cell1[i] = -1;
for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.cell2[i] = -1;
for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.latt1[i] = -1;
for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) pairs.latt2[i] = -1;

for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
  for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
    for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
      Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
      Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
      Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
      Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);

      if (Rsqrd < pair_p->cutoff * pair_p->cutoff) {

      lat_n1 = 0;
      keeper = 1;
      total_pairs = 0;
      pair_index_0 = lat_n2 * dim2 + atm_n1 * dim1 + atm_n2;
      //fprintf(file.out,"input %3d %3d %3d   %5d\n",atm_n1,atm_n2,lat_n2,pair_index_0); fflush(file.out);

      for (pm = 0; pm < number_of_permutations; pm++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {

          cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
          cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

          O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
          O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
          latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + \
          O1_temp];
          latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + \
          O2_temp];

     switch (pm) {

        case 0:
         rotated_cell1 = cell1_temp;
         rotated_cell2 = cell2_temp;
         rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
         rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
         break;

        case 1:
         rotated_cell1 = cell2_temp;
         rotated_cell2 = cell1_temp;
         rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
         rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
         break;

        } // close switch

         pair_index_1 = rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2;
         //pair_index_1 = rotated_cell1 * dim3 + rotated_cell2 * dim4 + rotated_latt2;
         //fprintf(file.out,"%3d %3d output %3d %3d %3d   %5d          %3d %3d\n",\
         latt1_temp,latt2_temp,rotated_cell1,rotated_cell2,rotated_latt2,pair_index_1,pm,k); fflush(file.out);

         if (pair_index_1 <  pair_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }
         //if (pair_index_1 <  pair_index_0) { keeper = 0; pm = 2; k = symmetry->number_of_operators; }

      } // close loop over k
     } // close loop over pm

         if (!keeper) continue;

      for (pm = 0; pm < number_of_permutations; pm++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {

          cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
          cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

          O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
          O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

          latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + \
          O1_temp];
          latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + \
          O2_temp];

     switch (pm) {

        case 0:
         rotated_cell1 = cell1_temp;
         rotated_cell2 = cell2_temp;
         rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
         rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
         break;

        case 1:
         rotated_cell1 = cell2_temp;
         rotated_cell2 = cell1_temp;
         rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
         rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
         break;

        } // close switch

        //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d %3d   %3d %3d\n",pm,k,pair_index_0,pair_index_1,\
        rotated_cell1,rotated_cell2,rotated_latt1,rotated_latt2,pair_p->tot,total_pairs);

   for (j = 0; j <= total_pairs; j++) {
     if (rotated_cell1 == pairs.cell1[j] && rotated_cell2 == pairs.cell2[j] && rotated_latt1 == pairs.latt1[j] && \
         rotated_latt2 == pairs.latt2[j])
         break;
     if (j == total_pairs && rotated_latt2 < R_tables->last_vector && (rotated_cell1 != pairs.cell1[j] || \
              rotated_cell2 != pairs.cell2[j] || rotated_latt1 != pairs.latt1[j] || rotated_latt2 != pairs.latt2[j])) {
         pairs.cell1[j] = rotated_cell1;
         pairs.cell2[j] = rotated_cell2;
         pairs.latt1[j] = rotated_latt1;
         pairs.latt2[j] = rotated_latt2;
         total_pairs++;
         //fprintf(file.out,"pair_p %3d %3d %3d %3d %3d   %3d %3d  %3d %3d\n",\
         j,pair_p->nump,total_pairs,pair_p->tot,rotated_cell1,rotated_cell2,rotated_latt1,rotated_latt2); fflush(file.out);
         break;
        }
       } // close loop over j
      } // close loop over k
     } // close loop over pm
    pair_p->tot += total_pairs;
   (pair_p->nump)++;

   } // close if Rsqrd
   } // close loop over lat_n2
   } // close loop over atm_n2
   } // close loop over atm_n1

  //fprintf(file.out,"unique_pairs %3d total pairs %3d\n",pair_p->nump,pair_p->tot);
  //printf("unique_pairs %3d total pairs %3d\n",pair_p->nump,pair_p->tot);
       
}

void generate_density_pairs4(PAIR_TRAN *pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm, number_of_permutations, p;
int atm_n1, atm_n2, lat_n1, lat_n2;
int cell1_temp, cell2_temp, O1_temp, O2_temp;
int latt1_temp, latt2_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int count4;
int keeper, offset, Offset;
int pair_index_0, pair_index_1;
int unique_pairs, total_pairs;

  double Rsqrd;
  VECTOR_DOUBLE Rvec_tmp;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

  for (i = 0; i < pair_p->nump; i++)
      pair_p->numb[i] = 0;
  //for (i = 0; i < pair_p->nump; i++)
      //pair_p->numi[i] = 0;
  for (i = 0; i < pair_p->nump; i++)
      pair_p->Off[i] = 0;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->cell1[i] = -1;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->latt1[i] = -1;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->cell2[i] = -1;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->off[i] = 0;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->latt2[i] = -1;
  for (i = 0; i < dim2 * R_tables->last_vector; i++)
      pair_p->ptr[i] = -1;
  for (i = 0; i < dim2 * R_tables->last_vector; i++)
      pair_p->Ptr[i] = -1;
  for (i = 0; i < dim2 * R_tables->last_vector; i++)
      pair_p->uniq[i] = 0;
  for (i = 0; i < pair_p->nump * number_of_permutations * symmetry->number_of_operators; i++)
      pair_p->rot[i] = -1;

    pair_p->tot = 0;
    unique_pairs = 0;
    total_pairs = 0;
    offset = 0;
    Offset = 0;
    count4 = 0;

  for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
    for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
      for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
        Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
        Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
        Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
        Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);

        if (Rsqrd < R->cutoff * R->cutoff) {

        count4++;
        lat_n1 = 0;
        keeper = 1;
        pair_index_0 = lat_n2 * dim2 + atm_n1 * dim1 + atm_n2;
        //pair_index_0 = atm_n1 * dim3 + atm_n2 * dim4 + lat_n2;

        for (pm = 0; pm < number_of_permutations; pm++) {
          for (k = 0; k < symmetry->number_of_operators; k++) {

            cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
            cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

            O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
            O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

            latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
            latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];

       switch (pm) {

          case 0:
           rotated_cell1 = cell1_temp;
           rotated_cell2 = cell2_temp;
           rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
           rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
           break;

          case 1:
           rotated_cell1 = cell2_temp;
           rotated_cell2 = cell1_temp;
           rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
           rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
           break;

          } // close switch

           pair_index_1 = rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2;
           //pair_index_1 = rotated_cell1 * dim3 + rotated_cell2 * dim4 + rotated_latt2;

           //fprintf(file.out,"%3d %3d cell %3d %3d latt %3d %3d index %5d %5d\n",\
           O1_temp,O2_temp,cell1_temp,cell2_temp,latt1_temp,latt2_temp,pair_index_0,pair_index_1);

           if (pair_index_1 < pair_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }
           //if (pair_index_1 < pair_index_0) { keeper = 0; pm = 2; k = symmetry->number_of_operators; }

        } // close loop over k
       } // close loop over pm

           if (!keeper) continue;

        for (pm = 0; pm < number_of_permutations; pm++) {
          for (k = 0; k < symmetry->number_of_operators; k++) {

            cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
            cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

            O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
            O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

            latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
            latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];

       switch (pm) {

          case 0:
           rotated_cell1 = cell1_temp;
           rotated_cell2 = cell2_temp;
           rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
           rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
           break;

          case 1:
           rotated_cell1 = cell2_temp;
           rotated_cell2 = cell1_temp;
           rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
           rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
           break;

          } // close switch

          //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d %3d   %3d %3d\n",pm,k,pair_index_0,pair_index_1,\
          rotated_cell1,rotated_cell2,rotated_latt1,rotated_latt2,pair_p->tot,total_pairs);
  
      for (j = pair_p->tot; j <= total_pairs; j++) {
        if (rotated_cell1 == pair_p->cell1[j] && rotated_cell2 == pair_p->cell2[j] && rotated_latt1 == pair_p->latt1[j] && \
            rotated_latt2 == pair_p->latt2[j])
            break;
        if (j == total_pairs && \
           (rotated_cell1 != pair_p->cell1[j] || rotated_cell2 != pair_p->cell2[j] || rotated_latt1 != pair_p->latt1[j] || \
            rotated_latt2 != pair_p->latt2[j])) {
            pair_p->cell1[j] = rotated_cell1;
            pair_p->cell2[j] = rotated_cell2;
            pair_p->latt1[j] = rotated_latt1;
            pair_p->latt2[j] = rotated_latt2;
            //pair_p->k[j] = symmetry->inverse[k];
            pair_p->k[j] = k;
            pair_p->p[j] = pm;
            offset += atoms->bfnnumb_sh[rotated_cell1] * atoms->bfnnumb_sh[rotated_cell2];
            pair_p->off[j + 1] = offset;
            pair_p->ptr[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2] = j;
            // if it is unique, set it to -1
            if (k == 0 && pm == 0) pair_p->uniq[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2] = -1 ; 
            if (k == 0 && pm == 0) pair_p->Ptr[ dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2] = unique_pairs ;
           (pair_p->numb[unique_pairs])++;
            total_pairs++;
            //fprintf(file.out,"pair_p %3d %3d %3d %3d pm %3d op %2d numb %3d   %3d %3d  %3d %3d   %3d %3d     %3d  %3d  off %3d\n",\
            j,unique_pairs,total_pairs,pair_p->tot,pair_p->p[j],k,pair_p->numb[unique_pairs],rotated_cell1,rotated_cell2,\
            pair_p->cell1[j],pair_p->cell2[j],pair_p->latt1[j],pair_p->latt2[j],\
            pair_p->ptr[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2],\
            pair_p->uniq[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2],pair_p->off[j]); 
            break;
           }
          } // close loop over j
         } // close loop over k
        } // close loop over pm
       pair_p->posn[unique_pairs] = pair_p->tot;
       Offset += atoms->bfnnumb_sh[rotated_cell1] * atoms->bfnnumb_sh[rotated_cell2];
       pair_p->Off[unique_pairs + 1]  =  Offset;
       //fprintf(file.out,"%3d %3d  %3d   %5d %4d\n",atm_n1,atm_n2,lat_n2,pair_p->Off[unique_pairs],pair_p->posn[unique_pairs]);
       pair_p->tot = total_pairs;
       unique_pairs++;

      } // close if Rsqrd
     } // close loop over lat_n2
    } // close loop over atm_n2
   } // close loop over atm_n1

      pair_p->nump = unique_pairs;
       
       for (i = 0; i < pair_p->nump; i++) {

           p = pair_p->posn[i];

           for (pm = 0; pm < number_of_permutations; pm++) { 
             for (k = 0; k < symmetry->number_of_operators; k++) {

               atm_n1 = pair_p->cell1[p];
               atm_n2 = pair_p->cell2[p];
               lat_n1 = 0;
               lat_n2 = pair_p->latt2[p];

               cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
               cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

               O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
               O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

               latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
               latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

      //for (j = pair_p_tot1; j <= total_pairs1; j++) {
      for (j = p; j < p + pair_p->numb[i]; j++) {
        if (rotated_cell1 == pair_p->cell1[j] && rotated_cell2 == pair_p->cell2[j] && rotated_latt1 == pair_p->latt1[j] && rotated_latt2 == pair_p->latt2[j]) {
        pair_p->rot[i * number_of_permutations * symmetry->number_of_operators + pm * symmetry->number_of_operators + k] = j;
        //pair_p->rot[i * 2 * symmetry->number_of_operators + pm * symmetry->number_of_operators + k] = j;
        //fprintf(file.out,"pair_p transformation %3d %3d  pm %3d op %2d   %3d %3d  %5d   %3d %3d  %5d  %3d\n",\
        i,j,pm,k,atm_n1,atm_n1,lat_n2,rotated_cell1,rotated_cell2,rotated_latt2,\
        pair_p->rot[i * 2 * symmetry->number_of_operators + pm * symmetry->number_of_operators + k]);
        break;
       }
       } // close loop over j
      } // close loop over k
     } // close loop over pm
    } // close loop over i

     //fprintf(file.out,"unique_pairs %5d total_pairs %5d input pairs %5d\n",pair_p->nump,pair_p->tot,count4);
     //for(j=0;j<dim2 * R_tables->last_vector;j++) fprintf(file.out,"pairs %3d  ptr %3d Ptr %3d uniq %3d\n",\
     j,pair_p->ptr[j],pair_p->Ptr[j],pair_p->uniq[j]);

}

void generate_range_selected_pairs(PAIR_TRAN *pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm, number_of_permutations, p, q;
int atm_n1, atm_n2, lat_n1, lat_n2;
int cell1_temp, cell2_temp, O1_temp, O2_temp;
int latt1_temp, latt2_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int count4;
int keeper, offset, Offset;
int pair_index_0, pair_index_1;
int unique_pairs, total_pairs;
double Rsqrd;
VECTOR_DOUBLE Rvec_tmp;

if (pair_p->cutoff > R->cutoff) {
if (job->taskid == 0) fprintf(file.out,"Pairs with cutoff %7.1f greater than R->cutoff %7.1f requested. Increase R->cutoff\n",\
pair_p->cutoff, R->cutoff);
MPI_Finalize();
exit(0);
}

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

  for (i = 0; i < pair_p->nump; i++)
      pair_p->numb[i] = 0;
  //for (i = 0; i < pair_p->nump; i++)
      //pair_p->numi[i] = 0;
  for (i = 0; i < pair_p->nump; i++)
      pair_p->Off[i] = 0;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->cell1[i] = -1;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->latt1[i] = -1;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->cell2[i] = -1;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->off[i] = 0;
  for (i = 0; i < pair_p->tot; i++)
      pair_p->latt2[i] = -1;
  for (i = 0; i < dim2 * R_tables->last_vector; i++)
      pair_p->ptr[i] = -1;
  for (i = 0; i < dim2 * R_tables->last_vector; i++)
      pair_p->Ptr[i] = -1;
  for (i = 0; i < dim2 * R_tables->last_vector; i++)
      pair_p->uniq[i] = 0;
  for (i = 0; i < pair_p->nump * number_of_permutations * symmetry->number_of_operators; i++)
      pair_p->rot[i] = -1;

    pair_p->tot = 0;
    unique_pairs = 0;
    total_pairs = 0;
    offset = 0;
    Offset = 0;
    count4 = 0;

   for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
     for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
       for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
         Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
         Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
         Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
         Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);

         if (Rsqrd < pair_p->cutoff * pair_p->cutoff) {

         count4++;
         lat_n1 = 0;
         keeper = 1;
         pair_index_0 = lat_n2 * dim2 + atm_n1 * dim1 + atm_n2;
         //pair_index_0 = atm_n1 * dim3 + atm_n2 * dim4 + lat_n2;

         for (pm = 0; pm < number_of_permutations; pm++) {
           for (k = 0; k < symmetry->number_of_operators; k++) {

             cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
             cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

             O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
             O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

             latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
             latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

            pair_index_1 = rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2;
            //pair_index_1 = rotated_cell1 * dim3 + rotated_cell2 * dim4 + rotated_latt2;

            //fprintf(file.out,"%3d %3d cell %3d %3d latt %3d %3d index %5d %5d\n",\
            O1_temp,O2_temp,cell1_temp,cell2_temp,latt1_temp,latt2_temp,pair_index_0,pair_index_1);

            if (pair_index_1 < pair_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }
            //if (pair_index_1 < pair_index_0) { keeper = 0; pm = 2; k = symmetry->number_of_operators; }

         } // close loop over k
        } // close loop over pm

            if (!keeper) continue;

         for (pm = 0; pm < number_of_permutations; pm++) {
           for (k = 0; k < symmetry->number_of_operators; k++) {

             cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
             cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

             O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
             O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

             latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
             latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d %3d   %3d %3d\n",pm,k,pair_index_0,pair_index_1,\
           rotated_cell1,rotated_cell2,rotated_latt1,rotated_latt2,pair_p->tot,total_pairs);
  
      for (j = pair_p->tot; j <= total_pairs; j++) {
        if (rotated_cell1 == pair_p->cell1[j] && rotated_cell2 == pair_p->cell2[j] && rotated_latt1 == pair_p->latt1[j] \
         && rotated_latt2 == pair_p->latt2[j])
            break;
        if (j == total_pairs && \
           (rotated_cell1 != pair_p->cell1[j] || rotated_cell2 != pair_p->cell2[j] || rotated_latt1 != pair_p->latt1[j] || \
            rotated_latt2 != pair_p->latt2[j])) {
            pair_p->cell1[j] = rotated_cell1;
            pair_p->cell2[j] = rotated_cell2;
            pair_p->latt1[j] = rotated_latt1;
            pair_p->latt2[j] = rotated_latt2;
            //pair_p->k[j] = symmetry->inverse[k];
            pair_p->k[j] = k;
            pair_p->p[j] = pm;
            offset += atoms->bfnnumb_sh[rotated_cell1] * atoms->bfnnumb_sh[rotated_cell2];
            pair_p->off[j + 1] = offset;
            pair_p->ptr[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2] = j;
            // if it is unique, set it to -1
            if (k == 0 && pm == 0) pair_p->uniq[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2] = -1 ; 
            if (k == 0 && pm == 0) pair_p->Ptr[ dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2] = unique_pairs ;
           (pair_p->numb[unique_pairs])++;
            total_pairs++;
            //fprintf(file.out,"pair_p %3d %3d %3d %3d pm %3d op %2d numb %3d   %3d %3d  %3d %3d   %3d %3d     %3d  %3d  off %3d\n",\
            j,unique_pairs,total_pairs,pair_p->tot,pair_p->p[j],k,pair_p->numb[unique_pairs],rotated_cell1,rotated_cell2,\
            pair_p->cell1[j],pair_p->cell2[j],pair_p->latt1[j],pair_p->latt2[j],\
            pair_p->ptr[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2],\
            pair_p->uniq[dim2 * rotated_latt2 + dim1 * rotated_cell1 + rotated_cell2],pair_p->off[j]); 
            break;
           }
          } // close loop over j
         } // close loop over k
        } // close loop over pm
       pair_p->posn[unique_pairs] = pair_p->tot;
       Offset += atoms->bfnnumb_sh[rotated_cell1] * atoms->bfnnumb_sh[rotated_cell2];
       pair_p->Off[unique_pairs + 1]  =  Offset;
       //fprintf(file.out,"%3d %3d  %3d   %5d %4d\n",atm_n1,atm_n2,lat_n2,pair_p->Off[unique_pairs],pair_p->posn[unique_pairs]);
       pair_p->tot = total_pairs;
       unique_pairs++;

      } // close if Rsqrd
     } // close loop over lat_n2
    } // close loop over atm_n2
   } // close loop over atm_n1

      pair_p->nump = unique_pairs;
       
       for (i = 0; i < pair_p->nump; i++) {

           p = pair_p->posn[i];

           for (pm = 0; pm < number_of_permutations; pm++) { 
             for (k = 0; k < symmetry->number_of_operators; k++) {

               atm_n1 = pair_p->cell1[p];
               atm_n2 = pair_p->cell2[p];
               lat_n1 = 0;
               lat_n2 = pair_p->latt2[p];

               cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
               cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];

               O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
               O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];

               latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
               latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

      //for (j = pair_p_tot1; j <= total_pairs1; j++) {
      for (j = p; j < p + pair_p->numb[i]; j++) {
        if (rotated_cell1 == pair_p->cell1[j] && rotated_cell2 == pair_p->cell2[j] && rotated_latt1 == pair_p->latt1[j] \
        && rotated_latt2 == pair_p->latt2[j]) {
        pair_p->rot[i * number_of_permutations * symmetry->number_of_operators + pm * symmetry->number_of_operators + k] = j;
        //pair_p->rot[i * 2 * symmetry->number_of_operators + pm * symmetry->number_of_operators + k] = j;
        //fprintf(file.out,"pair_p transformation %3d %3d  pm %3d op %2d   %3d %3d  %5d   %3d %3d  %5d  %3d\n",\
        i,j,pm,k,atm_n1,atm_n1,lat_n2,rotated_cell1,rotated_cell2,rotated_latt2,\
        pair_p->rot[i * 2 * symmetry->number_of_operators + pm * symmetry->number_of_operators + k]);
        break;
       }
       } // close loop over j
      } // close loop over k
     } // close loop over pm
    } // close loop over i

    for (i = 0; i < pair_p->nump; i++) {
      q  = pair_p->posn[i];
      for (j = 0; j < pair_p->numb[i]; j++) {
        atm_n1 = pair_p->cell1[q + j];
        atm_n2 = pair_p->cell2[q + j];
        lat_n1 = pair_p->latt1[q + j];
        lat_n2 = pair_p->latt2[q + j];
        Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 + R->vec_ai[lat_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - \
        R->vec_ai[lat_n2].comp1;
        Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 + R->vec_ai[lat_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - \
        R->vec_ai[lat_n2].comp2;
        Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 + R->vec_ai[lat_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - \
        R->vec_ai[lat_n2].comp3;
        pair_p->dist[q + j] = sqrt(double_vec_dot(&Rvec_tmp,&Rvec_tmp));
 //fprintf(file.out,"pair[%5d]   %3d %3d  %3d   transforms to pair[%5d]  %3d %3d  %5d under operator %2d permutation %2d %10.4lf\n",\
        q, pair_p->cell1[q], pair_p->cell2[q], pair_p->latt2[q], q + j, pair_p->cell1[q + j], pair_p->cell2[q + j], \
        pair_p->latt2[q + j],pair_p->k[q + j],pair_p->p[q + j],pair_p->dist * bohr_to_AA);
       }
       //fprintf(file.out,"\n");
      }

     //fprintf(file.out,"unique_pairs %5d total_pairs %5d input pairs %5d\n",pair_p->nump,pair_p->tot,count4);
     //for(j=0;j<dim2 * R_tables->last_vector;j++) fprintf(file.out,"pairs %3d  ptr %3d Ptr %3d uniq %3d\n",\
     j,pair_p->ptr[j],pair_p->Ptr[j],pair_p->uniq[j]);

}

void generate_triples(TRIPLE_TRAN *Triple, int atm_n1, int atm_n2, int atm_n3, int lat_n1, int lat_n2, int lat_n3, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)
//void generate_triples(TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
int dil1 = R_tables->last_vector;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dil2 = dim1 * dil1;
int dil3 = dil1 * dil2;
int dil4 = dil3 * dim1;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 6;

     Triple->tot = 0;

     for (i=0;i<8 * symmetry->number_of_operators;i++) Triple->cell1[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Triple->cell2[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Triple->cell3[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Triple->latt1[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Triple->latt2[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Triple->latt3[i] = -1;

     triple_index_0 = atm_n1 * dil4 + atm_n2 * dil3 + lat_n2 * dil2 + atm_n3 * dil1 + lat_n3;

     //fprintf(file.out,"INPUT  %5d %3d %3d %3d  %3d %3d %3d\n",triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);

    for (pm = 0; pm < number_of_permutations; pm++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {

      keeper = 0;

      cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
      cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
      cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

      O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
      O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
      O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

      latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
      latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
      latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];

     //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell3_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            break;

           } // close switch

            triple_index_1 = rotated_cell1 * dil4 + rotated_cell2 * dil3 + rotated_latt2 * dil2 + rotated_cell3 * dil1 + rotated_latt3;

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
           pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
           rotated_latt2,rotated_latt3);



           if (triple_index_1 >= triple_index_0) keeper = 1;
           if (triple_index_1 < triple_index_0) { Triple->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }

            //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d  latt  %3d %3d %3d\n",i,pm,k,rotated_cell1,rotated_cell2,\
            rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

            if (!keeper) continue;

            for (j = 0; j <= Triple->tot; j++) {
              if (rotated_cell1 == Triple->cell1[j] && rotated_latt1 == Triple->latt1[j] && \
                  rotated_cell2 == Triple->cell2[j] && rotated_latt2 == Triple->latt2[j] && \
                  rotated_cell3 == Triple->cell3[j] && rotated_latt3 == Triple->latt3[j])
                  break;
            if (j == Triple->tot && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                 (rotated_cell1 != Triple->cell1[j] || rotated_latt1 != Triple->latt1[j] || \
                  rotated_cell2 != Triple->cell2[j] || rotated_latt2 != Triple->latt2[j] || \
                  rotated_cell3 != Triple->cell3[j] || rotated_latt3 != Triple->latt3[j])) {
                  Triple->cell1[j] = rotated_cell1;
                  Triple->cell2[j] = rotated_cell2;
                  Triple->cell3[j] = rotated_cell3;
                  Triple->latt1[j] = rotated_latt1;
                  Triple->latt2[j] = rotated_latt2;
                  Triple->latt3[j] = rotated_latt3;
                  Triple->k[j] = k;
                  Triple->p[j] = pm;
                  (Triple->tot)++;
                  //fprintf(file.out,"Triples %3d %3d pm %3d op %2d      %3d %3d %3d       %3d %3d %3d     %3d %3d %3d  %3d %3d %3d  %3d\n",\
                  j,Triple->tot,Triple->p[j],k, \
                  rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  rotated_latt1,rotated_latt2,rotated_latt3,Triple->latt1[j],Triple->latt2[j],Triple->latt3[j],keeper); 
                  //printf("Triples %3d %3d pm %3d op %2d      %3d %3d %3d       %3d %3d %3d        %3d\n",j,Triple->tot,Triple->p[j],k, \
                  rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],keeper); 
                  break;
                 }
                } // close loop over j
               } // close loop over k
              } // close loop over pm

              //if (Triple->tot > 0) fprintf(file.out,"\n");
              //for(j=0;j<Triple->tot;j++) fprintf(file.out,"gathered e Triples %3d pm %3d op %3d   %3d %3d %3d   %3d %3d %3d\n",\
              j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
              Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);
              //printf("\n");
              //for(j=0;j<Triple->tot;j++) printf("gathered e Triples %3d pm %3d op %3d   %3d %3d %3d   %3d %3d %3d\n",\
              j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
              Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);

}

void count_triples1(TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, atm_n3, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
int dil1 = R_tables->last_vector;
int dim3 = dil1;
int dil2 = atoms->number_of_atoms_in_unit_cell * dil1;
int dim2 = dil1 * dil2;
int dim1 = dil2 * dil2;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples;
TRIPLE_TMP triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;

    for (i = 0; i < MXT; i++) triples.cell1[i] = -1;
    for (i = 0; i < MXT; i++) triples.cell2[i] = -1;
    for (i = 0; i < MXT; i++) triples.cell3[i] = -1;
    for (i = 0; i < MXT; i++) triples.latt1[i] = -1;
    for (i = 0; i < MXT; i++) triples.latt2[i] = -1;
    for (i = 0; i < MXT; i++) triples.latt3[i] = -1;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell1[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell2[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell3[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt1[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt2[i] = -1;
    for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt3[i] = -1;

     Triple->tot = 0;
     Triple->nump = 0;
     unique_triples = 0;
     total_triples = 0;

     for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
       for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
         for (atm_n3 = 0; atm_n3 < atoms->number_of_atoms_in_unit_cell; atm_n3++) {
           for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
             for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
               Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
               Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
               Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
               Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
               Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
               Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
               Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
               Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
               Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
               Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
               Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
               Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             //fprintf(file.out,"INPUT  %5d %3d %3d %3d  %3d %3d %3d\n",triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
             //fflush(file.out);
             if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
               lat_n1 = 0;
               keeper = 1;
               triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;

               for (pm = 0; pm < number_of_permutations; pm++) {
                 for (k = 0; k < symmetry->number_of_operators; k++) {

                   cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                   cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                   cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                   O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                   O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                   O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                   latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
                   latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
                   latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];
                   //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                   //fflush(file.out);

                   switch (pm) {

                     case 0:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 1:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 2:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 3:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 4:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     case 5:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     } // close switch

                      triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + rotated_cell3 * dim3 + rotated_latt3;
                      //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                      pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                      rotated_latt2,rotated_latt3);

                      if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                      //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d  latt  %3d %3d %3d\n",i,pm,k,rotated_cell1,rotated_cell2,\
                      rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                         } // close loop over k
                        } // close loop over pm

                      if (!keeper) continue;

               for (pm = 0; pm < number_of_permutations; pm++) {
                 for (k = 0; k < symmetry->number_of_operators; k++) {

                   cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                   cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                   cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                   O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                   O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                   O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                   latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
                   latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
                   latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];
                   //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                   //fflush(file.out);

                   switch (pm) {

                     case 0:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 1:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 2:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 3:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 4:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     case 5:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     } // close switch

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d   %3d %3d %3d   %3d %3d\n",pm,k,triple_index_0,triple_index_1,\
           rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);
  
                      for (j = 0; j <= total_triples; j++) {
                        if (rotated_cell1 == triples.cell1[j] && rotated_latt1 == triples.latt1[j] && \
                            rotated_cell2 == triples.cell2[j] && rotated_latt2 == triples.latt2[j] && \
                            rotated_cell3 == triples.cell3[j] && rotated_latt3 == triples.latt3[j])
                            break;
                      if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                           (rotated_cell1 != triples.cell1[j] || rotated_latt1 != triples.latt1[j] || \
                            rotated_cell2 != triples.cell2[j] || rotated_latt2 != triples.latt2[j] || \
                            rotated_cell3 != triples.cell3[j] || rotated_latt3 != triples.latt3[j])) {
                            triples.cell1[j] = rotated_cell1;
                            triples.cell2[j] = rotated_cell2;
                            triples.cell3[j] = rotated_cell3;
                            triples.latt1[j] = rotated_latt1;
                            triples.latt2[j] = rotated_latt2;
                            triples.latt3[j] = rotated_latt3;
                            total_triples++;
                            //fprintf(file.out,"triples %3d %3d  %3d pm %3d op %2d      %3d %3d %3d       %3d %3d %3d     %3d %3d %3d  %3d %3d %3d  %3d\n",\
                            j,unique_triples+1,total_triples, pm, k, \
                            rotated_cell1,rotated_cell2,rotated_cell3,triples.cell1[j],triples.cell2[j],triples.cell3[j],\
                            rotated_latt1,rotated_latt2,rotated_latt3,triples.latt1[j],triples.latt2[j],triples.latt3[j],keeper); 
                            break;
                           }
                          } // close loop over j
                         } // close loop over k
                        } // close loop over pm
                         unique_triples++;
                         Triple->tot = total_triples;
                        } // close if (Rsqrd > 
                       } // close loop over lat_n3
                      } // close loop over lat_n2
                     } // close loop over atm_n3
                    } // close loop over atm_n2
                   } // close loop over atm_n1

                     Triple->nump = unique_triples;

}

void generate_triples1(TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, atm_n3, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
int dil1 = R_tables->last_vector;
int dim3 = dil1;
int dil2 = atoms->number_of_atoms_in_unit_cell * dil1;
int dim2 = dil1 * dil2;
int dim1 = dil2 * dil2;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;


    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

  for (i = 0; i < Triple->nump; i++)
      Triple->numb[i] = 0;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell3[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt3[i] = -1;

     Triple->tot = 0;
     unique_triples = 0;
     total_triples = 0;

     for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell1[i] = -1;
     for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell2[i] = -1;
     for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell3[i] = -1;
     for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt1[i] = -1;
     for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt2[i] = -1;
     for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt3[i] = -1;

     for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
       for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
         for (atm_n3 = 0; atm_n3 < atoms->number_of_atoms_in_unit_cell; atm_n3++) {
           for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
             for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
               Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
               Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
               Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
               Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
               Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
               Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
               Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
               Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
               Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
               Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
               Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
               Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
               lat_n1 = 0;
               keeper = 1;
               triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;
               //fprintf(file.out,"INPUT %10.4lf %10.4lf %10.4lf %10.4lf %5d %3d %3d %3d  %3d %3d %3d\n", \
               Rsqrd12,Rsqrd13,Rsqrd23,R->cutoff * R->cutoff,triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
               //fflush(file.out);

               for (pm = 0; pm < number_of_permutations; pm++) {
                 for (k = 0; k < symmetry->number_of_operators; k++) {

                   cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                   cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                   cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                   O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                   O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                   O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                   latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
                   latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
                   latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];
                   //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                   //fflush(file.out);

                   switch (pm) {

                     case 0:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 1:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 2:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 3:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 4:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     case 5:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     } // close switch

                      triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + rotated_cell3 * dim3 + rotated_latt3;
                      //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                      pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                      rotated_latt2,rotated_latt3);

                      if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                      //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d  latt  %3d %3d %3d\n",i,pm,k,rotated_cell1,rotated_cell2,\
                      rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                         } // close loop over k
                        } // close loop over pm

                      if (!keeper) continue;

               for (pm = 0; pm < number_of_permutations; pm++) {
                 for (k = 0; k < symmetry->number_of_operators; k++) {

                   cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                   cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                   cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                   O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                   O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                   O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                   latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
                   latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
                   latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];
                   //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                   //fflush(file.out);

                   switch (pm) {

                     case 0:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 1:
                      rotated_cell1 = cell1_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                      break;

                     case 2:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell3_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 3:
                      rotated_cell1 = cell2_temp;
                      rotated_cell2 = cell3_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                      break;

                     case 4:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell1_temp;
                      rotated_cell3 = cell2_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     case 5:
                      rotated_cell1 = cell3_temp;
                      rotated_cell2 = cell2_temp;
                      rotated_cell3 = cell1_temp;
                      rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                      rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                      break;

                     } // close switch

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d   %3d %3d %3d   %3d %3d\n",pm,k,triple_index_0,triple_index_1,\
           rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);
  
                      for (j = 0; j <= total_triples; j++) {
                        if (rotated_cell1 == Triple->cell1[j] && rotated_latt1 == Triple->latt1[j] && \
                            rotated_cell2 == Triple->cell2[j] && rotated_latt2 == Triple->latt2[j] && \
                            rotated_cell3 == Triple->cell3[j] && rotated_latt3 == Triple->latt3[j])
                            break;
                      if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                           (rotated_cell1 != Triple->cell1[j] || rotated_latt1 != Triple->latt1[j] || \
                            rotated_cell2 != Triple->cell2[j] || rotated_latt2 != Triple->latt2[j] || \
                            rotated_cell3 != Triple->cell3[j] || rotated_latt3 != Triple->latt3[j])) {
                            Triple->cell1[j] = rotated_cell1;
                            Triple->cell2[j] = rotated_cell2;
                            Triple->cell3[j] = rotated_cell3;
                            Triple->latt1[j] = rotated_latt1;
                            Triple->latt2[j] = rotated_latt2;
                            Triple->latt3[j] = rotated_latt3;
                            Triple->k[j] = k;
                            Triple->p[j] = pm;
                           (Triple->numb[unique_triples])++;
                            total_triples++;
                            //fprintf(file.out,"Triples %3d %3d %3d %3d pm %3d op %2d      %3d %3d %3d     %3d %3d %3d     %3d %3d %3d  %3d %3d %3d  %3d\n",\
                            j,unique_triples+1,total_triples,Triple->numb[unique_triples],pm,k, \
                            rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                            rotated_latt1,rotated_latt2,rotated_latt3,Triple->latt1[j],Triple->latt2[j],Triple->latt3[j],keeper); 
                            //printf("Triples %3d %3d pm %3d op %2d      %3d %3d %3d       %3d %3d %3d        %3d\n",j,Triple->tot,Triple->p[j],k, \
                            rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],keeper); 
                            break;
                           }
                          } // close loop over j
                         } // close loop over k
                        } // close loop over pm
                         Triple->posn[unique_triples] = Triple->tot;
                         unique_triples++;
                         Triple->tot = total_triples;
                        } // close if (Rsqrd >
                       } // close loop over lat_n3
                      } // close loop over lat_n2
                     } // close loop over atm_n3
                    } // close loop over atm_n2
                   } // close loop over atm_n1

                     Triple->nump = unique_triples;

                     //if (Triple->tot > 0) fprintf(file.out,"\n");
                     //for(j=0;j<Triple->tot;j++) fprintf(file.out,"gathered e Triples %3d pm %3d op %3d   %3d %3d %3d   %3d %3d %3d\n",\
                     j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                     Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);
                     //printf("\n");
                     //for(j=0;j<Triple->tot;j++) printf("gathered e Triples %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d\n",\
                     j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                     Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);

}

void count_triples1_reversed(int atm_n3, TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, lat_n1, lat_n2, lat_n3;
//int atm_n1, atm_n2, atm_n3, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
int dim2 = R_tables->last_vector;
int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;
//int dil3 = atoms->number_of_atoms_in_unit_cell * dim1;
//int dim3 = dim2 * dil3;
//int dil1 = R_tables->last_vector;
//int dim3 = dil1;
//int dil2 = atoms->number_of_atoms_in_unit_cell * dil1;
//int dim2 = dil1 * dil2;
//int dim1 = dil2 * dil2;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples;
TRIPLE_TMP triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;

  for (i = 0; i < MXT; i++) triples.cell1[i] = -1;
  for (i = 0; i < MXT; i++) triples.cell2[i] = -1;
  for (i = 0; i < MXT; i++) triples.cell3[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt1[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt2[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt3[i] = -1;

  if (job->pms == 0) number_of_permutations = 1;
  else               number_of_permutations = 2;

  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell1[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell2[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell3[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt1[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt2[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt3[i] = -1;

   Triple->tot = 0;
   Triple->nump = 0;
   unique_triples = 0;
   total_triples = 0;

   lat_n3 = 0;

   for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
     for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
       //for (atm_n3 = 0; atm_n3 < atoms->number_of_atoms_in_unit_cell; atm_n3++) {
         for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
           //for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
             Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
             Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n3].comp1 - \
             R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n3].comp2 - \
             R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n3].comp3 - \
             R->vec_ai[lat_n3].comp3;
             Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
           //fprintf(file.out,"INPUT  %5d %3d %3d %3d  %3d %3d %3d\n",triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
           //fflush(file.out);
           if (Rsqrd12 < R->cutoff * R->cutoff) { // Rsqrd13 and Rsqrd23 are Coulombic separations and are not restricted
           //if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
           //if (Rsqrd12 < R->cutoff * R->cutoff) {
             lat_n1 = 0;
             keeper = 1;

             // modified for new order for atm_n1, atm_n2 and atm_n3
             triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
             //triple_index_0 = atm_n3 * dim3 + lat_n3 * dil3 + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
             //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;

             for (pm = 0; pm < number_of_permutations; pm++) {
               for (k = 0; k < symmetry->number_of_operators; k++) {

                 cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                 cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                 cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                 if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                 O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                 O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                 O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                 latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O1_temp];
                 latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O2_temp];
                 latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O3_temp];

                 //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                 latt3_temp);
                 //fflush(file.out);

                 switch (pm) {

                   case 0:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    break;
                    // for reversed routine swap cases 1 and 2 for permutations
                   case 1:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 2:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_cell1 = cell2_temp;
                    //rotated_cell2 = cell1_temp;
                    //rotated_cell3 = cell3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 3:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 4:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   case 5:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   } // close switch

                    // modified for new order for atm_n1, atm_n2 and atm_n3
                    triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2;
                    //triple_index_1 = rotated_cell3 * dim3 + rotated_latt3 * dil3 + rotated_cell1 * dim1 + rotated_cell2 * dim2 + \
                    rotated_latt2;
                    //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + rotated_cell3 * dim3 + \
                    rotated_latt3;
                    //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                    pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                    rotated_latt2,rotated_latt3);

                 if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                    //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell %3d %3d %3d latt  %3d %3d %3d\n",i,pm,k,rotated_cell1,\
                    rotated_cell2,\
                    rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                       } // close loop over k
                      } // close loop over pm

                    if (!keeper) continue;

             for (pm = 0; pm < number_of_permutations; pm++) {
               for (k = 0; k < symmetry->number_of_operators; k++) {

                 cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                 cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                 cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                 if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                 O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                 O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                 O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                 latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O1_temp];
                 latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O2_temp];
                 latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O3_temp];

                 //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                 latt3_temp);
                 //fflush(file.out);

                 switch (pm) {

                   case 0:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    break;

                    // for reversed routine swap cases 1 and 2 for permutations
                   case 1:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 2:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_cell1 = cell2_temp;
                    //rotated_cell2 = cell1_temp;
                    //rotated_cell3 = cell3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 3:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 4:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   case 5:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   } // close switch

         //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d    %3d %3d %3d   %3d %3d %3d   %3d %3d\n",\
         pm,k,triple_index_0,triple_index_1,\
         rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);

                    for (j = 0; j <= total_triples; j++) {
                      if (rotated_cell1 == triples.cell1[j] && rotated_latt1 == triples.latt1[j] && \
                          rotated_cell2 == triples.cell2[j] && rotated_latt2 == triples.latt2[j] && \
                          rotated_cell3 == triples.cell3[j] && rotated_latt3 == triples.latt3[j])
                          break;
                    if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                         (rotated_cell1 != triples.cell1[j] || rotated_latt1 != triples.latt1[j] || \
                          rotated_cell2 != triples.cell2[j] || rotated_latt2 != triples.latt2[j] || \
                          rotated_cell3 != triples.cell3[j] || rotated_latt3 != triples.latt3[j])) {
                          triples.cell1[j] = rotated_cell1;
                          triples.cell2[j] = rotated_cell2;
                          triples.cell3[j] = rotated_cell3;
                          triples.latt1[j] = rotated_latt1;
                          triples.latt2[j] = rotated_latt2;
                          triples.latt3[j] = rotated_latt3;
                          total_triples++;
                //fprintf(file.out,"triples %3d %3d  %3d pm %3d op %2d  %3d %3d %3d   %3d %3d %3d  %3d %3d %3d  %3d %3d %3d  %3d\n",\
                          j,unique_triples+1,total_triples, pm, k, \
                          rotated_cell1,rotated_cell2,rotated_cell3,triples.cell1[j],triples.cell2[j],triples.cell3[j],\
                          rotated_latt1,rotated_latt2,rotated_latt3,triples.latt1[j],triples.latt2[j],triples.latt3[j],keeper); 
                          break;
                         }
                        } // close loop over j
                       } // close loop over k
                      } // close loop over pm
                       unique_triples++;
                       Triple->tot = total_triples;
                      } // close if (Rsqrd > 
                     //} // close loop over lat_n3
                    } // close loop over lat_n2
                   //} // close loop over atm_n3
                  } // close loop over atm_n2
                 } // close loop over atm_n1

                 Triple->nump = unique_triples;

}

void generate_triples1_reversed(int atm_n3, TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, lat_n1, lat_n2, lat_n3;
//int atm_n1, atm_n2, atm_n3, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
int dim2 = R_tables->last_vector;
int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;
//int dil1 = R_tables->last_vector;
//int dim3 = dil1;
//int dil2 = atoms->number_of_atoms_in_unit_cell * dil1;
//int dim2 = dil1 * dil2;
//int dim1 = dil2 * dil2;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples; //
//int pos3, unique_pos3_triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;


    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

  for (i = 0; i < Triple->nump; i++)
      Triple->numb[i] = 0;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell3[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt3[i] = -1;

  Triple->tot = 0;
  unique_triples = 0;
  total_triples = 0;

  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell1[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell2[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell3[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt1[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt2[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt3[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell1[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell2[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell3[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt1[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt2[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt3[i] = -1;

  lat_n3 = 0;

      for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
        for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
        //for (atm_n3 = 0; atm_n3 < atoms->number_of_atoms_in_unit_cell; atm_n3++) {
          for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
          //for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
            Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
            Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
            Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
            Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
            Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
            Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
            Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
            Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
            Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n3].comp1 - \
            R->vec_ai[lat_n3].comp1;
            Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n3].comp2 - \
            R->vec_ai[lat_n3].comp2;
            Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n3].comp3 - \
            R->vec_ai[lat_n3].comp3;
            Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
            if (Rsqrd12 < R->cutoff * R->cutoff) { // Rsqrd13 and Rsqrd23 are Coulombic separations and are not restricted
            //if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
            //if (Rsqrd12 < R->cutoff * R->cutoff) {
            lat_n1 = 0;
            keeper = 1;

            // modified for new order for atm_n1, atm_n2 and atm_n3
            triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
            //triple_index_0 = atm_n3 * dim3 + lat_n3 * dil3 + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;

            //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;
            //fprintf(file.out,"INPUT %10.4lf %10.4lf %10.4lf %10.4lf %5d   %3d %3d %3d  %3d %3d %3d\n", \
            Rsqrd12,Rsqrd13,Rsqrd23,R->cutoff * R->cutoff,triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
            //fflush(file.out);

            for (pm = 0; pm < number_of_permutations; pm++) {
              for (k = 0; k < symmetry->number_of_operators; k++) {

                cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O1_temp];
                latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O2_temp];
                latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O3_temp];

                //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                latt3_temp);
                //fflush(file.out);

                switch (pm) {

                  case 0:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   break;

                   // for reversed routine swap cases 1 and 2 for permutations
                  case 1:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_cell1 = cell1_temp;
                   //rotated_cell2 = cell3_temp;
                   //rotated_cell3 = cell2_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 2:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 3:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   break;

                  case 4:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  case 5:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  } // close switch

                   // modified for new order for atm_n1, atm_n2 and atm_n3
                   triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2;
                   //triple_index_1 = rotated_cell3 * dim3 + rotated_latt3 * dil3 + rotated_cell1 * dim1 + \
                   rotated_cell2 * dim2 + rotated_latt2;
                   //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + \
                   rotated_cell3 * dim3 + rotated_latt3;
                   //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                   pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                   rotated_latt2,rotated_latt3);

                 if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                   //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d  latt  %3d %3d %3d\n",\
                   i,pm,k,rotated_cell1,rotated_cell2,\
                   rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                      } // close loop over k
                     } // close loop over pm

                   if (!keeper) continue;

            for (pm = 0; pm < number_of_permutations; pm++) {
              for (k = 0; k < symmetry->number_of_operators; k++) {

                cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O1_temp];
                latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O2_temp];
                latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O3_temp];
                //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                //fflush(file.out);

                switch (pm) {

                  case 0:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   break;

                   // for reversed routine swap cases 1 and 2 for permutations
                  case 1:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_cell1 = cell1_temp;
                   //rotated_cell2 = cell3_temp;
                   //rotated_cell3 = cell2_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 2:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 3:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   break;

                  case 4:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  case 5:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  } // close switch

        //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d   %3d %3d %3d   %3d %3d\n", \
        pm,k,triple_index_0,triple_index_1,\
        rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);

                   for (j = 0; j <= total_triples; j++) {
                     if (rotated_cell1 == Triple->cell1[j] && rotated_latt1 == Triple->latt1[j] && \
                         rotated_cell2 == Triple->cell2[j] && rotated_latt2 == Triple->latt2[j] && \
                         rotated_cell3 == Triple->cell3[j] && rotated_latt3 == Triple->latt3[j])
                         break;
                   if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                        (rotated_cell1 != Triple->cell1[j] || rotated_latt1 != Triple->latt1[j] || \
                         rotated_cell2 != Triple->cell2[j] || rotated_latt2 != Triple->latt2[j] || \
                         rotated_cell3 != Triple->cell3[j] || rotated_latt3 != Triple->latt3[j])) {
                         Triple->cell1[j] = rotated_cell1;
                         Triple->cell2[j] = rotated_cell2;
                         Triple->cell3[j] = rotated_cell3;
                         Triple->latt1[j] = rotated_latt1;
                         Triple->latt2[j] = rotated_latt2;
                         Triple->latt3[j] = rotated_latt3;
                         Triple->k[j] = k;
                         Triple->p[j] = pm;
                        (Triple->numb[unique_triples])++;
                         total_triples++;
              //fprintf(file.out,"Triples %3d %3d %3d %3d pm %3d op %2d  %3d %3d %3d  %3d %3d %3d %3d %3d %3d  %3d %3d %3d  %3d\n",\
                         j,unique_triples+1,total_triples,Triple->numb[unique_triples],pm,k, \
                         rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                         rotated_latt1,rotated_latt2,rotated_latt3,Triple->latt1[j],Triple->latt2[j],Triple->latt3[j],keeper); 
              //printf("Triples %3d %3d pm %3d op %2d    %3d %3d %3d       %3d %3d %3d        %3d\n",j,Triple->tot,Triple->p[j],k, \
                         rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],keeper); 
                         break;
                        }
                       } // close loop over j
                      } // close loop over k
                     } // close loop over pm
                      Triple->posn[unique_triples] = Triple->tot;
                      unique_triples++;
                      Triple->tot = total_triples;
                     } // close if (Rsqrd >
                    //} // close loop over lat_n3
                   } // close loop over lat_n2
                  //} // close loop over atm_n3
                 } // close loop over atm_n2
                } // close loop over atm_n1

                  Triple->nump = unique_triples;

                  //if (Triple->tot > 0) fprintf(file.out,"\n");
                  //for(j=0;j<Triple->tot;j++) fprintf(file.out,"gathered e Triples %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d\n",\
                  j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);
                  //for(j=0;j<Triple->nump;j++) fprintf(file.out,"%3d posn %3d numb %3d\n",\
                  j,Triple->posn[j],Triple->numb[j]);
                  //printf("\n");
                  //for(j=0;j<Triple->tot;j++) printf("gathered e Triples %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d\n",\
                  j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);

}

void count_triples2_reversed(int atm_n3, TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, lat_n1, lat_n2, lat_n3;
//int atm_n1, atm_n2, atm_n3, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
//int dim2 = R_tables->last_vector;
//int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;

int dil1 = R_tables->last_vector;
int dim2 = dil1 * R_tables->last_vector;
int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;

int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples;
TRIPLE_TMP triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;

  for (i = 0; i < MXT; i++) triples.cell1[i] = -1;
  for (i = 0; i < MXT; i++) triples.cell2[i] = -1;
  for (i = 0; i < MXT; i++) triples.cell3[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt1[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt2[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt3[i] = -1;

  if (job->pms == 0) number_of_permutations = 1;
  else               number_of_permutations = 2;

  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell1[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell2[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell3[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt1[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt2[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt3[i] = -1;

   Triple->tot = 0;
   Triple->nump = 0;
   unique_triples = 0;
   total_triples = 0;

   lat_n3 = 0;

   for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
     for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
       //for (atm_n3 = 0; atm_n3 < atoms->number_of_atoms_in_unit_cell; atm_n3++) {
       for (lat_n1 = 0; lat_n1 < R->max_vector; lat_n1++) {
         for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
           //for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 + R->vec_ai[lat_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - \
             R->vec_ai[lat_n2].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 + R->vec_ai[lat_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - \
             R->vec_ai[lat_n2].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 + R->vec_ai[lat_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - \
             R->vec_ai[lat_n2].comp3;
             Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 + R->vec_ai[lat_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - \
             R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 + R->vec_ai[lat_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - \
             R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 + R->vec_ai[lat_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - \
             R->vec_ai[lat_n3].comp3;
             Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n3].comp1 - \
             R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n3].comp2 - \
             R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n3].comp3 - \
             R->vec_ai[lat_n3].comp3;
             Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
           //fprintf(file.out,"INPUT  %5d %3d %3d %3d  %3d %3d %3d\n",triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
           //fflush(file.out);
           if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
             //lat_n1 = 0;
             keeper = 1;

             // modified for new order for atm_n1, atm_n2 and atm_n3
             //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
             triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n1 * dil1 + lat_n2;
             //triple_index_0 = atm_n3 * dim3 + lat_n3 * dil3 + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
             //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;

             for (pm = 0; pm < number_of_permutations; pm++) {
               for (k = 0; k < symmetry->number_of_operators; k++) {

                 cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                 cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                 cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                 if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                 O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                 O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                 O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                 latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O1_temp];
                 latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O2_temp];
                 latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O3_temp];

                 //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                 latt3_temp);
                 //fflush(file.out);

                 switch (pm) {

                   case 0:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = latt1_temp;
                    rotated_latt2 = latt2_temp;
                    rotated_latt3 = latt3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    break;
                    // for reversed routine swap cases 1 and 2 for permutations
                   case 1:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = latt2_temp;
                    rotated_latt2 = latt1_temp;
                    rotated_latt3 = latt3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 2:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_cell1 = cell2_temp;
                    //rotated_cell2 = cell1_temp;
                    //rotated_cell3 = cell3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 3:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 4:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   case 5:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   } // close switch

                    // modified for new order for atm_n1, atm_n2 and atm_n3
                    //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2;
                    triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt1 * dil1 + rotated_latt2;
                    //triple_index_1 = rotated_cell3 * dim3 + rotated_latt3 * dil3 + rotated_cell1 * dim1 + rotated_cell2 * dim2 + \
                    rotated_latt2;
                    //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + rotated_cell3 * dim3 + \
                    rotated_latt3;
                    //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                    pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                    rotated_latt2,rotated_latt3);

                 if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                    //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell %3d %3d %3d latt  %3d %3d %3d\n",i,pm,k,rotated_cell1,\
                    rotated_cell2,\
                    rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                       } // close loop over k
                      } // close loop over pm

                    if (!keeper) continue;

             for (pm = 0; pm < number_of_permutations; pm++) {
               for (k = 0; k < symmetry->number_of_operators; k++) {

                 cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                 cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                 cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                 if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                 O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                 O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                 O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                 latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O1_temp];
                 latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O2_temp];
                 latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O3_temp];

                 //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                 latt3_temp);
                 //fflush(file.out);

                 switch (pm) {

                   case 0:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = latt1_temp;
                    rotated_latt2 = latt2_temp;
                    rotated_latt3 = latt3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    break;

                    // for reversed routine swap cases 1 and 2 for permutations
                   case 1:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = latt2_temp;
                    rotated_latt2 = latt1_temp;
                    rotated_latt3 = latt3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 2:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_cell1 = cell2_temp;
                    //rotated_cell2 = cell1_temp;
                    //rotated_cell3 = cell3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 3:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 4:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   case 5:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   } // close switch

         //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d    %3d %3d %3d   %3d %3d %3d   %3d %3d\n",\
         pm,k,triple_index_0,triple_index_1,\
         rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);

                    for (j = 0; j <= total_triples; j++) {
                      if (rotated_cell1 == triples.cell1[j] && rotated_latt1 == triples.latt1[j] && \
                          rotated_cell2 == triples.cell2[j] && rotated_latt2 == triples.latt2[j] && \
                          rotated_cell3 == triples.cell3[j] && rotated_latt3 == triples.latt3[j])
                          break;
                    if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                         (rotated_cell1 != triples.cell1[j] || rotated_latt1 != triples.latt1[j] || \
                          rotated_cell2 != triples.cell2[j] || rotated_latt2 != triples.latt2[j] || \
                          rotated_cell3 != triples.cell3[j] || rotated_latt3 != triples.latt3[j])) {
                          triples.cell1[j] = rotated_cell1;
                          triples.cell2[j] = rotated_cell2;
                          triples.cell3[j] = rotated_cell3;
                          triples.latt1[j] = rotated_latt1;
                          triples.latt2[j] = rotated_latt2;
                          triples.latt3[j] = rotated_latt3;
                          total_triples++;
                //fprintf(file.out,"triples %3d %3d  %3d pm %3d op %2d  %3d %3d %3d   %3d %3d %3d  %3d %3d %3d  %3d %3d %3d  %3d\n",\
                          j,unique_triples+1,total_triples, pm, k, \
                          rotated_cell1,rotated_cell2,rotated_cell3,triples.cell1[j],triples.cell2[j],triples.cell3[j],\
                          rotated_latt1,rotated_latt2,rotated_latt3,triples.latt1[j],triples.latt2[j],triples.latt3[j],keeper); 
                          break;
                         }
                        } // close loop over j
                       } // close loop over k
                      } // close loop over pm
                       unique_triples++;
                       Triple->tot = total_triples;
                      } // close if (Rsqrd > 
                     //} // close loop over lat_n3
                    } // close loop over lat_n2
                   } // close loop over lat_n1
                   //} // close loop over atm_n3
                  } // close loop over atm_n2
                 } // close loop over atm_n1

                 Triple->nump = unique_triples;

}

void generate_triples2_reversed(int atm_n3, TRIPLE_TRAN *Triple, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, lat_n1, lat_n2, lat_n3;
//int atm_n1, atm_n2, atm_n3, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3;
//int dim2 = R_tables->last_vector;
//int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;

int dil1 = R_tables->last_vector;
int dim2 = dil1 * R_tables->last_vector;
int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;

int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples; //
//int pos3, unique_pos3_triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;


    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

  for (i = 0; i < Triple->nump; i++)
      Triple->numb[i] = 0;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell3[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt3[i] = -1;

  Triple->tot = 0;
  unique_triples = 0;
  total_triples = 0;

  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell1[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell2[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell3[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt1[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt2[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt3[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell1[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell2[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell3[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt1[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt2[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt3[i] = -1;

  lat_n3 = 0;

      for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
        for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
        //for (atm_n3 = 0; atm_n3 < atoms->number_of_atoms_in_unit_cell; atm_n3++) {

/*
          for (lat_n1 = 0; lat_n1 < R->max_vector; lat_n1++) {
            for (lat_n2 = 0; lat_n2 < 1; lat_n2++) { // FIX
*/
          for (lat_n1 = 0; lat_n1 < 1; lat_n1++) { // FIX
            for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {

            //for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 + R->vec_ai[lat_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - \
             R->vec_ai[lat_n2].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 + R->vec_ai[lat_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - \
             R->vec_ai[lat_n2].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 + R->vec_ai[lat_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - \
             R->vec_ai[lat_n2].comp3;
             Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 + R->vec_ai[lat_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - \
             R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 + R->vec_ai[lat_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - \
             R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 + R->vec_ai[lat_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - \
             R->vec_ai[lat_n3].comp3;
             Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n3].comp1 - \
             R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n3].comp2 - \
             R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n3].comp3 - \
             R->vec_ai[lat_n3].comp3;
             Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
          if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
            //lat_n1 = 0;
            keeper = 1;

            // modified for new order for atm_n1, atm_n2 and atm_n3
            ////triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
            triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n1 * dil1 + lat_n2;
            //triple_index_0 = atm_n3 * dim3 + lat_n3 * dil3 + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;

            //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;
            //fprintf(file.out,"INPUT %10.4lf %10.4lf %10.4lf %10.4lf %5d %3d %3d %3d  %3d %3d %3d\n", \
            Rsqrd12,Rsqrd13,Rsqrd23,R->cutoff * R->cutoff,triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
            //fflush(file.out);

            for (pm = 0; pm < number_of_permutations; pm++) {
              for (k = 0; k < symmetry->number_of_operators; k++) {

                cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O1_temp];
                latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O2_temp];
                latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O3_temp];

                //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                latt3_temp);
                //fflush(file.out);

                switch (pm) {

                  case 0:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = latt1_temp;
                   rotated_latt2 = latt2_temp;
                   rotated_latt3 = latt3_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   break;

                   // for reversed routine swap cases 1 and 2 for permutations
                  case 1:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = latt2_temp;
                   rotated_latt2 = latt1_temp;
                   rotated_latt3 = latt3_temp;
                   //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_cell1 = cell1_temp;
                   //rotated_cell2 = cell3_temp;
                   //rotated_cell3 = cell2_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 2:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 3:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   break;

                  case 4:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  case 5:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  } // close switch

                   // modified for new order for atm_n1, atm_n2 and atm_n3
                   //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2;
                   triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt1 * dil1 + rotated_latt2;
                   //triple_index_1 = rotated_cell3 * dim3 + rotated_latt3 * dil3 + rotated_cell1 * dim1 + \
                   rotated_cell2 * dim2 + rotated_latt2;
                   //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + \
                   rotated_cell3 * dim3 + rotated_latt3;
                   //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                   pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                   rotated_latt2,rotated_latt3);

                 if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                   //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d  latt  %3d %3d %3d\n",\
                   i,pm,k,rotated_cell1,rotated_cell2,\
                   rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                      } // close loop over k
                     } // close loop over pm

                   if (!keeper) continue;

            for (pm = 0; pm < number_of_permutations; pm++) {
              for (k = 0; k < symmetry->number_of_operators; k++) {

                cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];

                if (cell3_temp != atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];

                latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O1_temp];
                latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O2_temp];
                latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O3_temp];
                //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                //fflush(file.out);

                switch (pm) {

                  case 0:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = latt1_temp;
                   rotated_latt2 = latt2_temp;
                   rotated_latt3 = latt3_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   break;

                   // for reversed routine swap cases 1 and 2 for permutations
                  case 1:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = latt2_temp;
                   rotated_latt2 = latt1_temp;
                   rotated_latt3 = latt3_temp;
                   //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_cell1 = cell1_temp;
                   //rotated_cell2 = cell3_temp;
                   //rotated_cell3 = cell2_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 2:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 3:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   break;

                  case 4:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  case 5:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  } // close switch

        //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d   %3d %3d %3d   %3d %3d\n", \
        pm,k,triple_index_0,triple_index_1,\
        rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);

                   for (j = 0; j <= total_triples; j++) {
                     if (rotated_cell1 == Triple->cell1[j] && rotated_latt1 == Triple->latt1[j] && \
                         rotated_cell2 == Triple->cell2[j] && rotated_latt2 == Triple->latt2[j] && \
                         rotated_cell3 == Triple->cell3[j] && rotated_latt3 == Triple->latt3[j])
                         break;
                   if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                        (rotated_cell1 != Triple->cell1[j] || rotated_latt1 != Triple->latt1[j] || \
                         rotated_cell2 != Triple->cell2[j] || rotated_latt2 != Triple->latt2[j] || \
                         rotated_cell3 != Triple->cell3[j] || rotated_latt3 != Triple->latt3[j])) {
                         Triple->cell1[j] = rotated_cell1;
                         Triple->cell2[j] = rotated_cell2;
                         Triple->cell3[j] = rotated_cell3;
                         Triple->latt1[j] = rotated_latt1;
                         Triple->latt2[j] = rotated_latt2;
                         Triple->latt3[j] = rotated_latt3;
                         Triple->k[j] = k;
                         Triple->p[j] = pm;
                        (Triple->numb[unique_triples])++;
                         total_triples++;
              //fprintf(file.out,"Triples %3d %3d %3d %3d pm %3d op %2d  %3d %3d %3d  %3d %3d %3d %3d %3d %3d  %3d %3d %3d  %3d\n",\
                         j,unique_triples+1,total_triples,Triple->numb[unique_triples],pm,k, \
                         rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                         rotated_latt1,rotated_latt2,rotated_latt3,Triple->latt1[j],Triple->latt2[j],Triple->latt3[j],keeper); 
              //printf("Triples %3d %3d pm %3d op %2d    %3d %3d %3d       %3d %3d %3d        %3d\n",j,Triple->tot,Triple->p[j],k, \
                         rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],keeper); 
                         break;
                        }
                       } // close loop over j
                      } // close loop over k
                     } // close loop over pm
                      Triple->posn[unique_triples] = Triple->tot;
                      unique_triples++;
                      Triple->tot = total_triples;
                     } // close if (Rsqrd >
                    //} // close loop over lat_n3
                   } // close loop over lat_n2
                  } // close loop over lat_n1
                  //} // close loop over atm_n3
                 } // close loop over atm_n2
                } // close loop over atm_n1

                  Triple->nump = unique_triples;

                  //if (Triple->tot > 0) fprintf(file.out,"\n");
                  //for(j=0;j<Triple->tot;j++) fprintf(file.out,"gathered e Triples %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d\n",\
                  j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);
                  //for(j=0;j<Triple->nump;j++) fprintf(file.out,"%3d posn %3d numb %3d\n",\
                  j,Triple->posn[j],Triple->numb[j]);
                  //printf("\n");
                  //for(j=0;j<Triple->tot;j++) printf("gathered e Triples %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d\n",\
                  j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);

}

void count_triples_q(int *atm_n3, int *q1, int *q_tran, TRIPLE_TRAN *Triple, PAIR_TRAN *pair_p, KPOINT_TRAN *knet, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3, rotated_qvec;
int q, qvec;
int dil2 = atom_p->numb[*atm_n3];
int dim2 = R_tables->last_vector * dil2;
int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;
int dimq = knet->num[*q1] * dim1;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples;
int q_rot[MXT];
TRIPLE_TMP triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;

  for (i = 0; i < MXT; i++) triples.cell1[i] = -1;
  for (i = 0; i < MXT; i++) triples.cell2[i] = -1;
  for (i = 0; i < MXT; i++) triples.cell3[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt1[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt2[i] = -1;
  for (i = 0; i < MXT; i++) triples.latt3[i] = -1;

  if (job->pms == 0) number_of_permutations = 1;
  else               number_of_permutations = 2;

  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell1[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell2[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.cell3[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt1[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt2[i] = -1;
  for (i = 0;i < number_of_permutations * symmetry->number_of_operators; i++) triples.latt3[i] = -1;

   Triple->tot = 0;
   Triple->nump = 0;
   unique_triples = 0;
   total_triples = 0;

   lat_n3 = 0;
   for (qvec = 0; qvec < knet->num[*q1]; qvec++) {
     q = knet->bz[qvec];
   for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
     for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
       //for (*atm_n3 = 0; *atm_n3 < atoms->number_of_atoms_in_unit_cell; *atm_n3++) {
         for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
           //for (lat_n3 = 0; lat_n3 < R->max_vector; lat_n3++) {
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
             Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[*atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[*atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[*atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
             Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
             Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[*atm_n3].comp1 - \
             R->vec_ai[lat_n3].comp1;
             Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[*atm_n3].comp2 - \
             R->vec_ai[lat_n3].comp2;
             Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[*atm_n3].comp3 - \
             R->vec_ai[lat_n3].comp3;
             Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
           //fprintf(file.out,"INPUT  %5d %3d %3d %3d  %3d %3d %3d\n",triple_index_0,atm_n1,atm_n2,*atm_n3,lat_n1,lat_n2,lat_n3);
           //fflush(file.out);
           if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
           //if (Rsqrd12 < R->cutoff * R->cutoff) {
             lat_n1 = 0;
             keeper = 1;

             // modified for new order for atm_n1, atm_n2 and atm_n3
             ////triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
             //triple_index_0 = *atm_n3 * dim3 + lat_n3 * dil3 + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
             triple_index_0 = q * dimq + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + *atm_n3;
             //printf("qvec %3d n1 %3d n2 %3d l2 %3d n3 %3d q %3d\n",qvec,atm_n1,atm_n2,lat_n2,*atm_n3,q);

             for (pm = 0; pm < number_of_permutations; pm++) {
               for (k = 0; k < symmetry->number_of_operators; k++) {

                 cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                 cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                 cell3_temp = atom_p->K[*atm_n3 * symmetry->number_of_operators + k];
                 rotated_qvec = q_tran[qvec * symmetry->number_of_operators + k];
                 //printf("op %3d cell %3d %3d %3d q %3d %3d\n",k,cell1_temp,cell2_temp,cell3_temp,q,rotated_qvec);

                 if (cell3_temp != *atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                 O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                 O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                 O3_temp = atom_p->O[*atm_n3 * symmetry->number_of_operators + k];

                 latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O1_temp];
                 latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O2_temp];
                 latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O3_temp];

                 //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                 latt3_temp);
                 //fflush(file.out);

                 switch (pm) {

                   case 0:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    break;
                    // for reversed routine swap cases 1 and 2 for permutations
                   case 1:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 2:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_cell1 = cell2_temp;
                    //rotated_cell2 = cell1_temp;
                    //rotated_cell3 = cell3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 3:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 4:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   case 5:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   } // close switch

                    if (rotated_latt1 != 0 || rotated_latt3 != 0) continue;

                    // modified for new order for atm_n1, atm_n2 and atm_n3
                    ////triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2;
                    triple_index_1 = rotated_qvec * dimq + rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + *atm_n3;
                    //triple_index_1 = rotated_cell3 * dim3 + rotated_latt3 * dil3 + rotated_cell1 * dim1 + rotated_cell2 * dim2 + \
                    rotated_latt2;
                    //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + rotated_cell3 * dim3 + \
                    rotated_latt3;
                    //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d q q_rot  %5d %5d\n", \
                    pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                    rotated_latt2,rotated_latt3,qvec,rotated_qvec);

                 if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                    //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell %3d %3d %3d latt  %3d %3d %3d\n",i,pm,k,rotated_cell1,\
                    rotated_cell2,\
                    rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                       } // close loop over k
                      } // close loop over pm

                    if (!keeper) continue;

             for (pm = 0; pm < number_of_permutations; pm++) {
               for (k = 0; k < symmetry->number_of_operators; k++) {

                 cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                 cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                 cell3_temp = atom_p->K[*atm_n3 * symmetry->number_of_operators + k];
                 rotated_qvec = q_tran[qvec * symmetry->number_of_operators + k];

                 if (cell3_temp != *atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                 O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                 O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                 O3_temp = atom_p->O[*atm_n3 * symmetry->number_of_operators + k];

                 latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O1_temp];
                 latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O2_temp];
                 latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                 R_tables->margin_vector + O3_temp];

                 //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                 latt3_temp);
                 //fflush(file.out);

                 switch (pm) {

                   case 0:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    break;

                    // for reversed routine swap cases 1 and 2 for permutations
                   case 1:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell3_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 2:
                    rotated_cell1 = cell1_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                    //rotated_cell1 = cell2_temp;
                    //rotated_cell2 = cell1_temp;
                    //rotated_cell3 = cell3_temp;
                    //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 3:
                    rotated_cell1 = cell2_temp;
                    rotated_cell2 = cell3_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                    break;

                   case 4:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell1_temp;
                    rotated_cell3 = cell2_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   case 5:
                    rotated_cell1 = cell3_temp;
                    rotated_cell2 = cell2_temp;
                    rotated_cell3 = cell1_temp;
                    rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                    rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                    break;

                   } // close switch

                    if (rotated_latt1 != 0 || rotated_latt3 != 0) continue;
         //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d    %3d %3d %3d   %3d %3d %3d q %3d  %3d %3d\n",\
         pm,k,triple_index_0,triple_index_1,\
         rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,rotated_qvec,Triple->tot,total_triples);

                    for (j = 0; j <= total_triples; j++) {
                      if (rotated_cell1 == triples.cell1[j] && rotated_latt1 == triples.latt1[j] && \
                          rotated_cell2 == triples.cell2[j] && rotated_latt2 == triples.latt2[j] && \
                          rotated_cell3 == triples.cell3[j] && rotated_latt3 == triples.latt3[j] && rotated_qvec == q_rot[j])
                          break;
                    if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                         (rotated_cell1 != triples.cell1[j] || rotated_latt1 != triples.latt1[j] || \
                          rotated_cell2 != triples.cell2[j] || rotated_latt2 != triples.latt2[j] || \
                          rotated_cell3 != triples.cell3[j] || rotated_latt3 != triples.latt3[j] || rotated_qvec != q_rot[j])) {
                          triples.cell1[j] = rotated_cell1;
                          triples.cell2[j] = rotated_cell2;
                          triples.cell3[j] = rotated_cell3;
                          triples.latt1[j] = rotated_latt1;
                          triples.latt2[j] = rotated_latt2;
                          triples.latt3[j] = rotated_latt3;
                          q_rot[j] = rotated_qvec;
                          total_triples++;
      //fprintf(file.out,"triples %3d %3d  %3d pm %3d op %2d cell1-3  %3d %3d %3d   %3d %3d %3d latt 1-3 %3d %3d %3d  %3d %3d %3d q  %3d\n",\
                          j,unique_triples+1,total_triples, pm, k, \
                          rotated_cell1,rotated_cell2,rotated_cell3,triples.cell1[j],triples.cell2[j],triples.cell3[j],\
                          rotated_latt1,rotated_latt2,rotated_latt3,triples.latt1[j],triples.latt2[j],triples.latt3[j],\
                          rotated_qvec); 
                          break;
                         }
                        } // close loop over j
                       } // close loop over k
                      } // close loop over pm
                       unique_triples++;
                       //fprintf(file.out,"%3d %3d\n",unique_triples,total_triples);
                       Triple->tot = total_triples;
                      } // close if (Rsqrd > 
                     //} // close loop over lat_n3
                    } // close loop over lat_n2
                   //} // close loop over atm_n3
                  } // close loop over atm_n2
                 } // close loop over atm_n1
                } // close loop over qvec

                 Triple->nump = unique_triples;

}

void generate_triples_q(int *atm_n3, int *q1, int *q_tran, int *q_rot, TRIPLE_TRAN *Triple, PAIR_TRAN *pair_p, KPOINT_TRAN *knet, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int atm_n1, atm_n2, lat_n1, lat_n2, lat_n3;
int cell1_temp, cell2_temp, cell3_temp, O1_temp, O2_temp, O3_temp;
int latt1_temp, latt2_temp, latt3_temp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_latt3, rotated_qvec;
int q, qvec;
int dil2 = atom_p->numb[*atm_n3];
int dim2 = R_tables->last_vector * dil2;
int dim1 = atoms->number_of_atoms_in_unit_cell * dim2;
int dimq = knet->num[*q1] * dim1;
int keeper, number_of_permutations;
int triple_index_0, triple_index_1;
int unique_triples, total_triples;
double Rsqrd12, Rsqrd13, Rsqrd23;
VECTOR_DOUBLE Rvec_tmp;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 2;

  for (i = 0; i < Triple->nump; i++)
      Triple->numb[i] = 0;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt1[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt2[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->cell3[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      Triple->latt3[i] = -1;
  for (i = 0; i < Triple->tot; i++)
      q_rot[i] = -1;

  Triple->tot = 0;
  unique_triples = 0;
  total_triples = 0;

  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell1[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell2[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->cell3[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt1[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt2[i] = -1;
  //for (i=0;i<number_of_permutations * symmetry->number_of_operators;i++) Triple->latt3[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell1[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell2[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->cell3[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt1[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt2[i] = -1;
  for (i=0;i<Triple->tot;i++) Triple->latt3[i] = -1;

  lat_n3 = 0;

   for (qvec = 0; qvec < knet->num[*q1]; qvec++) {
     q = knet->bz[qvec];
      for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
        for (atm_n2 = 0; atm_n2 < atoms->number_of_atoms_in_unit_cell; atm_n2++) {
          for (lat_n2 = 0; lat_n2 < R->max_vector; lat_n2++) {
            Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
            Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
            Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
            Rsqrd12 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
            Rvec_tmp.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[*atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
            Rvec_tmp.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[*atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
            Rvec_tmp.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[*atm_n3].comp3 - R->vec_ai[lat_n3].comp3;
            Rsqrd13 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
            Rvec_tmp.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[*atm_n3].comp1 - \
            R->vec_ai[lat_n3].comp1;
            Rvec_tmp.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[*atm_n3].comp2 - \
            R->vec_ai[lat_n3].comp2;
            Rvec_tmp.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[*atm_n3].comp3 - \
            R->vec_ai[lat_n3].comp3;
            Rsqrd23 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
          if (Rsqrd12 < R->cutoff * R->cutoff && Rsqrd13 < R->cutoff * R->cutoff && Rsqrd23 < R->cutoff * R->cutoff) {
          //if (Rsqrd12 < R->cutoff * R->cutoff) {
            lat_n1 = 0;
            keeper = 1;

            // modified for new order for atm_n1, atm_n2 and atm_n3
            //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;
            triple_index_0 = q * dimq + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + *atm_n3;
            //triple_index_0 = atm_n3 * dim3 + lat_n3 * dil3 + atm_n1 * dim1 + atm_n2 * dim2 + lat_n2;

            //triple_index_0 = atm_n1 * dim1 + atm_n2 * dim2 + lat_n2 * dil2 + atm_n3 * dim3 + lat_n3;
            //fprintf(file.out,"INPUT %10.4lf %10.4lf %10.4lf %10.4lf %5d %3d %3d %3d  %3d %3d %3d\n", \
            Rsqrd12,Rsqrd13,Rsqrd23,R->cutoff * R->cutoff,triple_index_0,atm_n1,atm_n2,atm_n3,lat_n1,lat_n2,lat_n3);
            //fflush(file.out);

            for (pm = 0; pm < number_of_permutations; pm++) {
              for (k = 0; k < symmetry->number_of_operators; k++) {

                cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                cell3_temp = atom_p->K[*atm_n3 * symmetry->number_of_operators + k];
                rotated_qvec = q_tran[qvec * symmetry->number_of_operators + k];

                if (cell3_temp != *atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                O3_temp = atom_p->O[*atm_n3 * symmetry->number_of_operators + k];

                latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O1_temp];
                latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O2_temp];
                latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O3_temp];

                //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,\
                latt3_temp);
                //fflush(file.out);

                switch (pm) {

                  case 0:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   break;

                   // for reversed routine swap cases 1 and 2 for permutations
                  case 1:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_cell1 = cell1_temp;
                   //rotated_cell2 = cell3_temp;
                   //rotated_cell3 = cell2_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 2:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 3:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   break;

                  case 4:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  case 5:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  } // close switch

                   if (rotated_latt1 != 0 || rotated_latt3 != 0) continue;

                   // modified for new order for atm_n1, atm_n2 and atm_n3
                   //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2;
                   triple_index_1 = rotated_qvec * dimq + rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + *atm_n3;
                   //triple_index_1 = rotated_cell3 * dim3 + rotated_latt3 * dil3 + rotated_cell1 * dim1 + \
                   rotated_cell2 * dim2 + rotated_latt2;
                   //triple_index_1 = rotated_cell1 * dim1 + rotated_cell2 * dim2 + rotated_latt2 * dil2 + \
                   rotated_cell3 * dim3 + rotated_latt3;
                   //fprintf(file.out,"UNIQ pm %3d op %3d  indices %5d %5d      %5d %5d %5d    %5d %5d %5d\n", \
                   pm,k,triple_index_0,triple_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,\
                   rotated_latt2,rotated_latt3);

                 if (triple_index_1 < triple_index_0) { keeper = 0; pm = number_of_permutations; k = symmetry->number_of_operators; }

                   //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d  latt  %3d %3d %3d\n",\
                   i,pm,k,rotated_cell1,rotated_cell2,\
                   rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3); fflush(file.out);

                      } // close loop over k
                     } // close loop over pm

                   if (!keeper) continue;

            for (pm = 0; pm < number_of_permutations; pm++) {
              for (k = 0; k < symmetry->number_of_operators; k++) {

                cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
                cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
                cell3_temp = atom_p->K[*atm_n3 * symmetry->number_of_operators + k];
                rotated_qvec = q_tran[qvec * symmetry->number_of_operators + k];

                if (cell3_temp != *atm_n3) continue; // only retain rotations which leave atm_n3 fixed

                O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
                O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
                O3_temp = atom_p->O[*atm_n3 * symmetry->number_of_operators + k];

                latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O1_temp];
                latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O2_temp];
                latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * \
                R_tables->margin_vector + O3_temp];
                //fprintf(file.out,"OUTPUT %3d %3d %3d  %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,latt1_temp,latt2_temp,latt3_temp);
                //fflush(file.out);

                switch (pm) {

                  case 0:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   break;

                   // for reversed routine swap cases 1 and 2 for permutations
                  case 1:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell3_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   //rotated_cell1 = cell1_temp;
                   //rotated_cell2 = cell3_temp;
                   //rotated_cell3 = cell2_temp;
                   //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 2:
                   rotated_cell1 = cell1_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
                   break;

                  case 3:
                   rotated_cell1 = cell2_temp;
                   rotated_cell2 = cell3_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
                   break;

                  case 4:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell1_temp;
                   rotated_cell3 = cell2_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  case 5:
                   rotated_cell1 = cell3_temp;
                   rotated_cell2 = cell2_temp;
                   rotated_cell3 = cell1_temp;
                   rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
                   rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
                   break;

                  } // close switch

                   if (rotated_latt1 != 0 || rotated_latt3 != 0) continue;

        //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d   %3d %3d %3d   %3d %3d\n", \
        pm,k,triple_index_0,triple_index_1,\
        rotated_cell1,rotated_cell2,rotated_cell3,rotated_latt1,rotated_latt2,rotated_latt3,Triple->tot,total_triples);

                   for (j = 0; j <= total_triples; j++) {
                     if (rotated_cell1 == Triple->cell1[j] && rotated_latt1 == Triple->latt1[j] && \
                         rotated_cell2 == Triple->cell2[j] && rotated_latt2 == Triple->latt2[j] && \
                         rotated_cell3 == Triple->cell3[j] && rotated_latt3 == Triple->latt3[j] && rotated_qvec == q_rot[j])
                         break;
                   if (j == total_triples && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && \
                        (rotated_cell1 != Triple->cell1[j] || rotated_latt1 != Triple->latt1[j] || \
                         rotated_cell2 != Triple->cell2[j] || rotated_latt2 != Triple->latt2[j] || \
                         rotated_cell3 != Triple->cell3[j] || rotated_latt3 != Triple->latt3[j] || rotated_qvec != q_rot[j])) {
                         Triple->cell1[j] = rotated_cell1;
                         Triple->cell2[j] = rotated_cell2;
                         Triple->cell3[j] = rotated_cell3;
                         Triple->latt1[j] = rotated_latt1;
                         Triple->latt2[j] = rotated_latt2;
                         Triple->latt3[j] = rotated_latt3;
                         q_rot[j] = rotated_qvec;
                         Triple->k[j] = k;
                         Triple->p[j] = pm;
                        (Triple->numb[unique_triples])++;
                         total_triples++;
              //fprintf(file.out,"Triples %3d %3d %3d %3d pm %3d op %2d  %3d %3d %3d  %3d %3d %3d %3d %3d %3d  %3d %3d %3d  %3d\n",\
                         j,unique_triples+1,total_triples,Triple->numb[unique_triples],pm,k, \
                         rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                         rotated_latt1,rotated_latt2,rotated_latt3,Triple->latt1[j],Triple->latt2[j],Triple->latt3[j],keeper); 
              //printf("Triples %3d %3d pm %3d op %2d    %3d %3d %3d       %3d %3d %3d        %3d\n",j,total_triples,Triple->p[j],k, \
                         rotated_cell1,rotated_cell2,rotated_cell3,Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],q_rot[j]); 
                         break;
                        }
                       } // close loop over j
                      } // close loop over k
                     } // close loop over pm
                      Triple->posn[unique_triples] = Triple->tot;
                      unique_triples++;
                      Triple->tot = total_triples;
                     } // close if (Rsqrd >
                    //} // close loop over lat_n3
                   } // close loop over lat_n2
                  //} // close loop over atm_n3
                 } // close loop over atm_n2
                } // close loop over atm_n1
               } // close loop over qvec

                  Triple->nump = unique_triples;

                  //if (Triple->tot > 0) fprintf(file.out,"\n");
                  //for(j=0;j<Triple->tot;j++) fprintf(file.out,"gathered e Triples %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d\n",\
                  j,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  Triple->latt1[j],Triple->latt2[j],Triple->latt3[j]);
                  //for(j=0;j<Triple->nump;j++) fprintf(file.out,"%3d posn %3d numb %3d\n",\
                  j,Triple->posn[j],Triple->numb[j]);
                  //printf("\n");
                  j = 0;
                  for(i=0;i<Triple->nump;i++) { for(k=0;k<Triple->numb[i];k++) { 
                  printf("gathered e Triples %3d %3d pm %3d op %3d  %3d %3d %3d   %3d %3d %3d q_rot  %3d\n",\
                  i,k,Triple->p[j],Triple->k[j],Triple->cell1[j],Triple->cell2[j],Triple->cell3[j],\
                  Triple->latt1[j],Triple->latt2[j],Triple->latt3[j],q_rot[j]); j++; }
                  printf("\n"); }

}

void generate_molecule_quads(PAIR_TRAN *pair_p, QUAD_TRAN *Quad, int atm_n1, int atm_n2, int atm_n3, int atm_n4, int lat_n1, int lat_n2, int lat_n3, int lat_n4, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int cell1_temp, cell2_temp, cell3_temp, cell4_temp, O1_temp, O2_temp, O3_temp, O4_temp;
int latt_tmp, latt_temp, latt1_temp, latt2_temp, latt3_temp, latt4_temp, latt3_tmp, latt4_tmp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_cell4, rotated_latt3, rotated_latt4;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int dim3 = dim1 * dim2;
int keeper, number_of_permutations;
int quad_index_0, quad_index_1;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 8;

     Quad->tot = 0;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell1[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell2[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell3[i] = -1;
     for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell4[i] = -1;

     quad_index_0 = atm_n1 * dim3 + atm_n2 * dim2 + atm_n3 * dim1 + atm_n4;

     //fprintf(file.out,"INPUT  %3d %3d %3d %3d  %3d %3d %3d %3d\n",atm_n1, atm_n2,atm_n3,atm_n4,lat_n1,lat_n2,lat_n3,lat_n4);
     //fflush(file.out);

    //for (pm = 0; pm < 8; pm++) {
    for (pm = 0; pm < number_of_permutations; pm++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {

      keeper = 0;

      cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
      cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
      cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];
      cell4_temp = atom_p->K[atm_n4 * symmetry->number_of_operators + k];

      O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
      O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
      O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];
      O4_temp = atom_p->O[atm_n4 * symmetry->number_of_operators + k];

      //latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->max_vector+O1_temp];
      //latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->max_vector+O2_temp];
      //latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->max_vector+O3_temp];
      //latt4_temp = R_tables->diffvec[R_tables->lattvec[lat_n4 * symmetry->number_of_operators + k] * R_tables->max_vector+O4_temp];
      latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->max_vector + O1_temp];
      latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->max_vector + O2_temp];
      latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->max_vector + O3_temp];
      latt4_temp = R_tables->diffvec[R_tables->lattvec[lat_n4 * symmetry->number_of_operators + k] * R_tables->max_vector + O4_temp];

     //fprintf(file.out,"OUTPUT %3d %3d %3d %3d  %3d %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,cell4_temp,latt1_temp,latt2_temp,latt3_temp,latt4_temp);
     //fflush(file.out);

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell3_temp;
            rotated_cell4 = cell4_temp;
            //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt1_temp];
            //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt2_temp];
            //rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt3_temp];
            //rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt4_temp];
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt1_temp];
            rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell4_temp;
            rotated_cell4 = cell3_temp;
            //rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt1_temp];
            //rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt2_temp];
            //rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt4_temp];
            //rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt3_temp];
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt1_temp];
            rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt1_temp];
            break;

           case 2:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_cell3 = cell3_temp;
            rotated_cell4 = cell4_temp;
            //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt2_temp];
            //rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt1_temp];
            //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt3_temp];
            //rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt4_temp];
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt2_temp];
            rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt2_temp];
            break;

           case 3:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_cell3 = cell4_temp;
            rotated_cell4 = cell3_temp;
            //rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt2_temp];
            //rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt1_temp];
            //rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt4_temp];
            //rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt3_temp];
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt2_temp];
            rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt2_temp];
            break;

           case 4:
            rotated_cell1 = cell3_temp;
            rotated_cell2 = cell4_temp;
            rotated_cell3 = cell1_temp;
            rotated_cell4 = cell2_temp;
            //rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt3_temp];
            //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt4_temp];
            //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt1_temp];
            //rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt2_temp];
            rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt3_temp];
            rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt3_temp];
            rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt3_temp];
            rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt3_temp];
            break;

           case 5:
            rotated_cell1 = cell3_temp;
            rotated_cell2 = cell4_temp;
            rotated_cell3 = cell2_temp;
            rotated_cell4 = cell1_temp;
            //rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt3_temp];
            //rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt4_temp];
            //rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt2_temp];
            //rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt1_temp];
            rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt3_temp];
            rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt3_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt3_temp];
            rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt3_temp];
            break;

           case 6:
            rotated_cell1 = cell4_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell1_temp;
            rotated_cell4 = cell2_temp;
            //rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt4_temp];
            //rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt3_temp];
            //rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt1_temp];
            //rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt2_temp];
            rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt4_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt4_temp];
            rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt4_temp];
            rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt4_temp];
            break;

           case 7:
            rotated_cell1 = cell4_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell2_temp;
            rotated_cell4 = cell1_temp;
            //rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt4_temp];
            //rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt3_temp];
            //rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt2_temp];
            //rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt1_temp];
            rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->max_vector + latt4_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->max_vector + latt4_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->max_vector + latt4_temp];
            rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->max_vector + latt4_temp];
            break;

           } // close switch

            quad_index_1 = rotated_cell1 * dim3 + rotated_cell2 * dim2 + rotated_cell3 * dim1 + rotated_cell4;

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %3d %3d      %3d %3d %3d %3d    %3d %3d %3d %3d    %3d\n",\
           pm,k,quad_index_0,quad_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,rotated_latt1,\
           rotated_latt2,rotated_latt3,rotated_latt4,pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2]);
           //fflush(file.out);
           // needed for F_ijg
           if ((pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2] == -1 || \
           pair_p->uniq[rotated_latt3 * dim2 + rotated_cell1 * dim1 + rotated_cell3] == -1) && quad_index_1 >= quad_index_0) 
           keeper = 1; 
  
           if (job->type == 1) keeper = 1; // do all integrals for BSE calculation

           if (quad_index_1 < quad_index_0) { Quad->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }
           //if (quad_index_1 < quad_index_0) { Quad->tot = -1; pm = 8; k = symmetry->number_of_operators; }

            //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d %3d latt %3d %3d %3d %3d\n",i,pm,k,rotated_cell1,\
            rotated_cell2,rotated_cell3,rotated_cell4,rotated_latt1,rotated_latt2,rotated_latt3,rotated_latt4); fflush(file.out);
            //fflush(file.out);

            if (!keeper) continue;

            for (j = 0; j <= Quad->tot; j++) {
              if (rotated_cell1 == Quad->cell1[j] && rotated_cell2 == Quad->cell2[j] && \
                  rotated_latt1 == Quad->latt1[j] && rotated_latt2 == Quad->latt2[j] && \
                  rotated_cell3 == Quad->cell3[j] && rotated_cell4 == Quad->cell4[j] && \
                  rotated_latt3 == Quad->latt3[j] && rotated_latt4 == Quad->latt4[j])
                  break;
              if (j == Quad->tot && \
                 (rotated_cell1 != Quad->cell1[j] || rotated_cell2 != Quad->cell2[j] || \
                  rotated_latt1 != Quad->latt1[j] || rotated_latt2 != Quad->latt2[j] || \
                  rotated_cell3 != Quad->cell3[j] || rotated_cell4 != Quad->cell4[j] || \
                  rotated_latt3 != Quad->latt3[j] || rotated_latt4 != Quad->latt4[j])) {
                  Quad->cell1[j] = rotated_cell1;
                  Quad->cell2[j] = rotated_cell2;
                  Quad->cell3[j] = rotated_cell3;
                  Quad->cell4[j] = rotated_cell4;
                  Quad->latt1[j] = rotated_latt1;
                  Quad->latt2[j] = rotated_latt2;
                  Quad->latt3[j] = rotated_latt3;
                  Quad->latt4[j] = rotated_latt4;
                  Quad->k[j] = k;
                  Quad->p[j] = pm;
                  (Quad->tot)++;
                  //printf("Quads %3d %3d pm %3d op %2d      %3d %3d %3d %3d       %3d %3d %3d %3d        %3d\n",\
                  j,Quad->tot,Quad->p[j],k, \
                  rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,\
                  Quad->cell1[j],Quad->cell2[j],Quad->cell1[j],Quad->cell2[j],keeper); 
                  break;
                 }
                } // close loop over j
               } // close loop over k
              } // close loop over pm

              //if (Quad->tot > 0) fprintf(file.out,"\n");
              //for(j = 0; j < Quad->tot;j++) 
              //fprintf(file.out,"gathered m Quads %3d pm %3d op %3d   %3d %3d %3d %3d   %3d %3d %3d %3d\n",\
              j,Quad->p[j],Quad->k[j],Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],\
              Quad->latt1[j],Quad->latt2[j],Quad->latt3[j],Quad->latt4[j]);

}

void generate_c_quads(PAIR_TRAN *pair_p, QUAD_TRAN *Quad, int atm_n1, int atm_n2, int atm_n3, int atm_n4, int lat_n1, int lat_n2, int lat_n3, int lat_n4, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int cell1_temp, cell2_temp, cell3_temp, cell4_temp, O1_temp, O2_temp, O3_temp, O4_temp;
int latt_tmp, latt_temp, latt1_temp, latt2_temp, latt3_temp, latt4_temp, latt3_tmp, latt4_tmp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_cell4, rotated_latt3, rotated_latt4;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int keeper, number_of_permutations;
long long dil1 = R_tables->last_vector;
long long dil2 = dim1 * dil1;
long long dil3 = dil1 * dil2;
long long dil4 = dil3 * dim1;
long long dil5 = dil4 * dil1;
long long dil6 = dil5 * dim1;
long long dim2 = dim1 * dim1;
long long quad_index_0, quad_index_1;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 8;

     Quad->tot = 0;

     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->p[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->k[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell1[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell2[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell3[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell4[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt1[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt2[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt3[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt4[i] = -1;

     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell1[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell2[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell3[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->cell4[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->latt1[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->latt2[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->latt3[i] = -1;
     //for (i=0;i<8 * symmetry->number_of_operators;i++) Quad->latt4[i] = -1;

     quad_index_0 = atm_n1 * dil6 + atm_n2 * dil5 + lat_n2 * dil4 + atm_n3 * dil3 + lat_n3 * dil2 + atm_n4 * dil1 + lat_n4;

     //fprintf(file.out,"INPUT  %5d %3d %3d %3d %3d  %3d %3d %3d %3d\n",quad_index_0,atm_n1,atm_n2,atm_n3,atm_n4,lat_n1,lat_n2,lat_n3,lat_n4);

    for (pm = 0; pm < number_of_permutations; pm++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {

      keeper = 0;

      cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
      cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
      cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];
      cell4_temp = atom_p->K[atm_n4 * symmetry->number_of_operators + k];

      O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
      O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
      O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];
      O4_temp = atom_p->O[atm_n4 * symmetry->number_of_operators + k];

   latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
   latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
   latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];
   latt4_temp = R_tables->diffvec[R_tables->lattvec[lat_n4 * symmetry->number_of_operators + k] * R_tables->margin_vector + O4_temp];

     //fprintf(file.out,"OUTPUT %3d %3d %3d %3d  %3d %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,cell4_temp,latt1_temp,latt2_temp,latt3_temp,latt4_temp);

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell3_temp;
            rotated_cell4 = cell4_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell4_temp;
            rotated_cell4 = cell3_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 2:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_cell3 = cell3_temp;
            rotated_cell4 = cell4_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt2_temp];
            break;

           case 3:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_cell3 = cell4_temp;
            rotated_cell4 = cell3_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
            break;

           case 4:
            rotated_cell1 = cell3_temp;
            rotated_cell2 = cell4_temp;
            rotated_cell3 = cell1_temp;
            rotated_cell4 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
            break;

           case 5:
            rotated_cell1 = cell3_temp;
            rotated_cell2 = cell4_temp;
            rotated_cell3 = cell2_temp;
            rotated_cell4 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
            break;

           case 6:
            rotated_cell1 = cell4_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell1_temp;
            rotated_cell4 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt4_temp];
            break;

           case 7:
            rotated_cell1 = cell4_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell2_temp;
            rotated_cell4 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt4_temp];
            break;

           } // close switch

            quad_index_1 = rotated_cell1 * dil6 + rotated_cell2 * dil5 + rotated_latt2 * dil4 + rotated_cell3 * dil3 + \
            rotated_latt3 * dil2 + rotated_cell4 * dil1 + rotated_latt4;

           //fprintf(file.out,"UNIQ pm %3d op %3d  indices %11lli %11lli      %3d %3d %3d %3d    %3d %3d %3d %3d    %3d\n",\
           pm,k,quad_index_0,quad_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,rotated_latt1,\
           rotated_latt2,rotated_latt3,rotated_latt4,pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2]);

           // use these two conditions alone to test algorithm - number of quads should be square of number of pairs
           if (quad_index_1 >= quad_index_0 || rotated_latt3 != 0) keeper = 1;
           if (rotated_latt3 == 0 && quad_index_1 < quad_index_0) { 
           Quad->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }
           //if (rotated_latt3 == 0 && quad_index_1 < quad_index_0) { Quad->tot = -1; pm = 8; k = symmetry->number_of_operators; }

           if ((pm > 0 || k > 0) && pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2] != -1) continue;

           // use this additional condition for quads needed for SCF
           ////keep if (pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2] != -1) continue;

            //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d %3d  latt  %3d %3d %3d %3d keep %3d\n",\
            i,pm,k,rotated_cell1,rotated_cell2,\
            rotated_cell3,rotated_cell4,rotated_latt1,rotated_latt2,rotated_latt3,rotated_latt4,keeper); fflush(file.out);

            //if ((pm > 0 || k > 0) && pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2] != -1) continue;

            if (!keeper) continue;

            if ((rotated_latt4 >= R_tables->last_vector || rotated_latt2 >= R_tables->last_vector) && (pm == 0 && k == 0)) {
                 Quad->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }

            if (rotated_latt4 >= R_tables->last_vector || rotated_latt2 >= R_tables->last_vector) continue;

            for (j = 0; j <= Quad->tot; j++) {
              if (rotated_cell1 == Quad->cell1[j] && rotated_cell2 == Quad->cell2[j] && \
                  rotated_latt1 == Quad->latt1[j] && rotated_latt2 == Quad->latt2[j] && \
                  rotated_cell3 == Quad->cell3[j] && rotated_cell4 == Quad->cell4[j] && \
                  rotated_latt3 == Quad->latt3[j] && rotated_latt4 == Quad->latt4[j])
                  break;
              //if (j == Quad->tot && rotated_latt3 == 0 && 
              if (j == Quad->tot && rotated_latt3 == 0 && rotated_latt2 < R_tables->last_vector && \
              rotated_latt4 < R_tables->last_vector &&
                 (rotated_cell1 != Quad->cell1[j] || rotated_cell2 != Quad->cell2[j] || \
                  rotated_latt1 != Quad->latt1[j] || rotated_latt2 != Quad->latt2[j] || \
                  rotated_cell3 != Quad->cell3[j] || rotated_cell4 != Quad->cell4[j] || \
                  rotated_latt3 != Quad->latt3[j] || rotated_latt4 != Quad->latt4[j])) {
                  Quad->cell1[j] = rotated_cell1;
                  Quad->cell2[j] = rotated_cell2;
                  Quad->cell3[j] = rotated_cell3;
                  Quad->cell4[j] = rotated_cell4;
                  Quad->latt1[j] = rotated_latt1;
                  Quad->latt2[j] = rotated_latt2;
                  Quad->latt3[j] = rotated_latt3;
                  Quad->latt4[j] = rotated_latt4;
                  Quad->k[j] = k;
                  Quad->p[j] = pm;
                  (Quad->tot)++;
                  //fprintf(file.out,"Quads %3d %3d pm %3d op %2d      %3d %3d %3d %3d       %3d %3d %3d %3d     %3d %3d %3d %3d   %3d %3d %3d %3d  %3d\n",\
                  j,Quad->tot,Quad->p[j],k, \
                  rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],\
                  rotated_latt1,rotated_latt2,rotated_latt3,rotated_latt4,Quad->latt1[j],Quad->latt2[j],Quad->latt3[j],Quad->latt4[j],keeper); 
                  break;
                 }
                } // close loop over j
               } // close loop over k
              } // close loop over pm

              //if (Quad->tot > 0) fprintf(file.out,"\n");
      //for (j = 0; j < Quad->tot; j++) fprintf(file.out,"gathered c Quads %3d pm %3d op %3d   %3d %3d %3d %3d   %3d %3d %3d %3d\n",\
              j,Quad->p[j],Quad->k[j],Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],\
              Quad->latt1[j],Quad->latt2[j],Quad->latt3[j],Quad->latt4[j]);

              //if (Quad->tot > 0) printf("\n");
              //for (j = 0; j < Quad->tot; j++) \
              printf("%3d gathered c Quads %3d pm %3d op %3d   %3d %3d %3d %3d   %3d %3d %3d %3d\n",\
              job->taskid,j,Quad->p[j],Quad->k[j],Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],\
              Quad->latt1[j],Quad->latt2[j],Quad->latt3[j],Quad->latt4[j]);

}

void generate_e_quads(PAIR_TRAN *pair_p, QUAD_TRAN *Quad, int atm_n1, int atm_n2, int atm_n3, int atm_n4, int lat_n1, int lat_n2, int lat_n3, int lat_n4, ATOM *atoms, ATOM_TRAN *atom_p, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
  
int i, j, k, pm;
int cell1_temp, cell2_temp, cell3_temp, cell4_temp, O1_temp, O2_temp, O3_temp, O4_temp;
int latt_tmp, latt_temp, latt1_temp, latt2_temp, latt3_temp, latt4_temp, latt3_tmp, latt4_tmp;
int rotated_cell1, rotated_cell2, rotated_latt1, rotated_latt2;
int rotated_cell3, rotated_cell4, rotated_latt3, rotated_latt4;
int keeper, number_of_permutations;
long long dil1 = R_tables->last_vector;
long long dim1 = atoms->number_of_atoms_in_unit_cell;
long long dil2 = dim1 * dil1;
long long dil3 = dil1 * dil2;
long long dil4 = dil3 * dim1;
long long dil5 = dil4 * dil1;
long long dil6 = dil5 * dim1;
long long dim2 = dim1 * dim1;
long long quad_index_0, quad_index_1;

    if (job->pms == 0) number_of_permutations = 1;
    else               number_of_permutations = 8;

     Quad->tot = 0;

     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->p[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->k[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell1[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell2[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell3[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->cell4[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt1[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt2[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt3[i] = -1;
     for (i = 0; i < number_of_permutations * symmetry->number_of_operators; i++) Quad->latt4[i] = -1;

     quad_index_0 = atm_n1 * dil6 + atm_n2 * dil5 + lat_n2 * dil4 + atm_n3 * dil3 + lat_n3 * dil2 + atm_n4 * dil1 + lat_n4;

     //fprintf(file.out,"INPUT  %11lli %3d %3d %3d %3d  %3d %3d %3d %3d\n",quad_index_0,atm_n1,atm_n2,atm_n3,atm_n4,\
     lat_n1,lat_n2,lat_n3,lat_n4);
     //fprintf(file.out,"%lli %lli %lli %lli %lli %lli\n",dil1,dil2,dil3,dil4,dil5,dil6);

    //for (pm = 0; pm < 8; pm++) {
    for (pm = 0; pm < number_of_permutations; pm++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {

      keeper = 0;

      cell1_temp = atom_p->K[atm_n1 * symmetry->number_of_operators + k];
      cell2_temp = atom_p->K[atm_n2 * symmetry->number_of_operators + k];
      cell3_temp = atom_p->K[atm_n3 * symmetry->number_of_operators + k];
      cell4_temp = atom_p->K[atm_n4 * symmetry->number_of_operators + k];

      O1_temp = atom_p->O[atm_n1 * symmetry->number_of_operators + k];
      O2_temp = atom_p->O[atm_n2 * symmetry->number_of_operators + k];
      O3_temp = atom_p->O[atm_n3 * symmetry->number_of_operators + k];
      O4_temp = atom_p->O[atm_n4 * symmetry->number_of_operators + k];

   latt1_temp = R_tables->diffvec[R_tables->lattvec[lat_n1 * symmetry->number_of_operators + k] * R_tables->margin_vector + O1_temp];
   latt2_temp = R_tables->diffvec[R_tables->lattvec[lat_n2 * symmetry->number_of_operators + k] * R_tables->margin_vector + O2_temp];
   latt3_temp = R_tables->diffvec[R_tables->lattvec[lat_n3 * symmetry->number_of_operators + k] * R_tables->margin_vector + O3_temp];
   latt4_temp = R_tables->diffvec[R_tables->lattvec[lat_n4 * symmetry->number_of_operators + k] * R_tables->margin_vector + O4_temp];

      //if (latt1_temp > R_tables->last_vector || latt2_temp > R_tables->last_vector || latt3_temp > R_tables->last_vector || \
      latt4_temp > R_tables->last_vector) continue;
      //fprintf(file.out,"OUTPUT %3d %3d %3d %3d  %3d %3d %3d %3d\n",cell1_temp,cell2_temp,cell3_temp,cell4_temp,latt1_temp,\
      latt2_temp,latt3_temp,latt4_temp);

        switch (pm) {

           case 0:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell3_temp;
            rotated_cell4 = cell4_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 1:
            rotated_cell1 = cell1_temp;
            rotated_cell2 = cell2_temp;
            rotated_cell3 = cell4_temp;
            rotated_cell4 = cell3_temp;
            rotated_latt1 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt2 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt1_temp];
            rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt1_temp];
            break;

           case 2:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_cell3 = cell3_temp;
            rotated_cell4 = cell4_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt4 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt2_temp];
            break;

           case 3:
            rotated_cell1 = cell2_temp;
            rotated_cell2 = cell1_temp;
            rotated_cell3 = cell4_temp;
            rotated_cell4 = cell3_temp;
            rotated_latt1 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt2 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt3 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt2_temp];
            rotated_latt4 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt2_temp];
            break;

           case 4:
            rotated_cell1 = cell3_temp;
            rotated_cell2 = cell4_temp;
            rotated_cell3 = cell1_temp;
            rotated_cell4 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
            break;

           case 5:
            rotated_cell1 = cell3_temp;
            rotated_cell2 = cell4_temp;
            rotated_cell3 = cell2_temp;
            rotated_cell4 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt2 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt3_temp];
            rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt3_temp];
            break;

           case 6:
            rotated_cell1 = cell4_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell1_temp;
            rotated_cell4 = cell2_temp;
            rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt3 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt4 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt4_temp];
            break;

           case 7:
            rotated_cell1 = cell4_temp;
            rotated_cell2 = cell3_temp;
            rotated_cell3 = cell2_temp;
            rotated_cell4 = cell1_temp;
            rotated_latt1 = R_tables->diffvec[latt4_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt2 = R_tables->diffvec[latt3_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt3 = R_tables->diffvec[latt2_temp * R_tables->margin_vector + latt4_temp];
            rotated_latt4 = R_tables->diffvec[latt1_temp * R_tables->margin_vector + latt4_temp];
            break;

           } // close switch

            quad_index_1 = rotated_cell1 * dil6 + rotated_cell2 * dil5 + rotated_latt2 * dil4 + rotated_cell3 * dil3 + \
            rotated_latt3 * dil2 + rotated_cell4 * dil1 + rotated_latt4;

           //fprintf(file.out,"UNIQ pm %3d op %3d indices %11lli %11lli  %3d %3d %3d %3d    %3d %3d %3d %3d    %4d %4d %4d   %3d\n",\
           pm,k,quad_index_0,quad_index_1,rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,rotated_latt1,\
           rotated_latt2,rotated_latt3,rotated_latt4,R_tables->diffvec[rotated_latt4 * R_tables->last_vector + rotated_latt2],\
           rotated_latt3,R_tables->last_vector, pair_p->uniq[rotated_latt3 * dim2 + rotated_cell1 * dim1 + rotated_cell3]);

/*
           //generate_molecule_quads conditions
           
           // don't need this if ((pair_p->uniq[rotated_latt2 * dim2 + rotated_cell1 * dim1 + rotated_cell2] == -1 ||
           if ((pair_p->uniq[rotated_latt3 * dim2 + rotated_cell1 * dim1 + rotated_cell3] == -1) && quad_index_1 >= quad_index_0) \
           keeper = 1;
           // needed for F_ijg
  
           if (job->type == 1) keeper = 1; // do all integrals for BSE calculation

           if (quad_index_1 < quad_index_0) { Quad->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }
           //if (quad_index_1 < quad_index_0) { Quad->tot = -1; pm = 8; k = symmetry->number_of_operators; }

            //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d %3d  latt  %3d %3d %3d %3d\n",\
            i,pm,k,rotated_cell1,rotated_cell2,\\
            rotated_cell3,rotated_cell4,rotated_latt1,rotated_latt2,rotated_latt3,rotated_latt4); fflush(file.out);

            if (!keeper) continue;

           // in comments - these were previously conditions for e_quads - use molecule quads above instead
*/

            if (quad_index_1 >= quad_index_0) keeper = 1;
            if (quad_index_1 < quad_index_0) { Quad->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }
            //if (quad_index_1 < quad_index_0) { Quad->tot = -1; pm = 8; k = symmetry->number_of_operators; }

            if ((pm > 0 || k > 0) && pair_p->uniq[rotated_latt3 * dim2 + rotated_cell1 * dim1 + rotated_cell3] != -1) continue;
            //2019 if (pair_p->uniq[rotated_latt3 * dim2 + rotated_cell1 * dim1 + rotated_cell3] != -1) continue;
            //if (pair_p->uniq[rotated_latt3 * dim2 + rotated_cell1 * dim1 + rotated_cell3] == -1 && quad_index_1 >= quad_index_0) \
            keeper = 1;

            //fprintf(file.out,"rotated i %3d pm %3d op %3d  cell  %3d %3d %3d %3d  latt  %3d %3d %3d %3d\n",\
            i,pm,k,rotated_cell1,rotated_cell2,\
            rotated_cell3,rotated_cell4,rotated_latt1,rotated_latt2,rotated_latt3,rotated_latt4); fflush(file.out);

            if (!keeper) continue;
/* temporary
          if ((R_tables->diffvec[rotated_latt4 * R_tables->margin_vector + rotated_latt2] >= R_tables->last_vector || \
               R_tables->diffvec[rotated_latt3 * R_tables->margin_vector + rotated_latt1] >= R_tables->last_vector || \
               rotated_latt2 >= R_tables->last_vector || \
               rotated_latt3 >= R_tables->last_vector || \
               rotated_latt4 >= R_tables->last_vector) && (pm == 0 && k == 0)) {
           //fprintf(file.out,"diffvec %3d %3d    %3d %3d %3d   %3d\n", \
           R_tables->diffvec[rotated_latt4 * R_tables->margin_vector + rotated_latt2], \
           R_tables->diffvec[rotated_latt3 * R_tables->margin_vector + rotated_latt1], \
               rotated_latt2, rotated_latt3, rotated_latt4,R_tables->last_vector);
           Quad->tot = -1; pm = number_of_permutations; k = symmetry->number_of_operators; }

           if (R_tables->diffvec[rotated_latt4 * R_tables->margin_vector + rotated_latt2] >= R_tables->last_vector || \
               R_tables->diffvec[rotated_latt3 * R_tables->margin_vector + rotated_latt1] >= R_tables->last_vector) continue; 
temporary */

           if (R_tables->diffvec[rotated_latt4 * R_tables->margin_vector + rotated_latt2] >= R_tables->last_vector || \
               R_tables->diffvec[rotated_latt3 * R_tables->margin_vector + rotated_latt1] >= R_tables->last_vector || \
               rotated_latt2 >= R_tables->last_vector || \
               rotated_latt3 >= R_tables->last_vector || \
               rotated_latt4 >= R_tables->last_vector) continue;


            for (j = 0; j <= Quad->tot; j++) {
              if (rotated_cell1 == Quad->cell1[j] && rotated_cell2 == Quad->cell2[j] && \
                  rotated_latt1 == Quad->latt1[j] && rotated_latt2 == Quad->latt2[j] && \
                  rotated_cell3 == Quad->cell3[j] && rotated_cell4 == Quad->cell4[j] && \
                  rotated_latt3 == Quad->latt3[j] && rotated_latt4 == Quad->latt4[j])
                  break;
              //if (j == Quad_1->tot && 
            if (j == Quad->tot && rotated_latt2 < R_tables->last_vector && rotated_latt3 < R_tables->last_vector && rotated_latt4 < \
                R_tables->last_vector && \
                 (rotated_cell1 != Quad->cell1[j] || rotated_cell2 != Quad->cell2[j] || \
                  rotated_latt1 != Quad->latt1[j] || rotated_latt2 != Quad->latt2[j] || \
                  rotated_cell3 != Quad->cell3[j] || rotated_cell4 != Quad->cell4[j] || \
                  rotated_latt3 != Quad->latt3[j] || rotated_latt4 != Quad->latt4[j])) {
                  Quad->cell1[j] = rotated_cell1;
                  Quad->cell2[j] = rotated_cell2;
                  Quad->cell3[j] = rotated_cell3;
                  Quad->cell4[j] = rotated_cell4;
                  Quad->latt1[j] = rotated_latt1;
                  Quad->latt2[j] = rotated_latt2;
                  Quad->latt3[j] = rotated_latt3;
                  Quad->latt4[j] = rotated_latt4;
                  Quad->k[j] = k;
                  Quad->p[j] = pm;
                  (Quad->tot)++;
      //fprintf(file.out,"Quads %3d %3d pm %3d op %2d  %3d %3d %3d %3d   %3d %3d %3d %3d  %3d %3d %3d %3d  %3d %3d %3d %3d  %3d\n",\
                  j,Quad->tot,Quad->p[j],k,rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,\
                  Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],rotated_latt1,rotated_latt2,rotated_latt3,\
                  rotated_latt4,Quad->latt1[j],Quad->latt2[j],Quad->latt3[j],Quad->latt4[j],keeper); 
         //printf("Quads %3d %3d pm %3d op %2d      %3d %3d %3d %3d       %3d %3d %3d %3d        %3d\n",j,Quad->tot,Quad->p[j],k, \
                  rotated_cell1,rotated_cell2,rotated_cell3,rotated_cell4,\
                  Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],keeper); 
                  break;
                 }
                } // close loop over j
               } // close loop over k
              } // close loop over pm
            //if (Quad->tot > 0 && (Quad->cell1[0] > 16 || Quad->cell1[0] < -1)) { printf("%3d\n",job->taskid); fprintf(file.out,"\n");

          //for(j=0;j<Quad->tot;j++) fprintf(file.out,"gathered e Quads %3d pm %3d op %3d  %3d %3d %3d %3d  %3d %3d %3d %3d  %3d\n",\
              j,Quad->p[j],Quad->k[j],Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],Quad->latt1[j],Quad->latt2[j],\
              Quad->latt3[j],Quad->latt4[j],pair_p->uniq[Quad->latt3[j] * dim2 + Quad->cell1[j] * dim1 + Quad->cell3[j]]);
              //if (Quad->tot > 0) printf("\n"); \
              for(j=0;j<Quad->tot;j++) printf("gathered e Quads %3d pm %3d op %3d   %3d %3d %3d %3d   %3d %3d %3d %3d\n",\
              j,Quad->p[j],Quad->k[j],Quad->cell1[j],Quad->cell2[j],Quad->cell3[j],Quad->cell4[j],\
              Quad->latt1[j],Quad->latt2[j],Quad->latt3[j],Quad->latt4[j]);

}

void print_pairs(PAIR_TRAN *pair_p, ATOM *atoms, REAL_LATTICE *R, JOB_PARAM *job, FILES file)
    
{   
   
int i, j, q, qi;
int atm1, atm2, latt1, latt2;
double RAB2;
VECTOR_DOUBLE Rvec_tmp;

  if (job->taskid >= 0 && (job->verbosity > 1 || job->print_pairs == 1)) {
    for (i = 0; i < pair_p->nump; i++) {
      q  = pair_p->posn[i];
      for (j = 0; j < pair_p->numb[i]; j++) {
        atm1 = pair_p->cell1[q + j];
        atm2 = pair_p->cell2[q + j];
        latt1 = pair_p->latt1[q + j];
        latt2 = pair_p->latt2[q + j];
        Rvec_tmp.comp1 = atoms->cell_vector[atm1].comp1 + R->vec_ai[latt1].comp1 - atoms->cell_vector[atm2].comp1 - \
        R->vec_ai[latt2].comp1;
        Rvec_tmp.comp2 = atoms->cell_vector[atm1].comp2 + R->vec_ai[latt1].comp2 - atoms->cell_vector[atm2].comp2 - \
        R->vec_ai[latt2].comp2;
        Rvec_tmp.comp3 = atoms->cell_vector[atm1].comp3 + R->vec_ai[latt1].comp3 - atoms->cell_vector[atm2].comp3 - \
        R->vec_ai[latt2].comp3;
        RAB2 = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
   fprintf(file.out,"pair[%5d]   %3d %3d  %3d   transforms to pair[%5d]  %3d %3d  %5d under operator %2d permutation %2d %10.4lf\n",\
        q, pair_p->cell1[q], pair_p->cell2[q], pair_p->latt2[q], q + j, pair_p->cell1[q + j], pair_p->cell2[q + j], \
        pair_p->latt2[q + j],pair_p->k[q + j],pair_p->p[q + j],sqrt(RAB2) * bohr_to_AA);
       }
       fprintf(file.out,"\n");
      }
     }

}
