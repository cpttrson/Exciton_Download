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
#include <cstdlib>
#include "conversion_factors.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "LIMITS.h"
#include "MATRIX_UTIL.h"
#include "SETUP_CRYSTAL.h"

using namespace std;

  // ******************************************************************************************
  // * CRYSTAL structure contains conventional, primitive and reciprocal space lattices       *
  // ******************************************************************************************

void conventional_unit_cell(CRYSTAL *crystal, FILES file)

{

  // ******************************************************************************************
  // * Calculate conventional cell vectors from lattice parameters                            *
  // ******************************************************************************************

  int i ;
  //double det;

  for (i = 0; i < 3; i++) {
    crystal->conventional_cell[i].comp1 = k_zero;
    crystal->conventional_cell[i].comp2 = k_zero;
    crystal->conventional_cell[i].comp3 = k_zero;
  }

  switch (crystal->type[0]) {

    case 'C':

      crystal->conventional_cell[0].comp1 = crystal->lattice_a * sin(crystal->gamma);
      crystal->conventional_cell[0].comp2 = crystal->lattice_a * cos(crystal->gamma);

      crystal->conventional_cell[1].comp2 = crystal->lattice_b;

      crystal->conventional_cell[2].comp1 = crystal->lattice_c * (cos(crystal->beta) - cos(crystal->alpha) * \
      cos(crystal->gamma)) / sin(crystal->gamma);
      crystal->conventional_cell[2].comp2 = crystal->lattice_c * cos(crystal->alpha);
      crystal->conventional_cell[2].comp3 = crystal->lattice_c * sqrt((k_one - (cos(crystal->alpha) * cos(crystal->alpha) \
      + cos(crystal->beta) * cos(crystal->beta) - 2.0 * cos(crystal->alpha) * cos(crystal->beta) * cos(crystal->gamma)) \
      / sin(crystal->gamma) / sin(crystal->gamma)));

      break;

    case 'S':

      crystal->conventional_cell[0].comp1 = crystal->lattice_a * sin(crystal->gamma);
      crystal->conventional_cell[0].comp2 = crystal->lattice_a * cos(crystal->gamma);

      crystal->conventional_cell[1].comp2 = crystal->lattice_b;

      break;

    case 'P':

      crystal->conventional_cell[0].comp3 = crystal->lattice_c;

      break;

    case 'M':

      break;

  } // close switch{crystal->type[0]

}

void primitive_unit_cell(CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Convert conventional cell to primitive cell using crystal->centring  Table 5.1.3.1    *
  // ******************************************************************************************

  int i, j, k, l ;
  int inr1[9], *p_inr1, inr1_inv[9], *p_inr1_inv, inr_temp[9], inr_fac;
  double det ;
  double irr1[9], *p_irr1, irr1_inv[9], *p_irr1_inv;
  double irr_temp[9];

  p_inr1 = inr1;
  p_inr1_inv = inr1_inv;
  p_irr1 = irr1;
  p_irr1_inv = irr1_inv;

  switch (crystal->centring) {

    case 'P': // P -> P

      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;

      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;

      inr_fac = 1;

      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;

      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;

      break;

    case 'A': // A -> P

      *p_inr1 = 2;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;

      *p_inr1_inv = 2;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;

      inr_fac = 2;

      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = -half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;

      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      //*p_irr1_inv = half;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      //*p_irr1_inv = half;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      //*p_irr1_inv = -half;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      //*p_irr1_inv = half;
      p_irr1_inv++;

      break;

    case 'B': // B -> P

      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = -half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;

      inr_fac = 2;

      if (job->taskid == 0)
      fprintf(file.out, "Insert correct p_irr1_inv array line 1157 in MAIN.cpp\n");
      exit(1);

      break;

    case 'C': // C -> P

//modified to agree with crystal on Cc Fe3O4 10/7/2013

      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = -1;
      //*p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      //*p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;

      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;

      inr_fac = 2;

      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = -half;
      //*p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      //*p_irr1 = -half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;

      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = -k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;

      // fprintf(file_out,"Insert correct p_irr1_inv array line 1157 in MAIN.cpp\n") ; exit(1) ;

      break;

    case 'F': // F -> P

      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;

      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;

      inr_fac = 2;

      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;

      *p_irr1_inv = -k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = -k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = -k_one;
      p_irr1_inv++;

      break;

    case 'I': // I -> P

      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;

      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;

      inr_fac = 2;

      *p_irr1 = -half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = -half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = -half;
      p_irr1++;

      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;

      break;

    case 'R': // H -> R

      *p_inr1 = 2;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = -2;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;

      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 0;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;

      inr_fac = 3;

      *p_irr1 = two_thirds;
      p_irr1++;
      *p_irr1 = -third;
      p_irr1++;
      *p_irr1 = -third;
      p_irr1++;
      *p_irr1 = third;
      p_irr1++;
      *p_irr1 = third;
      p_irr1++;
      *p_irr1 = -two_thirds;
      p_irr1++;
      *p_irr1 = third;
      p_irr1++;
      *p_irr1 = third;
      p_irr1++;
      *p_irr1 = third;
      p_irr1++;

      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = -k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = -k_one;
      p_irr1_inv++;
      *p_irr1_inv = k_one;
      p_irr1_inv++;

      //fprintf(file.out, "Insert correct p_irr1_inv array line 1157 in MAIN.cpp\n");
      //exit(1);

/*
//03/11/20 
      // Convert Triple Hexagonal Cell obverse setting to Primitive Rhombohedral Cell R

      *p_inr1_inv = 2;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = -1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = -2;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;
      *p_inr1_inv = 1;
      p_inr1_inv++;

      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;
      *p_inr1 = 0;
      p_inr1++;
      *p_inr1 = -1;
      p_inr1++;
      *p_inr1 = 1;
      p_inr1++;

      inr_fac = 3;

      *p_irr1_inv = two_thirds;
      p_irr1_inv++;
      *p_irr1_inv = -third;
      p_irr1_inv++;
      *p_irr1_inv = -third;
      p_irr1_inv++;
      *p_irr1_inv = third;
      p_irr1_inv++;
      *p_irr1_inv = third;
      p_irr1_inv++;
      *p_irr1_inv = -two_thirds;
      p_irr1_inv++;
      *p_irr1_inv = third;
      p_irr1_inv++;
      *p_irr1_inv = third;
      p_irr1_inv++;
      *p_irr1_inv = third;
      p_irr1_inv++;

      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = -k_one;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = -k_one;
      p_irr1++;
      *p_irr1 = k_one;
      p_irr1++;

      //fprintf(file.out, "Insert correct p_irr1_inv array line 1157 in MAIN.cpp\n");
      //exit(1);
*/

      break;
  }

  for (i = 0; i < 3; i++) {
    crystal->primitive_cell[i].comp1 = k_zero;
    crystal->primitive_cell[i].comp2 = k_zero;
    crystal->primitive_cell[i].comp3 = k_zero;
  }

    for (l = 0; l < symmetry->number_of_operators; l++) {
      for (i = 0; i < 9; i++)
        inr_temp[i] = 0; 
          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                inr_temp[i * 3 + k] += symmetry->inr[9 * l + i * 3 + j] * inr1_inv[j * 3 + k];  
               }
              }
             }
      for (i = 0; i < 9; i++)
        symmetry->inr[9 * l + i] = 0; 
          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                symmetry->inr[9 * l + i * 3 + k] += inr1[i * 3 + j] * inr_temp[j * 3 + k]; 
               }
              }
             }
            for (i=0;i<9;i++) 
            symmetry->inr[9 * l + i] /= inr_fac;
           }

    // Cartesian operators need to be transformed only for M, R, H
//03/11/20   
/*
    if (crystal->system == 'R') {
    for (l = 0; l < symmetry->number_of_operators; l++) {
      for (i = 0; i < 9; i++)
        irr_temp[i] = 0; 
          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                irr_temp[i * 3 + k] += symmetry->irr[9 * l + i * 3 + j] * irr1_inv[j * 3 + k];  
               }
              }
             }
      for (i = 0; i < 9; i++)
        symmetry->irr[9 * l + i] = 0; 
          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                symmetry->irr[9 * l + i * 3 + k] += irr1[i * 3 + j] * irr_temp[j * 3 + k]; 
               }
              }
             }
           }
          }
*/

    if (job->taskid == 0 && job->verbosity >= 1) {
    fprintf(file.out,"Symmetry Operators transformed to primitive cell lattice basis\n");

    fprintf(file.out, "                                %2d SYMMETRY OPERATORS\n\n", symmetry->number_of_operators);
    for (i = 0; i < symmetry->number_of_operators; i++)
          fprintf(file.out,"%3d\n %5d %5d %5d\n %5d %5d %5d\n %5d %5d %5d\n",
          i + 1, symmetry->inr[i * 9 + 0], symmetry->inr[i * 9 + 1], symmetry->inr[i * 9 + 2], \
                 symmetry->inr[i * 9 + 3], symmetry->inr[i * 9 + 4], symmetry->inr[i * 9 + 5], \
                 symmetry->inr[i * 9 + 6], symmetry->inr[i * 9 + 7], symmetry->inr[i * 9 + 8]);
    fprintf(file.out, "\n");

    fprintf(file.out, "                                %2d SYMMETRY OPERATORS\n\n", symmetry->number_of_operators);
    for (i = 0; i < symmetry->number_of_operators; i++)
          fprintf(file.out,"%3d\n %5.2f %5.2f %5.2f\n %5.2f %5.2f %5.2f\n %5.2f %5.2f %5.2f\n",
          i + 1, symmetry->irr[i * 9 + 0], symmetry->irr[i * 9 + 1], symmetry->irr[i * 9 + 2], \
                 symmetry->irr[i * 9 + 3], symmetry->irr[i * 9 + 4], symmetry->irr[i * 9 + 5], \
                 symmetry->irr[i * 9 + 6], symmetry->irr[i * 9 + 7], symmetry->irr[i * 9 + 8]);
    fprintf(file.out, "\n");
   }

  switch (crystal->type[0]) {

    case 'C':

/*
    for (l = 0; l < symmetry->number_of_operators; l++) {
      for (i = 0; i < 9; i++)
        inr_temp[i] = 0; 
          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                inr_temp[i * 3 + k] += symmetry->inr[9 * l + i * 3 + j] * inr1_inv[j * 3 + k];
               }
              }
             }
      for (i = 0; i < 9; i++)
        symmetry->inr[9 * l + i] = 0; 
          for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
              for (k = 0; k < 3; k++) {
                symmetry->inr[9 * l + i * 3 + k] += inr1[i * 3 + j] * inr_temp[j * 3 + k];
               }
              }
             }
            for (i=0;i<9;i++) 
            symmetry->inr[9 * l + i] /= inr_fac;
           }
*/

      p_irr1 = irr1;
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
          crystal->primitive_cell[j].comp1 += *p_irr1 * crystal->conventional_cell[i].comp1;
          crystal->primitive_cell[j].comp2 += *p_irr1 * crystal->conventional_cell[i].comp2;
          crystal->primitive_cell[j].comp3 += *p_irr1 * crystal->conventional_cell[i].comp3;
          // change order 12/05/2011  see p78 in IT vol A
          //crystal->primitive_cell[i].comp1 += *p_irr1 * crystal->conventional_cell[j].comp1;
          //crystal->primitive_cell[i].comp2 += *p_irr1 * crystal->conventional_cell[j].comp2;
          //crystal->primitive_cell[i].comp3 += *p_irr1 * crystal->conventional_cell[j].comp3;
          p_irr1++;
        }
      }

      double_vec_cross(&crystal->primitive_cell[1], &crystal->primitive_cell[2], &crystal->reciprocal_cell[0]);
      double_vec_cross(&crystal->primitive_cell[2], &crystal->primitive_cell[0], &crystal->reciprocal_cell[1]);
      double_vec_cross(&crystal->primitive_cell[0], &crystal->primitive_cell[1], &crystal->reciprocal_cell[2]);

      crystal->primitive_cell_volume = fabs(double_vec_dot(&crystal->reciprocal_cell[0], &crystal->primitive_cell[0]));

      crystal->reciprocal_cell[0].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[0].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[0].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[2].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[2].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[2].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;

    if (job->taskid == 0) {
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| LATTICE CONSTANTS         |            a = %8.5lf |            b = %8.5lf |            c = %8.5lf |\n", \
        crystal->lattice_a, crystal->lattice_b, crystal->lattice_c);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| ANGLES                    |       alpha = %9.5lf |        beta = %9.5lf |       gamma = %9.5lf |\n", \
        crystal->alpha / deg_rad, crystal->beta / deg_rad, crystal->gamma / deg_rad);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"|                                                    CONVENTIONAL UNIT CELL VECTORS (ANGS)                |\n");
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->conventional_cell[0].comp1, crystal->conventional_cell[0].comp2,crystal->conventional_cell[0].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->conventional_cell[1].comp1, crystal->conventional_cell[1].comp2,crystal->conventional_cell[1].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->conventional_cell[2].comp1, crystal->conventional_cell[2].comp2,crystal->conventional_cell[2].comp3);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"| PRIMITIVE CELL VOL. (ANGS)^3 %12.7lf            PRIMITIVE UNIT CELL VECTORS (ANGS)                 |\n",
        crystal->primitive_cell_volume);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[0].comp1, crystal->primitive_cell[0].comp2,crystal->primitive_cell[0].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[1].comp1, crystal->primitive_cell[1].comp2,crystal->primitive_cell[1].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[2].comp1, crystal->primitive_cell[2].comp2,crystal->primitive_cell[2].comp3);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"|                                                   RECIPROCAL UNIT CELL VECTORS (ANGS)^(-1)              |\n");
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[0].comp1, crystal->reciprocal_cell[0].comp2,crystal->reciprocal_cell[0].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[1].comp1, crystal->reciprocal_cell[1].comp2,crystal->reciprocal_cell[1].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[2].comp1, crystal->reciprocal_cell[2].comp2,crystal->reciprocal_cell[2].comp3);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

     }

      crystal->primitive_cell_volume /= (bohr_to_AA * bohr_to_AA * bohr_to_AA);

      break;

    case 'S':

      // 2-D periodicity

      p_irr1 = irr1;
      printf("Check trans to primitive cell\n");
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
          crystal->primitive_cell[j].comp1 += *p_irr1 * crystal->conventional_cell[i].comp1;
          crystal->primitive_cell[j].comp2 += *p_irr1 * crystal->conventional_cell[i].comp2;
          crystal->primitive_cell[j].comp3 = k_zero;
          // change order 12/05/2011  see p78 in IT vol A
          //crystal->primitive_cell[i].comp1 += *p_irr1 * crystal->conventional_cell[j].comp1;
          //crystal->primitive_cell[i].comp2 += *p_irr1 * crystal->conventional_cell[j].comp2;
          p_irr1++;
        }
      }

      det = (crystal->primitive_cell[0].comp1 * crystal->primitive_cell[1].comp2 - crystal->primitive_cell[0].comp2 \
      * crystal->primitive_cell[1].comp1);

      crystal->primitive_cell_volume = fabs(det);

      if (fabs(det) < 0.001) {
        if (job->taskid == 0)
        fprintf(file.out, "Primitive lattice vectors linearly dependent\n");
        exit(1);
      }

      crystal->reciprocal_cell[0].comp1 =  crystal->primitive_cell[1].comp2 / det * two * pi;
      crystal->reciprocal_cell[0].comp2 = -crystal->primitive_cell[1].comp1 / det * two * pi;
      crystal->reciprocal_cell[0].comp3 = k_zero;
      crystal->reciprocal_cell[1].comp1 = -crystal->primitive_cell[0].comp2 / det * two * pi;
      crystal->reciprocal_cell[1].comp2 =  crystal->primitive_cell[0].comp1 / det * two * pi;
      crystal->reciprocal_cell[1].comp3 = k_zero;
      crystal->reciprocal_cell[2].comp1 = k_zero;
      crystal->reciprocal_cell[2].comp2 = k_zero;
      crystal->reciprocal_cell[2].comp3 = k_zero;
  
/*
    if (job->taskid == 0) {
      fprintf(file.out, "LATTICE CONSTANTS a     = %10.5lf, b    = %10.5lf\n", crystal->lattice_a, crystal->lattice_b);
      fprintf(file.out, "                  gamma = %10.5lf\n", crystal->gamma / deg_rad);
      fprintf(file.out, "RECIPROCAL LATTICE CELL(ANGS^(-1))\n");
      fprintf(file.out, " %12.7lf %12.7lf %12.7lf\n", crystal->reciprocal_cell[0].comp1, crystal->reciprocal_cell[0].comp2,
          crystal->reciprocal_cell[0].comp3);
      fprintf(file.out, " %12.7lf %12.7lf %12.7lf\n", crystal->reciprocal_cell[1].comp1, crystal->reciprocal_cell[1].comp2,
          crystal->reciprocal_cell[1].comp3);
      fprintf(file.out, "PRIMITIVE UNIT CELL(ANGS)\n");
      fprintf(file.out, " %12.7lf %12.7lf %12.7lf\n", crystal->primitive_cell[0].comp1, crystal->primitive_cell[0].comp2,
          crystal->primitive_cell[0].comp3);
      fprintf(file.out, " %12.7lf %12.7lf %12.7lf\n", crystal->primitive_cell[1].comp1, crystal->primitive_cell[1].comp2,
          crystal->primitive_cell[1].comp3);

      fprintf(file.out, "PRIMITIVE UNIT CELL VOLUME (ANGS)^2 %10.4lf\n", crystal->primitive_cell_volume);
     }

      crystal->primitive_cell_volume /= (bohr_to_AA * bohr_to_AA);
*/

      // 3-D periodicity

/*
      p_irr1 = irr1;
      printf("Check trans to primitive cell\n");
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
          crystal->primitive_cell[j].comp1 += *p_irr1 * crystal->conventional_cell[i].comp1;
          crystal->primitive_cell[j].comp2 += *p_irr1 * crystal->conventional_cell[i].comp2;
          // change order 12/05/2011  see p78 in IT vol A
          //crystal->primitive_cell[i].comp1 += *p_irr1 * crystal->conventional_cell[j].comp1;
          //crystal->primitive_cell[i].comp2 += *p_irr1 * crystal->conventional_cell[j].comp2;
          p_irr1++;
        }
      }

      crystal->primitive_cell[2].comp3 = slab_real_space_z_period;

      double_vec_cross(&crystal->primitive_cell[1], &crystal->primitive_cell[2], &crystal->reciprocal_cell[0]);
      double_vec_cross(&crystal->primitive_cell[2], &crystal->primitive_cell[0], &crystal->reciprocal_cell[1]);
      double_vec_cross(&crystal->primitive_cell[0], &crystal->primitive_cell[1], &crystal->reciprocal_cell[2]);

      crystal->primitive_cell_volume = fabs(double_vec_dot(&crystal->reciprocal_cell[0], &crystal->primitive_cell[0]));

      crystal->reciprocal_cell[0].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[0].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[0].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[2].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[2].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[2].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;

*/

    if (job->taskid == 0) {
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| LATTICE CONSTANTS         |            a = %8.5lf |            b = %8.5lf |                         |\n", \
        crystal->lattice_a, crystal->lattice_b);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| ANGLE                     |       gamma = %9.5lf |                         |                         |\n", \
        crystal->gamma / deg_rad);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"|                                                    CONVENTIONAL UNIT CELL VECTORS (ANGS)                |\n");
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->conventional_cell[0].comp1, crystal->conventional_cell[0].comp2,crystal->conventional_cell[0].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->conventional_cell[1].comp1, crystal->conventional_cell[1].comp2,crystal->conventional_cell[1].comp3);
        if (job->verbosity > 1) {
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->conventional_cell[2].comp1, crystal->conventional_cell[2].comp2,crystal->conventional_cell[2].comp3);
       }
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"| PRIMITIVE CELL VOL. (ANGS)^3 %12.7lf            PRIMITIVE UNIT CELL VECTORS (ANGS)                 |\n",
        crystal->primitive_cell_volume / slab_real_space_z_period);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[0].comp1, crystal->primitive_cell[0].comp2,crystal->primitive_cell[0].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[1].comp1, crystal->primitive_cell[1].comp2,crystal->primitive_cell[1].comp3);
        if (job->verbosity > 1) {
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[2].comp1, crystal->primitive_cell[2].comp2,crystal->primitive_cell[2].comp3);
       }
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"|                                                   RECIPROCAL UNIT CELL VECTORS (ANGS)^(-1)              |\n");
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[0].comp1, crystal->reciprocal_cell[0].comp2,crystal->reciprocal_cell[0].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[1].comp1, crystal->reciprocal_cell[1].comp2,crystal->reciprocal_cell[1].comp3);
        if (job->verbosity > 1) {
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[2].comp1, crystal->reciprocal_cell[2].comp2,crystal->reciprocal_cell[2].comp3);
       }
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

     }

      crystal->primitive_cell_volume /= (bohr_to_AA * bohr_to_AA);

      break;

    case 'P':

      for (i = 0; i < 3; i++) {
        crystal->primitive_cell[i].comp1 = k_zero;
        crystal->primitive_cell[i].comp2 = k_zero;
        crystal->primitive_cell[i].comp3 = k_zero;
        //crystal->primitive_cell[i] = crystal->conventional_cell[i];
        crystal->reciprocal_cell[i].comp1 = k_zero;
        crystal->reciprocal_cell[i].comp2 = k_zero;
        crystal->reciprocal_cell[i].comp3 = k_zero;
      }

/*
      // use only for true 1-D periodicity

      crystal->primitive_cell_volume = crystal->lattice_c;

      crystal->reciprocal_cell[0].comp1 = k_zero;
      crystal->reciprocal_cell[0].comp2 = k_zero;
      //crystal->reciprocal_cell[0].comp3 = k_zero;
      crystal->reciprocal_cell[0].comp3 = two * pi / crystal->lattice_c;
      crystal->reciprocal_cell[1].comp1 = k_zero;
      crystal->reciprocal_cell[1].comp2 = k_zero;
      crystal->reciprocal_cell[1].comp3 = k_zero;
      crystal->reciprocal_cell[2].comp1 = k_zero;
      crystal->reciprocal_cell[2].comp2 = k_zero;
      crystal->reciprocal_cell[2].comp3 = k_zero;
      //crystal->reciprocal_cell[2].comp3 = 2.0 * pi / crystal->lattice_c;

      crystal->primitive_cell_volume /= bohr_to_AA;

*/

      crystal->primitive_cell[0].comp3 = crystal->lattice_c;
      crystal->primitive_cell[1].comp1 = rod_real_space_x_y_period;
      crystal->primitive_cell[2].comp2 = rod_real_space_x_y_period;
      crystal->primitive_cell[1].comp1 = 20.0;
      crystal->primitive_cell[2].comp2 = 20.0;

      //double_vec_cross(&crystal->primitive_cell[1], &crystal->primitive_cell[2], &crystal->reciprocal_cell[0]);
      //double_vec_cross(&crystal->primitive_cell[2], &crystal->primitive_cell[0], &crystal->reciprocal_cell[1]);
      //double_vec_cross(&crystal->primitive_cell[0], &crystal->primitive_cell[1], &crystal->reciprocal_cell[2]);

      //crystal->primitive_cell_volume = fabs(double_vec_dot(&crystal->reciprocal_cell[0], &crystal->primitive_cell[0]));
      crystal->primitive_cell_volume = crystal->primitive_cell[0].comp3 * crystal->primitive_cell[1].comp1 * crystal->primitive_cell[2].comp2;

      crystal->reciprocal_cell[0].comp3 = two * pi / crystal->lattice_c;
      //crystal->reciprocal_cell[0].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      //crystal->reciprocal_cell[0].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      //crystal->reciprocal_cell[0].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;
      crystal->reciprocal_cell[1].comp1 = two * pi / rod_real_space_x_y_period;
      crystal->reciprocal_cell[2].comp2 = two * pi / rod_real_space_x_y_period;
      crystal->reciprocal_cell[1].comp1 = two * pi / 20.0;
      crystal->reciprocal_cell[2].comp2 = two * pi / 20.0;
      //crystal->reciprocal_cell[1].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      //crystal->reciprocal_cell[1].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;
      //crystal->reciprocal_cell[2].comp1 *= 2.0 * pi / crystal->primitive_cell_volume;
      //crystal->reciprocal_cell[2].comp2 *= 2.0 * pi / crystal->primitive_cell_volume;
      //crystal->reciprocal_cell[2].comp3 *= 2.0 * pi / crystal->primitive_cell_volume;

    if (job->taskid == 0) {
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"| LATTICE CONSTANT          |            a = %8.5lf |                         |                         |\n", \
        crystal->lattice_c);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"| PRIMITIVE CELL VOL. (ANGS) %12.7lf              PRIMITIVE UNIT CELL VECTOR (ANGS)                  |\n",
        crystal->lattice_c);
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[0].comp1, crystal->primitive_cell[0].comp2,crystal->primitive_cell[0].comp3);
        if (job->verbosity > 1) {
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[1].comp1, crystal->primitive_cell[1].comp2,crystal->primitive_cell[1].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->primitive_cell[2].comp1, crystal->primitive_cell[2].comp2,crystal->primitive_cell[2].comp3);
       }
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

        fprintf(file.out,"|                                                   RECIPROCAL UNIT CELL VECTORS (ANGS)^(-1)              |\n");
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[0].comp1, crystal->reciprocal_cell[0].comp2,crystal->reciprocal_cell[0].comp3);
        if (job->verbosity > 1) {
	fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[1].comp1, crystal->reciprocal_cell[1].comp2,crystal->reciprocal_cell[1].comp3);
        fprintf(file.out,"|                           |            %12.7lf |            %12.7lf |            %12.7lf |\n",\
        crystal->reciprocal_cell[2].comp1, crystal->reciprocal_cell[2].comp2,crystal->reciprocal_cell[2].comp3);
       }
        fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

     }

      crystal->primitive_cell_volume /= bohr_to_AA;
      //crystal->primitive_cell_volume /= (bohr_to_AA * bohr_to_AA * bohr_to_AA);

      break;

    case 'M':

      break;

  } // close switch{crystal->type[0]

  for (i = 0; i < 3; i++) {
    crystal->primitive_cell[i].comp1 /= bohr_to_AA;
    crystal->primitive_cell[i].comp2 /= bohr_to_AA;
    crystal->primitive_cell[i].comp3 /= bohr_to_AA;
    crystal->conventional_cell[i].comp1 /= bohr_to_AA;
    crystal->conventional_cell[i].comp2 /= bohr_to_AA;
    crystal->conventional_cell[i].comp3 /= bohr_to_AA;
    crystal->reciprocal_cell[i].comp1 *= bohr_to_AA;
    crystal->reciprocal_cell[i].comp2 *= bohr_to_AA;
    crystal->reciprocal_cell[i].comp3 *= bohr_to_AA;
  }

}

