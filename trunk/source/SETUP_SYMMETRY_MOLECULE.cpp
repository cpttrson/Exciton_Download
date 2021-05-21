  /*! \mainpage Exciton08 Manual
 this is what appears on the main page of the Exciton08 manual
 */

/*! \file MAIN.cpp
 \brief main file
 \details a proper description
 should be added later
 \remarks  here just remarks

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include "mycomplex.h"
#include "LIMITS.h"
#include "CRYSTAL09.h"
 */
#include <mpi.h>
#include "mylogical.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "PRINT_UTIL.h"
#include "ALLOCATE_MEMORY.h"
#include "TOOLS.h"
#include "SETUP_SYMMETRY.h"

using namespace std;

/*!
 \brief main function
 \details controls the flow of the programe and does the reading of the input files
 */

void read_symmetry_group(CRYSTAL *crystal, SYMMETRY *symmetry, int print, JOB_PARAM *job, FILES file)

{

  int i, j, k, n;
  int space_group, layer_group, rod_group, point_group;
  int group;
  int number_of_latt_vectors;
  int number_of_generators;
  int line_number;
  char job_title[80];
  char title[150];
  char junk, t1, t2, t3, u1, u2, u3, v1, v2, v3;
  char t[4], u[4], v[4], r[4], PQR[9];
  const char *explicit_symbol, *Hermann_Mauguin_symbol, *Schoenflies_symbol;
  char Schoenflies[9], Hermann_Mauguin[18];
  double irr1[9], gen_irr[9 * 72], *p_gen_irr;
  double angle;
  VECTOR_DOUBLE tau[8];
  VECTOR_DOUBLE f_tmp;

  /******************************************************************************************
   *  Read data from input file 'INPUT'                                                     *
   ******************************************************************************************/

  line_number = 0;

  read_line(file.in, title, 150);
  line_number++;
  sscanf(title, "%s", job_title);
  //if (job->taskid == 0 && print == 1)
  //fprintf(file.out, "Job Description       %s\n", job_title);
  read_line(file.in, title, 150);
  line_number++;
  sscanf(title, "%s", crystal->type);

  /***********************************************************************************************************************
   *  Parse explicit symbol for this space group. See Table A1.4.2.1 in Vol. B of International Tables for original data *
   **********************************************************************************************************************/

  switch (crystal->type[0]) {

    case 'C':

      break;

    case 'S': // Layer groups

      break;

    case 'P': // Rod groups

      break;

    case 'M': // Point groups

      read_line(file.in, title, 8);
      line_number++;
      sscanf(title, "%d", &group);

      number_of_latt_vectors = 0;

  /***********************************************************************************************************************
   *  Explicit symbols for point groups obtained from first incidence of that group in Table A1.4.2.1 in Vol. B of       *
   *  International Tables.                                                                                              *
   **********************************************************************************************************************/

      if (group >=  1 && group <=  8 || group == 17 || group == 18 || group == 25 || group == 26 || group == 33 || group == 34)
        number_of_generators = 1;
      if (group >=  9 && group <= 15 || group >= 19 && group <= 23 || group >= 27 && group <= 32 || group >= 35 || group == 39)
        number_of_generators = 2;
      if (group == 16 || group == 24 || group >= 40 && group <= 45)
        number_of_generators = 3;

      switch (group) {

        case 1:
          explicit_symbol = "PAN$P1A000";
          Schoenflies_symbol = "C_1";
          sprintf(Schoenflies,"C_1");
          break; 
        case 2:
          explicit_symbol = "PAC$I1A000";
          Schoenflies_symbol = "C_i";
          sprintf(Schoenflies,"C_i");
          break; 
        case 3:
          explicit_symbol = "PAN$P2A000";
          Schoenflies_symbol = "C_2 (x)";
          sprintf(Schoenflies,"C_2 (x)");
          break; 
        case 4:
          explicit_symbol = "PAN$P2B000";
          Schoenflies_symbol = "C_2 (y)";
          sprintf(Schoenflies,"C_2 (y)");
          break; 
        case 5:
          explicit_symbol = "PAN$P2C000";
          Schoenflies_symbol = "C_2 (z)";
          sprintf(Schoenflies,"C_2 (z)");
          break; 
        case 6:
          explicit_symbol = "PMN$I2A000";
          Schoenflies_symbol = "C_s (x)";
          sprintf(Schoenflies,"C_s (x)");
          break;
        case 7:
          explicit_symbol = "PAN$I2B000";
          Schoenflies_symbol = "C_s (y)";
          sprintf(Schoenflies,"C_s (y)");
          break; 
        case 8:
          explicit_symbol = "PAN$I2C000";
          Schoenflies_symbol = "C_s (z)";
          sprintf(Schoenflies,"C_s (z)");
          break;
        case 9:
          //explicit_symbol = "PAC$P2A000";
          explicit_symbol = "PAC$I1A000$P2A000";
          Schoenflies_symbol = "C_2h (x)";
          sprintf(Schoenflies,"C_2h (x)");
          break;
        case 10:
          //explicit_symbol = "PAC$P2B000";
          explicit_symbol = "PAC$I1A000$P2B000";
          Schoenflies_symbol = "C_2h (y)";
          sprintf(Schoenflies,"C_2h (y)");
          break;
        case 11:
          //explicit_symbol = "PAC$P2C000";
          explicit_symbol = "PAC$I1A000$P2C000";
          Schoenflies_symbol = "C_2h (z)";
          sprintf(Schoenflies,"C_2h (z)");
          break;
        case 12:
          explicit_symbol = "PAN$P2B000$P2C000";
          Schoenflies_symbol = "D_2";
          sprintf(Schoenflies,"D_2");
          break;
        case 13:
          explicit_symbol = "PAN$P2A000$I2B000";
          Schoenflies_symbol = "C_2v (x)";
          sprintf(Schoenflies,"C_2v (x)");
          break;
        case 14:
          explicit_symbol = "PAN$P2B000$I2C000";
          Schoenflies_symbol = "C_2v (y)";
          sprintf(Schoenflies,"C_2v (y)");
          break;
        case 15:
          explicit_symbol = "PAN$P2C000$I2A000";
          Schoenflies_symbol = "C_2v (z)";
          sprintf(Schoenflies,"C_2v (z)");
          break;
        case 16:
          //explicit_symbol = "PAC$P2B000$P2C000";
          explicit_symbol = "PAC$I1A000$P2C000$P2A000";
          Schoenflies_symbol = "D_2h";
          sprintf(Schoenflies,"D_2h");
          break;
        case 17:
          //explicit_symbol = "PAN$P2C000$P4C000";
          explicit_symbol = "PAN$P4C000";
          Schoenflies_symbol = "C_4";
          sprintf(Schoenflies,"C_4");
          break;
        case 18:
          //explicit_symbol = "PAN$P2C000$I4C000";
          explicit_symbol = "PAN$I4C000";
          Schoenflies_symbol = "S_4";
          sprintf(Schoenflies,"S_4");
          break;
        case 19:
          //explicit_symbol = "PAC$P2C000$P4C000";
          explicit_symbol = "PAC$I1A000$P4C000";
          Schoenflies_symbol = "C_4h";
          break;
        case 20:
          //explicit_symbol = "PAN$P2B000$P2C000$P4C000";
          explicit_symbol = "PAN$P4C000$P2A000";
          Schoenflies_symbol = "D_4";
          sprintf(Schoenflies,"D_4");
          break;
        case 21:
          //explicit_symbol = "PAN$I2B000$P2C000$P4C000";
          explicit_symbol = "PAN$P4C000$I2A000";
          Schoenflies_symbol = "C_4v";
          sprintf(Schoenflies,"C_4v");
          break;
        case 22:
          //explicit_symbol = "PAN$P2B000$P2C000$I4C000";
          explicit_symbol = "PAN$I4C000$P2A000";
          Schoenflies_symbol = "D_2d";
          sprintf(Schoenflies,"D_2d");
          break;
        case 23:
          //explicit_symbol = "PAN$P2B000$P2C000$I4C000";
          explicit_symbol = "PAN$I4C000$P2A000";
          Schoenflies_symbol = "D_2d";
          sprintf(Schoenflies,"D_2d");
          break;
        case 24:
          //explicit_symbol = "PAC$P2B000$P2C000$P4C000";
          explicit_symbol = "PAC$I1A000$P4C000$P2A000";
          Schoenflies_symbol = "D_4h";
          sprintf(Schoenflies,"D_4h");
          break;
        case 25:
          //explicit_symbol = "PAN$P3C000";
          explicit_symbol = "PHN$P3C000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_3";
          sprintf(Schoenflies,"C_3");
          break;
        case 26:
          //explicit_symbol = "PHC$I3C000";
          explicit_symbol = "PHC$I3C000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_3i";
          sprintf(Schoenflies,"C_3i");
          break;
        case 27:
          //explicit_symbol = "PAN$P3C000$P2E000";
          explicit_symbol = "PHN$P3C000$P2G000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_3";
          sprintf(Schoenflies,"D_3");
          break;
        case 28:
          //explicit_symbol = "PHN$P3C000$P2E000";
          explicit_symbol = "PHN$P3C000$P2G000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_3";
          sprintf(Schoenflies,"D_3");
          break;
        case 29:
          //explicit_symbol = "PHN$P3C000$I2D000";
          explicit_symbol = "PHN$P3C000$I2F000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_3v";
          sprintf(Schoenflies,"C_3v");
          break;
        case 30:
          //explicit_symbol = "PHN$P3C000$I2D000";
          explicit_symbol = "PHN$P3C000$I2F000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_3v";
          sprintf(Schoenflies,"C_3v");
          break;
        case 31:
          //explicit_symbol = "PHC$P3C000$P2E000";
          explicit_symbol = "PHC$I3C000$P2G000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_3d";
          sprintf(Schoenflies,"D_3d");
          break;
        case 32:
          //explicit_symbol = "PHC$P3C000$P2E000";
          explicit_symbol = "PHC$I3C000$P2G000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_3d";
          sprintf(Schoenflies,"D_3d");
          break;
        case 33:
          //explicit_symbol = "PHN$P2C000$P3C000";
          explicit_symbol = "PHN$P6C000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_6";
          sprintf(Schoenflies,"C_6");
          break;
        case 34:
          //explicit_symbol = "PHN$I2C000$P3C000";
          explicit_symbol = "PHN$I6C000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_3h";
          sprintf(Schoenflies,"C_3h");
          break;
        case 35:
          //explicit_symbol = "PHC$P2C000$P3C000";
          explicit_symbol = "PHC$I1A000$P6C000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_6h";
          sprintf(Schoenflies,"C_6h");
          break;
        case 36:
          //explicit_symbol = "PHN$P2C000$P3C000$P2D000";
          explicit_symbol = "PHN$P6C000$P2F000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_6";
          sprintf(Schoenflies,"D_6");
          break;
        case 37:
          //explicit_symbol = "PHN$P2C000$P3C000$I2D000";
          explicit_symbol = "PHN$P6C000$I2F000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "C_6v";
          sprintf(Schoenflies,"C_6v");
          break;
        case 38:
          //explicit_symbol = "PHN$I2C000$P3C000$I2D000";
          explicit_symbol = "PHN$I6C000$P2G000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_3h";
          sprintf(Schoenflies,"D_3h");
          break;
        case 39:
          //explicit_symbol = "PHN$I2C000$P3C000$I2D000";
          explicit_symbol = "PHN$I6C000$P2G000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_3h needs attention";
          sprintf(Schoenflies,"D_3h");
          break;
        case 40:
          //explicit_symbol = "PAC$P2C000$P3C000$P2D000";
          explicit_symbol = "PHC$I1A000$P6C000$P2F000";
          crystal->gamma = 120.0 * deg_rad;
          Schoenflies_symbol = "D_6h";
          sprintf(Schoenflies,"D_6h");
          break;
        case 41:
          //explicit_symbol = "PAN$P2B000$P2C000$P3Q000";
          explicit_symbol = "PAN$P3Q000$P2C000$P2A000";
          Schoenflies_symbol = "T";
          sprintf(Schoenflies,"T");
          break;
        case 42:
          //explicit_symbol = "PAC$P2B000$P2C000$P3Q000";
          explicit_symbol = "PAC$I3Q000$P2C000$P2A000";
          Schoenflies_symbol = "T_h";
          sprintf(Schoenflies,"T_h");
          break;
        case 43:
          //explicit_symbol = "PAN$P3Q000$P4C000$P2D000";
          explicit_symbol = "PAN$P3Q000$P4C000$P2D000";
          Schoenflies_symbol = "O";
          sprintf(Schoenflies,"O");
          break;
        case 44:
          explicit_symbol = "PAN$P3Q000$I4C000$I2D000";
          Schoenflies_symbol = "T_d";
          sprintf(Schoenflies,"T_d");
          break;
        case 45:
          explicit_symbol = "PAC$I3Q000$P4C000$P2D000";
          Schoenflies_symbol = "O_h";
          sprintf(Schoenflies,"O_h");
          break;

          break;

      } // close switch(group)

      if (group < 1 || group > 45) {
        if (job->taskid == 0 && print == 1)
        fprintf(file.out,"ERROR: Molecular point group must be in range 1 to 45. Value entered = %5d\n",group);
        MPI_Finalize();
        exit(1); }

      //fprintf(file.out,"%s\n",explicit_symbol) ;
      if (number_of_generators == 1) {
        sscanf(explicit_symbol, "%c%c%c%c%c%2s%c%c%c", &crystal->centring, &crystal->system, &crystal->centrosymm, &junk,\
            &r[0],&PQR[0], &t[0], &u[0], &v[0]);
        fprintf(file.out,"%c%c%c%c%c%2s%c%c%c\n\n",crystal->centring,crystal->system,crystal->centrosymm,junk,r[0],PQR,t[0],u[0],v[0]) ;
      }

      if (number_of_generators == 2) {
        sscanf(explicit_symbol, "%c%c%c%c%c%2s%c%c%c%c%c%2s%c%c%c", &crystal->centring, &crystal->system, &crystal->centrosymm,
            &junk, &r[0], &PQR[0], &t[0], &u[0], &v[0], &junk, &r[1], &PQR[2], &t[1], &u[1], &v[1]);
        //fprintf(file.out,"%c%c%c%c%c%2s%c%c%c%c%c%2s%c%c%c\n\n",crystal->centring,crystal->system,crystal->centrosymm,\
        junk,r[0],PQR,t[0],u[0],v[0],junk,r[1],PQR,t[1],u[1],v[1]) ;
      }

      if (number_of_generators == 3) {
        sscanf(explicit_symbol, "%c%c%c%c%c%2s%c%c%c%c%c%2s%c%c%c%c%c%2s%c%c%c", &crystal->centring, &crystal->system,
            &crystal->centrosymm, &junk, &r[0], &PQR[0], &t[0], &u[0], &v[0], &junk, &r[1], &PQR[2], &t[1], &u[1], &v[1],
            &junk, &r[2], &PQR[4], &t[2], &u[2], &v[2]);
        //fprintf(file.out,"%c%c%c%c%c%2s%c%c%c%c%c%2s%c%c%c%c%c%2s%c%c%c\n\n",crystal->centring,crystal->system,crystal->centrosymm,    \
        junk,r[0],PQR,t[0],u[0],v[0],junk,r[1],PQR,t[1],u[1],v[1],junk,r[2],PQR,t[2],u[2],v[2]) ;
      }

  } // close switch(crystal->type)

  /******************************************************************************************
   *  Generate full set of symmetry operators from generators                               *
   ******************************************************************************************/

  p_gen_irr = gen_irr;
  double dum, *p_dum;
  p_dum = &dum;
  symmetry->number_of_operators = 1;
  symmetry->irr[0] = k_one;
  symmetry->irr[1] = k_zero;
  symmetry->irr[2] = k_zero;
  symmetry->irr[3] = k_zero;
  symmetry->irr[4] = k_one;
  symmetry->irr[5] = k_zero;
  symmetry->irr[6] = k_zero;
  symmetry->irr[7] = k_zero;
  symmetry->irr[8] = k_one;
  symmetry->inr[0] = 1;
  symmetry->inr[1] = 0;
  symmetry->inr[2] = 0;
  symmetry->inr[3] = 0;
  symmetry->inr[4] = 1;
  symmetry->inr[5] = 0;
  symmetry->inr[6] = 0;
  symmetry->inr[7] = 0;
  symmetry->inr[8] = 1;
  symmetry->taur[0].comp1 = k_zero;
  symmetry->taur[0].comp2 = k_zero;
  symmetry->taur[0].comp3 = k_zero;

  for (i = 0; i < number_of_generators; i++) {

    switch (t[i]) {
      case '0':
        tau[i].comp1 = k_zero;
        break;
      case '2':
        tau[i].comp1 = sixth;
        break;
      case '3':
        tau[i].comp1 = quarter;
        break;
      case '4':
        tau[i].comp1 = third;
        break;
      case '5':
        tau[i].comp1 = five_sixths;
        break;
      case '6':
        tau[i].comp1 = half;
        break;
      case '8':
        tau[i].comp1 = two_thirds;
        break;
      case '9':
        tau[i].comp1 = three_quarters;
        break;
    }
    switch (u[i]) {
      case '0':
        tau[i].comp2 = k_zero;
        break;
      case '2':
        tau[i].comp2 = sixth;
        break;
      case '3':
        tau[i].comp2 = quarter;
        break;
      case '4':
        tau[i].comp2 = third;
        break;
      case '5':
        tau[i].comp2 = five_sixths;
        break;
      case '6':
        tau[i].comp2 = half;
        break;
      case '8':
        tau[i].comp2 = two_thirds;
        break;
      case '9':
        tau[i].comp2 = three_quarters;
        break;
    }
    switch (v[i]) {
      case '0':
        tau[i].comp3 = k_zero;
        break;
      case '2':
        tau[i].comp3 = sixth;
        break;
      case '3':
        tau[i].comp3 = quarter;
        break;
      case '4':
        tau[i].comp3 = third;
        break;
      case '5':
        tau[i].comp3 = five_sixths;
        break;
      case '6':
        tau[i].comp3 = half;
        break;
      case '8':
        tau[i].comp3 = two_thirds;
        break;
      case '9':
        tau[i].comp3 = three_quarters;
        break;
    }

    switch (PQR[2 * i]) {

      case '1': //fprintf(file.out,"its 1A ") ;

        // operator 1A identity

        *p_gen_irr = k_one;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_one;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_one;
        p_gen_irr++;

        break;

      case '2': //fprintf(file.out,"its 2") ;

        switch (PQR[2 * i + 1]) {

          case 'A': //fprintf(file.out,"A ") ;

            // operator 2A 2 x,0,0

            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;

            break;

          case 'B': //fprintf(file.out,"B ") ;

            // operator 2B 2 0,y,0

            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;

            break;

          case 'C': //fprintf(file.out,"C ") ;

            // operator 2C 2 0,0,z

            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;

            break;

          case 'D': //fprintf(file.out,"D ") ;

            // operator 2D 2 x,x,0

            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;

            break;

          case 'E': //fprintf(file.out,"E ") ;

            // operator 2E 2 x,-x,0

            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;

            break;

          case 'F': //fprintf(file.out,"F ") ;

            // operator 2F 2 x,0,0

            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;

            break;

          case 'G': //fprintf(file.out,"G ") ;

            // operator 2G 2 2x,0,0

            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;

            break;

        }
        break;

      case '3': //fprintf(file.out,"its 3") ;

        switch (PQR[2 * i + 1]) {

          case 'C': //fprintf(file.out,"C ") ;

            // operator 3C -6- 0,0,z

            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = -k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            //*p_gen_irr = -k_one; 17/04/12
            p_gen_irr++;

            break;

          case 'Q': //fprintf(file.out,"Q ") ;

            // operator 3Q 3+ x,x,x

            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;
            *p_gen_irr = k_one;
            p_gen_irr++;
            *p_gen_irr = k_zero;
            p_gen_irr++;

            break;

        }
        break;

      case '4': //fprintf(file.out,"its 4C ") ;

        // operator 4C 4+ 0,0,z

        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = -k_one;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_one;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_one;
        p_gen_irr++;

        break;

      case '6': //fprintf(file.out,"its 6C ") ;

        // operator 6C 6+ 0,0,z

        //*p_gen_irr =  k_zero ; p_gen_irr++ ;
        *p_gen_irr = k_one;
        p_gen_irr++;
        *p_gen_irr = -k_one;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_one;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_zero;
        p_gen_irr++;
        *p_gen_irr = k_one;
        p_gen_irr++;

        break;

    }

    p_gen_irr = gen_irr;
    if (r[i] == 'I') {
          for (n = 0; n < 9; n++) {
            *p_gen_irr *= -k_one;
            p_gen_irr++;
          }
         }

    p_gen_irr = gen_irr;
    for (j = 0; j < symmetry->number_of_operators; j++) {
      double_mat_dot(p_gen_irr, symmetry->irr + 9 * j, irr1);
      rotate_vector3(p_gen_irr, &symmetry->taur[j], &f_tmp);
      f_tmp.comp1 += tau[i].comp1;
      f_tmp.comp1 = modf(f_tmp.comp1, p_dum);
      if (f_tmp.comp1 < k_zero)
        f_tmp.comp1 += k_one;
      f_tmp.comp2 += tau[i].comp2;
      f_tmp.comp2 = modf(f_tmp.comp2, p_dum);
      if (f_tmp.comp2 < k_zero)
        f_tmp.comp2 += k_one;
      f_tmp.comp3 += tau[i].comp3;
      f_tmp.comp3 = modf(f_tmp.comp3, p_dum);
      if (f_tmp.comp3 < k_zero)
        f_tmp.comp3 += k_one;
      for (k = 0; k < symmetry->number_of_operators; k++) {
        if (check_mat(&symmetry->irr[9 * k], irr1) == 1)
          break;
        if (check_mat(&symmetry->irr[9 * k], irr1) == 0 && k == symmetry->number_of_operators - 1) {

          symmetry->irr[9 * symmetry->number_of_operators + 0] = irr1[0];
          symmetry->irr[9 * symmetry->number_of_operators + 1] = irr1[1];
          symmetry->irr[9 * symmetry->number_of_operators + 2] = irr1[2];
          symmetry->irr[9 * symmetry->number_of_operators + 3] = irr1[3];
          symmetry->irr[9 * symmetry->number_of_operators + 4] = irr1[4];
          symmetry->irr[9 * symmetry->number_of_operators + 5] = irr1[5];
          symmetry->irr[9 * symmetry->number_of_operators + 6] = irr1[6];
          symmetry->irr[9 * symmetry->number_of_operators + 7] = irr1[7];
          symmetry->irr[9 * symmetry->number_of_operators + 8] = irr1[8];

          symmetry->inr[9 * symmetry->number_of_operators + 0] = (int)irr1[0];
          symmetry->inr[9 * symmetry->number_of_operators + 1] = (int)irr1[1];
          symmetry->inr[9 * symmetry->number_of_operators + 2] = (int)irr1[2];
          symmetry->inr[9 * symmetry->number_of_operators + 3] = (int)irr1[3];
          symmetry->inr[9 * symmetry->number_of_operators + 4] = (int)irr1[4];
          symmetry->inr[9 * symmetry->number_of_operators + 5] = (int)irr1[5];
          symmetry->inr[9 * symmetry->number_of_operators + 6] = (int)irr1[6];
          symmetry->inr[9 * symmetry->number_of_operators + 7] = (int)irr1[7];
          symmetry->inr[9 * symmetry->number_of_operators + 8] = (int)irr1[8];

       // Note that symmetry operators in integer representation are transformed to primitive cell basis in void primitive_unit_cell() below

          symmetry->taur[symmetry->number_of_operators].comp1 = f_tmp.comp1;
          symmetry->taur[symmetry->number_of_operators].comp2 = f_tmp.comp2;
          symmetry->taur[symmetry->number_of_operators].comp3 = f_tmp.comp3;
          symmetry->number_of_operators++;

        } // close if (check_mat ...
      } // close loop on k
    } // close loop on j
  } // close loop on i

      if (job->taskid == 0 && print == 1) {
  
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| JOB DESCRIPTION %-79s         |\n", job_title);
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");

      if (crystal->type[0] == 'C') {
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| SPACE GROUP           %3d | HERMANN MAUGUIN %7s | ORIGIN              %3d | Cell                %3d |\n", \
      space_group,Hermann_Mauguin_symbol,crystal->cell_choice, crystal->origin);
     }
      else if (crystal->type[0] == 'S') {
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| LAYER GROUP           %3d | HERMANN MAUGUIN  %8s |                           | Cell                    |\n", \
      //fprintf(file.out,"| LAYER GROUP           %3d | HERMANN MAUGUIN  %6s | ORIGIN              %3d | Cell                %3d |\n", 
      layer_group,Hermann_Mauguin_symbol);
     }
      else if (crystal->type[0] == 'P') {
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| ROD GROUP             %3d | HERMANN MAUGUIN  %6s | ORIGIN              %3d | Cell                %3d |\n", 
      rod_group,Hermann_Mauguin,crystal->cell_choice, crystal->origin);
      //fprintf(file.out,"| ROD GROUP             %3d | HERMANN MAUGUIN  %11s                      |                         |\n", \
      rod_group,Hermann_Mauguin,crystal->cell_choice, crystal->origin);
     }
      else if (crystal->type[0] == 'M') {
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| SYMMETRY GROUP        %3d | SYMBOL            %5s | GENERATORS          %3d | SYMMETRY OPERATORS  %3d |\n", \
      group,Schoenflies,number_of_generators,symmetry->number_of_operators);
     }
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    for (j=0;j<symmetry->number_of_operators;j++) fprintf(file.out,\
     "|                        %2d |   %3d %3d %3d   %5.2lf   |   %3d %3d %3d   %5.2lf   |   %3d %3d %3d   %5.2lf   |\n", \
      j+1,symmetry->inr[j*9+0],symmetry->inr[j*9+1],symmetry->inr[j*9+2],symmetry->taur[j].comp1, \
          symmetry->inr[j*9+3],symmetry->inr[j*9+4],symmetry->inr[j*9+5],symmetry->taur[j].comp2, \
          symmetry->inr[j*9+6],symmetry->inr[j*9+7],symmetry->inr[j*9+8],symmetry->taur[j].comp3) ;
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");
  
     }

}

void generate_cartesian_symmetry_operators(CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Convert symmetry operators to Cartesian basis                                         *
  // ******************************************************************************************

  int i, j, k ;
  double irr1[9], *p_irr1, irr1_inv[9], *p_irr1_inv;
  double *p_irr, *q_irr, irr_tmp[9], irr_tmp1[9];
  VECTOR_DOUBLE tau[48];
  VECTOR_DOUBLE Rvec_tmp[2];

  // first convert translations to primitive Cartesian basis

  p_irr1 = irr1;
  p_irr1_inv = irr1_inv;

  switch (crystal->centring) {

    case 'P': // P -> P

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
      *p_irr1_inv = half;
      p_irr1_inv++;
      *p_irr1_inv = half;
      p_irr1_inv++;
      *p_irr1_inv = k_zero;
      p_irr1_inv++;
      *p_irr1_inv = -half;
      p_irr1_inv++;
      *p_irr1_inv = half;
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

      if (job->taskid == 0)
      fprintf(file.out, "Insert correct p_irr1_inv array line 1157 in MAIN.cpp\n");
      exit(1);

      break;

    case 'C': // C -> P

      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = half;
      p_irr1++;
      *p_irr1 = k_zero;
      p_irr1++;
      *p_irr1 = -half;
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

    case 'R': // R -> H

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

      break;
}

  for (i = 0; i < symmetry->number_of_operators; i++) {
    tau[0].comp1 = symmetry->taur[i].comp1;
    tau[0].comp2 = symmetry->taur[i].comp2;
    tau[0].comp3 = symmetry->taur[i].comp3;
    //fprintf(file.out,"%d %lf %lf %lf\n",i,tau[0].comp1,tau[0].comp2,tau[0].comp3) ;
    p_irr1_inv = irr1_inv;
    tau[1].comp1 = fmod((*(p_irr1_inv + 0) * tau[0].comp1 + *(p_irr1_inv + 1) * tau[0].comp2 + *(p_irr1_inv + 2)
        * tau[0].comp3), k_one);
    tau[1].comp2 = fmod((*(p_irr1_inv + 3) * tau[0].comp1 + *(p_irr1_inv + 4) * tau[0].comp2 + *(p_irr1_inv + 5)
        * tau[0].comp3), k_one);
    tau[1].comp3 = fmod((*(p_irr1_inv + 6) * tau[0].comp1 + *(p_irr1_inv + 7) * tau[0].comp2 + *(p_irr1_inv + 8)
        * tau[0].comp3), k_one);
    //fprintf(file.out,"%d %lf %lf %lf\n",i,tau[1].comp1,tau[1].comp2,tau[1].comp3) ;
  symmetry->taur[i].comp1 = tau[1].comp1 * crystal->primitive_cell[0].comp1 + tau[1].comp2 * crystal->primitive_cell[1].comp1 +\
  tau[1].comp3 * crystal->primitive_cell[2].comp1;
  symmetry->taur[i].comp2 = tau[1].comp1 * crystal->primitive_cell[0].comp2 + tau[1].comp2 * crystal->primitive_cell[1].comp2 +\
  tau[1].comp3 * crystal->primitive_cell[2].comp2;
  symmetry->taur[i].comp3 = tau[1].comp1 * crystal->primitive_cell[0].comp3 + tau[1].comp2 * crystal->primitive_cell[1].comp3 +\
  tau[1].comp3 * crystal->primitive_cell[2].comp3;
    //fprintf(file.out,"%d %lf %lf %lf\n",i,symmetry->taur[i].comp1,symmetry->taur[i].comp2,symmetry->taur[i].comp3) ;
  }

  // then transform symmetry operators in conventional cell vector basis to Cartesian basis
  // this is only necessary for monoclinic, rhombohedral and hexagonal lattices

  for (i = 0; i < symmetry->number_of_operators; i++) {
    switch (crystal->system) {
      case 'M':
        break;
      case 'R':
// CHP 03/11/2020
        irr_tmp[0] =  rtthree / two;
        irr_tmp[1] =  k_zero;
        irr_tmp[2] =  k_zero;
        irr_tmp[3] = -half;
        irr_tmp[4] =  k_one;
        irr_tmp[5] =  k_zero;
        irr_tmp[6] =  k_zero;
        irr_tmp[7] =  k_zero;
        irr_tmp[8] =  k_one;
        double_mat_dot(irr_tmp, symmetry->irr + 9 * i, irr_tmp1);
        irr_tmp[0] =  two / rtthree;
        irr_tmp[1] =  k_zero;
        irr_tmp[2] =  k_zero;
        irr_tmp[3] =  k_one / rtthree;
        irr_tmp[4] =  k_one;
        irr_tmp[5] =  k_zero;
        irr_tmp[6] =  k_zero;
        irr_tmp[7] =  k_zero;
        irr_tmp[8] =  k_one;
        double_mat_dot(irr_tmp1, irr_tmp, symmetry->irr + 9 * i);
        break;
      case 'H': {
        irr_tmp[0] = k_one;
        irr_tmp[1] = k_zero;
        irr_tmp[2] = k_zero;
        irr_tmp[3] = cos(crystal->gamma) / sin(crystal->gamma);
        irr_tmp[4] = k_one / sin(crystal->gamma);
        irr_tmp[5] = k_zero;
        irr_tmp[6] = k_zero;
        irr_tmp[7] = k_zero;
        irr_tmp[8] = k_one;
        double_mat_dot(irr_tmp, symmetry->irr + 9 * i, irr_tmp1);
//fprintf(file.out,"irrtmp2 %3d %f %f %f     %f %f %f     %f %f %f\n",i,irr_tmp1[0],irr_tmp1[1],irr_tmp1[2],irr_tmp1[3],irr_tmp1[4],irr_tmp1[5],irr_tmp1[6],irr_tmp1[7],irr_tmp1[8]);
        irr_tmp[0] = k_one;
        irr_tmp[1] = k_zero;
        irr_tmp[2] = k_zero;
        irr_tmp[3] = -cos(crystal->gamma);
        irr_tmp[4] = sin(crystal->gamma);
        irr_tmp[5] = k_zero;
        irr_tmp[6] = k_zero;
        irr_tmp[7] = k_zero;
        irr_tmp[8] = k_one;
//fprintf(file.out,"irrtmp1 %3d %f %f %f     %f %f %f     %f %f %f\n",i,irr_tmp[0],irr_tmp[1],irr_tmp[2],irr_tmp[3],irr_tmp[4],irr_tmp[5],irr_tmp[6],irr_tmp[7],irr_tmp[8]);
        double_mat_dot(irr_tmp1, irr_tmp, symmetry->irr + 9 * i);
      }
        break;
    } // close switch(crystal->system)
  }

    //if (job->C09 == 1) {
    //fprintf(file.out,"reading09\n");
    //read_SYMMOP_crystal_09(symmetry,job,file);
   //}

   int print = 0;
   generate_operator_inverses(symmetry, job, file);
   generate_group_multiplication_table(symmetry, print, file);
   generate_group_conjugacy_classes(symmetry, print, file);
   if (job->C09 != 1) // reorder operators by conjugacy classes and recompute group multiplication table and inverses
   reorder_symmetry_operators_by_class(symmetry,job,file);
   generate_Dirac_characters(symmetry,job,file);
   generate_permutation_group_table(symmetry,job,file);
   print_symmetry_operators(symmetry,job,file);

/*
   if (job->C09 != 1) { // reorder operators by conjugacy classes and recompute group multiplication table and inverses

   int new_operator_order[symmetry->number_of_operators];
   int inr_temp[symmetry->number_of_operators][9];
   int count1 = 0;
   double irr_temp[symmetry->number_of_operators][9];
   VECTOR_DOUBLE tau_temp[symmetry->number_of_operators];
   for (i = 0; i < symmetry->number_of_operators; i++) {
     for (j = 0; j < 9; j++) {
     irr_temp[i][j] = symmetry->irr[i * 9 + j];
     tau_temp[i].comp1 = symmetry->taur[i].comp1;
     tau_temp[i].comp2 = symmetry->taur[i].comp2;
     tau_temp[i].comp3 = symmetry->taur[i].comp3;
     inr_temp[i][j] = symmetry->inr[i * 9 + j];
    }
   }
   for (i = 0; i < symmetry->number_of_classes; i++) {
     for (j = 0; j < symmetry->cls_num_k[i]; j++)  {
       new_operator_order[count1 + j] = symmetry->cls_k[symmetry->cls_pos_k[i] + j];
      }
       count1 += symmetry->cls_num_k[i];
      }
   for (i = 0; i < symmetry->number_of_operators; i++) {
     for (j = 0; j < 9; j++) {
       symmetry->irr[i * 9 + j] = irr_temp[new_operator_order[i]][j];
       symmetry->taur[i].comp1 = tau_temp[new_operator_order[i]].comp1;
       symmetry->taur[i].comp2 = tau_temp[new_operator_order[i]].comp2;
       symmetry->taur[i].comp3 = tau_temp[new_operator_order[i]].comp3;
       symmetry->inr[i * 9 + j] = inr_temp[new_operator_order[i]][j];
      }
     }

   print = 1;
   if (job->taskid == 0) print = 1;
   generate_operator_inverses(symmetry, job, file);
   generate_group_multiplication_table(symmetry, print, file);
   generate_group_conjugacy_classes(symmetry, print, file);

  } // end if (job->C09

   generate_Dirac_characters(symmetry,job,file);
   generate_permutation_group_table(symmetry,job,file);
   print_symmetry_operators(symmetry,job,file);
*/

}

void reorder_symmetry_operators_by_class(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j;
int new_operator_order[symmetry->number_of_operators];
int inr_temp[symmetry->number_of_operators][9];
int print, count1 = 0;
double irr_temp[symmetry->number_of_operators][9];
VECTOR_DOUBLE tau_temp[symmetry->number_of_operators];

  for (i = 0; i < symmetry->number_of_operators; i++) {
    for (j = 0; j < 9; j++) {
    irr_temp[i][j] = symmetry->irr[i * 9 + j];
    tau_temp[i].comp1 = symmetry->taur[i].comp1;
    tau_temp[i].comp2 = symmetry->taur[i].comp2;
    tau_temp[i].comp3 = symmetry->taur[i].comp3;
    inr_temp[i][j] = symmetry->inr[i * 9 + j];
   }
  }
  for (i = 0; i < symmetry->number_of_classes; i++) {
    for (j = 0; j < symmetry->cls_num_k[i]; j++)  {
      new_operator_order[count1 + j] = symmetry->cls_k[symmetry->cls_pos_k[i] + j];
     }
      count1 += symmetry->cls_num_k[i];
     }
  for (i = 0; i < symmetry->number_of_operators; i++) {
    for (j = 0; j < 9; j++) {
      symmetry->irr[i * 9 + j] = irr_temp[new_operator_order[i]][j];
      symmetry->taur[i].comp1 = tau_temp[new_operator_order[i]].comp1;
      symmetry->taur[i].comp2 = tau_temp[new_operator_order[i]].comp2;
      symmetry->taur[i].comp3 = tau_temp[new_operator_order[i]].comp3;
      symmetry->inr[i * 9 + j] = inr_temp[new_operator_order[i]][j];
     }
    }

  print = 0;
  generate_operator_inverses(symmetry, job, file);
  generate_group_multiplication_table(symmetry, print, file);
  generate_group_conjugacy_classes(symmetry, print, file);

}

void print_symmetry_operators(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
    
int i;

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"Symmetry Operators transformed to Cartesian basis\n");

    fprintf(file.out, "                                %2d SYMMETRY OPERATORS\n\n", symmetry->number_of_operators);
    for (i = 0; i < symmetry->number_of_operators; i++)
          fprintf(file.out,\
          "%3d %3d\n %5.2lf %5.2lf %5.2lf     %5.2lf\n %5.2lf %5.2lf %5.2lf     %5.2lf \n %5.2lf %5.2lf %5.2lf     %5.2lf\n",
          i + 1, symmetry->inverse[i] + 1, \
                 symmetry->irr[i * 9 + 0], symmetry->irr[i * 9 + 1], symmetry->irr[i * 9 + 2], symmetry->taur[i].comp1, \
                 symmetry->irr[i * 9 + 3], symmetry->irr[i * 9 + 4], symmetry->irr[i * 9 + 5], symmetry->taur[i].comp2, \
                 symmetry->irr[i * 9 + 6], symmetry->irr[i * 9 + 7], symmetry->irr[i * 9 + 8], symmetry->taur[i].comp3);
          fprintf(file.out, "\n");

    fprintf(file.out, "                                %2d SYMMETRY OPERATORS\n\n", symmetry->number_of_operators);
    for (i = 0; i < symmetry->number_of_operators; i++)
          fprintf(file.out,"%3d %3d\n %5d %5d %5d\n %5d %5d %5d\n %5d %5d %5d\n",
          i + 1, symmetry->inverse[i] + 1, symmetry->inr[i * 9 + 0], symmetry->inr[i * 9 + 1], symmetry->inr[i * 9 + 2], \
                 symmetry->inr[i * 9 + 3], symmetry->inr[i * 9 + 4], symmetry->inr[i * 9 + 5], \
                 symmetry->inr[i * 9 + 6], symmetry->inr[i * 9 + 7], symmetry->inr[i * 9 + 8]);
    fprintf(file.out, "\n");
   }

}

void count_little_k_group_operators(int k, SYMMETRY *symmetry_little_k_group, SYMMETRY *symmetry, CRYSTAL *crystal, KPOINT_TRAN *knet, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int l, t;
int *p_inr;
VECTOR_INT kvec[3];

    symmetry_little_k_group->number_of_operators = 0;
    for (t = 0; t <= job->trs; t++) {
      for (l = 0; l < symmetry->number_of_operators; l++) {
        kvec[0].comp1 =  knet->ibz[k] / (fermi->is[2] * fermi->is[1]);
        kvec[0].comp2 = (knet->ibz[k] - kvec[0].comp1 * fermi->is[2] * fermi->is[1]) / fermi->is[2];
        kvec[0].comp3 =  knet->ibz[k] - kvec[0].comp1 * fermi->is[2] * fermi->is[1] - kvec[0].comp2 * fermi->is[2];
        p_inr = symmetry->inr + symmetry->inverse[l] * 9;
        rotate_vector_int_latt(p_inr, fermi->is, t, &kvec[0], &kvec[1], crystal, job, file);
        //March2019 rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
        //fprintf(file.out,"trs %3d op %3d kvec[0] %3d %3d %3d  kvec[1] %3d %3d %3d\n",t,l,kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,\
        kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
        if (kvec[1].comp1 == kvec[0].comp1 && kvec[1].comp2 == kvec[0].comp2 && kvec[1].comp3 == kvec[0].comp3) {
        symmetry_little_k_group->number_of_operators++;
       }
      }
     }

}

void generate_little_k_group_operators(int k, SYMMETRY *symmetry_little_k_group, SYMMETRY *symmetry, CRYSTAL *crystal,KPOINT_TRAN *knet, FERMI *fermi, JOB_PARAM *job, FILES file)
//void generate_little_k_group_operators(int k, SYMMETRY *symmetry_little_k_group, SYMMETRY *symmetry, KPOINT_TRAN *knet, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, j, l, t, count, count_operators;
int *p_inr, *p_ind_i, *p_ind_j, *p_num, *p_ost;
int *q_inr, *q_ind_i, *q_ind_j, *q_num, *q_ost;
double *p_rot, *q_rot;
VECTOR_INT kvec[3];

  p_ind_i = symmetry->ind_i;
  p_ind_j = symmetry->ind_j;
  p_rot = symmetry->rot;
  p_ost = symmetry->op_shift;
  p_num = symmetry->num_ij;

  q_ind_i = symmetry_little_k_group->ind_i;
  q_ind_j = symmetry_little_k_group->ind_j;
  q_rot = symmetry_little_k_group->rot;
  q_ost = symmetry_little_k_group->op_shift;
  q_num = symmetry_little_k_group->num_ij;

  count = 0;
  count_operators = 0;
  for (t = 0; t <= job->trs; t++) {
  for (l = 0; l < symmetry->number_of_operators; l++) {
    kvec[0].comp1 =  knet->ibz[k] / (fermi->is[2] * fermi->is[1]);
    kvec[0].comp2 = (knet->ibz[k] - kvec[0].comp1 * fermi->is[2] * fermi->is[1]) / fermi->is[2];
    kvec[0].comp3 =  knet->ibz[k] - kvec[0].comp1 * fermi->is[2] * fermi->is[1] - kvec[0].comp2 * fermi->is[2];
    p_inr = symmetry->inr + symmetry->inverse[l] * 9;
    rotate_vector_int_latt(p_inr, fermi->is, t, &kvec[0], &kvec[1], crystal, job, file);
    //rotate_vector_int(p_inr, &kvec[0], &kvec[1]);
    if (kvec[1].comp1 == kvec[0].comp1 && kvec[1].comp2 == kvec[0].comp2 && kvec[1].comp3 == kvec[0].comp3) {
    for (i = 0; i < 9; i++) {
    symmetry_little_k_group->inr[count_operators * 9 + i] = symmetry->inr[l * 9 + i]; 
    symmetry_little_k_group->irr[count_operators * 9 + i] = symmetry->irr[l * 9 + i];
   }
    symmetry_little_k_group->taur[count_operators].comp1 = symmetry->taur[l].comp1;
    symmetry_little_k_group->taur[count_operators].comp2 = symmetry->taur[l].comp2;
    symmetry_little_k_group->taur[count_operators].comp3 = symmetry->taur[l].comp3;
    for (i = 0; i < job->l_max + 2; i++) {
      *q_num = *p_num;
      *q_ost = count;
      for (j = 0; j < *p_num; j++) {
        *q_rot   = *p_rot;
        *q_ind_i = *p_ind_i;
        *q_ind_j = *p_ind_j;
         p_ind_i++;
         p_ind_j++;
         p_rot++;
         q_ind_i++;
         q_ind_j++;
         q_rot++;
        }
       count += *p_num;
       p_num++;
       q_ost++;
       q_num++;
      }
     count_operators++;
    } // end if
    else {
    for (i = 0; i < job->l_max + 2; i++) {
      for (j = 0; j < *p_num; j++) {
        //fprintf(file.out, "symm %3d ind_i %3d %3d num_ij %3d %3d rot  %10.4lf\n", i, *p_ind_i, *p_ind_j, *p_num, *p_ost, *p_rot);
        p_ind_i++;
        p_ind_j++;
        p_rot++;
      }
      //fprintf(file.out, "\n");
      p_ost++;
      p_num++;
     }
    //fprintf(file.out, "\n");
   }
  } // close loop on l
 } // close loop on t
  //fprintf(file.out, "\n");

}
    
void generate_Dirac_characters(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
    
    // generate Dirac characters

    // See Group Theory with Applications in Chemical Physics, P. Jacobs, CUP (2005) A2 Class Algebra
    // or Symmetry Principles in Solid State and Molecular Physics, M. Lax, Wiley (1974) Sec 1.3 for details

    //fprintf(file.out,"DIRAC CHARACTERS\n\n");

  int dim = symmetry->number_of_classes;
  int i, j, k, l, pos_i, pos_j, count;
  int info = 0;
  int inv_cls_k[symmetry->number_of_operators];
  int degeneracy[symmetry->number_of_classes];
  int character_product[symmetry->number_of_operators];
  int class_constants[symmetry->number_of_classes][symmetry->number_of_classes][symmetry->number_of_classes];
  char trans = 'T', no_trans = 'N';
  char uplo = 'U';
  char jobz = 'V';
  double *eigenvalues_imag;
  double alpha = k_one, beta = k_zero;
  double sum[symmetry->number_of_classes];
  DoubleMatrix *eigenvalues_real, *eigval, *tmp, *tmp1;
  DoubleMatrix *matrix, *eigenvectors, *eigvec;

  AllocateDoubleMatrix(&matrix,&dim,&dim,job);
  AllocateDoubleMatrix(&eigenvalues_real,&dim,&dim,job);
  AllocateDoubleMatrix(&eigval,&dim,&dim,job);
  AllocateDoubleArray(&eigenvalues_imag,&dim,job);
  AllocateDoubleMatrix(&eigenvectors,&dim,&dim,job);
  AllocateDoubleMatrix(&eigvec,&dim,&dim,job);
  AllocateDoubleMatrix(&tmp,&dim,&dim,job);
  AllocateDoubleMatrix(&tmp1,&dim,&dim,job);

  for (i = 0; i < symmetry->number_of_operators; i++) { 
    inv_cls_k[i] = symmetry->inverse[symmetry->cls_k[i]];
   }
  for (i = 0; i < symmetry->number_of_classes; i++) { 
    for (j = 0; j < symmetry->number_of_classes; j++) { 
      for (k = 0; k < symmetry->number_of_classes; k++) { 
        class_constants[i][j][k] = 0;
       }
      }
     }

  for (i = 0; i < symmetry->number_of_classes; i++) { 
    pos_i = symmetry->cls_pos_k[i];
    for (j = 0; j < symmetry->number_of_classes; j++) { 
      pos_j = symmetry->cls_pos_k[j];
      for (k = 0; k < symmetry->number_of_operators; k++) 
      character_product[k] = 0;
      for (k = 0; k < symmetry->cls_num_k[i]; k++) { 
        for (l = 0; l < symmetry->cls_num_k[j]; l++) { 
          (character_product[symmetry->grp_k[symmetry->cls_k[pos_i + k] * symmetry->number_of_operators + inv_cls_k[pos_j + l]]])++;
          //fprintf(file.out,"%3d %3d  %3d %3d   %3d\n",\
          i+1,j+1,k,l,symmetry->grp_k[symmetry->cls_k[pos_i + k] * symmetry->number_of_operators + inv_cls_k[pos_j + l]]);
         }
        }
       //fprintf(file.out,"\n");
       count = 0;
       for (k = 0; k < symmetry->number_of_classes; k++) { 
         for (l = 0; l < symmetry->cls_num_k[k]; l++) { 
           class_constants[i][j][k] += character_product[count];
           count++;
          }
         }
        }
       }
  for (i = 0; i < symmetry->number_of_classes; i++) { 
    for (j = 0; j < symmetry->number_of_classes; j++) { 
      for (k = 0; k < symmetry->number_of_classes; k++) { 
        class_constants[i][j][k] /= symmetry->cls_num_k[k];
       }
      }
     }

  if (job->verbosity > 1 && job->taskid == 0) {
  fprintf(file.out,"GROUP CLASS CONSTANTS\n\n");
  for (k = 0; k < symmetry->number_of_classes; k++) { 
    fprintf(file.out,"k = %3d\n\n",k);
    for (i = 0; i < symmetry->number_of_classes; i++) { 
      for (j = 0; j < symmetry->number_of_classes; j++) { 
        fprintf(file.out,"%3d ",class_constants[i][j][k]);
       }
       fprintf(file.out,"\n");
      }
       fprintf(file.out,"\n");
      }
     }

  // Identify which eigenvalues belong to each representation using eigenvectors with no eigenvalue degeneracy

    ResetDoubleMatrix(matrix);
  for (i = 0; i < symmetry->number_of_classes; i++) {
    degeneracy[i] = 0;
    for (j = 0; j < symmetry->number_of_classes; j++) { 
      for (k = 0; k < symmetry->number_of_classes; k++) { 
        //matrix->a[j][k] = (double) class_constants[i][j][k];
        matrix->a[j][k] += sqrt (i + 3.4) * (double) class_constants[i][j][k];
        //matrix1[i]->a[j][k] += 0.277 * sqrt(i + 1) * (double) class_constants[i][j][k];
        //matrix->a[j][k] += class_constants[i][j][k] / two;
        //matrix->a[k][j] += class_constants[i][j][k] / two;
        ResetDoubleMatrix(eigenvectors);
       }
      }
       //print_double_matrix(matrix, file); fflush(file.out);
       //DiagonaliseSymmetrical(&matrix, &eigenvalues_real->a[i], &eigenvectors, &jobz, &uplo, &info);
       if (i == symmetry->number_of_classes - 1) 
       DiagonaliseRealGeneral(&matrix, &eigenvalues_real->a[i], &eigenvalues_imag, &eigenvectors, &info);
       //print_real_eigenvector_matrix(eigenvectors, eigenvalues_real->a[i], file); fflush(file.out);
       for (j = 0; j < symmetry->number_of_classes; j++) { 
         for (k = 0; k < j; k++) { 
           //fprintf(file.out,"%3d %3d %f %f\n",j,k,eigenvalues_real->a[i][j] , eigenvalues_real->a[i][k]); fflush(file.out);
           if (fabs(eigenvalues_real->a[i][j] - eigenvalues_real->a[i][k]) < 0.00001) (degeneracy[i])++;
          }
         }
          //fprintf(file.out,"DEGENERACY %3d\n",degeneracy[i]); fflush(file.out);
          if (degeneracy[i] == 0) {
          for (j = 0; j < symmetry->number_of_classes; j++) {
            for (k = 0; k < symmetry->number_of_classes; k++) {
              eigvec->a[j][k] = eigenvectors->a[j][k];
             }
            }
            //print_real_eigenvector_matrix(eigvec, eigenvalues_real->a[i], file);
           }
          } // close loop on i

  for (i = 0; i < symmetry->number_of_classes; i++) sum[i] = k_zero;
    ResetDoubleMatrix(eigval);
    for (i = 0; i < symmetry->number_of_classes; i++) {
      for (j = 0; j < symmetry->number_of_classes; j++) { 
        for (k = 0; k < symmetry->number_of_classes; k++) { 
          matrix->a[j][k] = class_constants[i][j][k];
         }
        }
       ResetDoubleMatrix(tmp);
       ResetDoubleMatrix(tmp1);
       //print_double_matrix(matrix, file);
       //print_double_matrix(eigvec, file);
       DoubleGEMM(&no_trans, &no_trans, &alpha, &eigvec, &matrix, &beta, &tmp);
       //print_double_matrix(tmp, file);
       DoubleGEMM(&no_trans, &trans, &alpha, &tmp, &eigvec, &beta, &tmp1);
       //print_double_matrix(tmp1, file);
       for (j = 0; j < symmetry->number_of_classes; j++) { 
         eigval->a[i][j] = tmp1->a[j][j];
         sum[j] += tmp1->a[j][j] * tmp1->a[j][j] / (double) symmetry->cls_num_k[i];
        }
         //print_real_eigenvector_matrix(eigvec, eigval->a[i], file);
        } // close loop on i

   for (i = 0; i < symmetry->number_of_classes; i++) {
     symmetry->irp_dim_k[i] = (int) sqrt((double) symmetry->grp_dim / sum[i] + 1e-04);
     //fprintf(file.out,"%3d %10.4lf %3d %20.14f %3d\n",i,sum[i],symmetry->grp_dim,sqrt((double) symmetry->grp_dim / sum[i]), \
     symmetry->irp_dim_k[i]);
    }
   
  // generate character table

  for (i = 0; i < symmetry->number_of_classes; i++) { 
  //fprintf(file.out,"|  %2d   |    %2d     |",i + 1, symmetry->irp_dim_k[i]);
  for (j = 0; j < symmetry->number_of_classes; j++) { 
  symmetry->character_table[i * symmetry->number_of_classes + j] = (double) symmetry->irp_dim_k[i] / \
  (double) symmetry->cls_num_k[j] * eigval->a[j][i];
  //fprintf(file.out,"%3d %3d %f\n",symmetry->irp_dim_k[i] , symmetry->cls_num_k[j] , eigval->a[j][i]);
  //fprintf(file.out,"%4.0lf  |",symmetry->character_table[i * symmetry->number_of_classes + j]);
   }
  //fprintf(file.out,"\n");
   }

  if (job->taskid >= 0 && job->verbosity > 1) {
  fprintf(file.out,"---------------------");
  for (i = 0; i < symmetry->number_of_classes; i++) fprintf(file.out,"-------");
  fprintf(file.out,"\n");
  fprintf(file.out,"| CHARACTER TABLE   |");
  for (i = 0; i < symmetry->number_of_classes - 1; i++) fprintf(file.out,"       ");
  fprintf(file.out,"      |\n");
  fprintf(file.out,"---------------------");
  for (i = 0; i < symmetry->number_of_classes; i++) fprintf(file.out,"-------");
  fprintf(file.out,"\n");
  fprintf(file.out,"| IRREP | DIMENSION |");
  for (i = 0; i < symmetry->number_of_classes; i++) fprintf(file.out,"  %2d  |",symmetry->cls_num_k[i]);
  fprintf(file.out,"\n");
  fprintf(file.out,"---------------------");
  for (i = 0; i < symmetry->number_of_classes; i++) fprintf(file.out,"-------");
  fprintf(file.out,"\n");
  for (i = 0; i < symmetry->number_of_classes; i++) { 
  fprintf(file.out,"|  %2d   |    %2d     |",i + 1, symmetry->irp_dim_k[i]);
  for (j = 0; j < symmetry->number_of_classes; j++) { 
  ////symmetry->character_table[i * symmetry->number_of_classes + j] = (double) symmetry->irp_dim_k[i] / \
  (double) symmetry->cls_num_k[j] * eigval->a[j][i];
  //fprintf(file.out,"%3d %3d %f\n",symmetry->irp_dim_k[i] , symmetry->cls_num_k[j] , eigval->a[j][i]);
  fprintf(file.out,"%4.0lf  |",symmetry->character_table[i * symmetry->number_of_classes + j]);
   }
  fprintf(file.out,"\n");
   }
  fprintf(file.out,"---------------------");
  for (i = 0; i < symmetry->number_of_classes; i++) fprintf(file.out,"-------");
  fprintf(file.out,"\n");
  fprintf(file.out,"\n\n");
  //fflush(file.out);
 }

  DestroyDoubleMatrix(&matrix,job);
  DestroyDoubleMatrix(&eigenvalues_real,job);
  DestroyDoubleArray(&eigenvalues_imag,&dim,job);
  DestroyDoubleMatrix(&eigval,job);
  DestroyDoubleMatrix(&eigenvectors,job);
  DestroyDoubleMatrix(&eigvec,job);
  DestroyDoubleMatrix(&tmp,job);
  DestroyDoubleMatrix(&tmp1,job);

}

void generate_permutation_group_table(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
    
    // populate permutation group multiplication table

      symmetry->grp_pm[0]  = 0;
      symmetry->grp_pm[1]  = 1;
      symmetry->grp_pm[2]  = 2;
      symmetry->grp_pm[3]  = 3;
      symmetry->grp_pm[4]  = 4;
      symmetry->grp_pm[5]  = 5;
      symmetry->grp_pm[6]  = 6;
      symmetry->grp_pm[7]  = 7;

      symmetry->grp_pm[8]  = 1;
      symmetry->grp_pm[9]  = 0;
      symmetry->grp_pm[10] = 3;
      symmetry->grp_pm[11] = 2;
      symmetry->grp_pm[12] = 6;
      symmetry->grp_pm[13] = 7;
      symmetry->grp_pm[14] = 4;
      symmetry->grp_pm[15] = 5;

      symmetry->grp_pm[16] = 2;
      symmetry->grp_pm[17] = 3;
      symmetry->grp_pm[18] = 0;
      symmetry->grp_pm[19] = 1;
      symmetry->grp_pm[20] = 5;
      symmetry->grp_pm[21] = 4;
      symmetry->grp_pm[22] = 7;
      symmetry->grp_pm[23] = 6;

      symmetry->grp_pm[24] = 3;
      symmetry->grp_pm[25] = 2;
      symmetry->grp_pm[26] = 1;
      symmetry->grp_pm[27] = 0;
      symmetry->grp_pm[28] = 7;
      symmetry->grp_pm[29] = 6;
      symmetry->grp_pm[30] = 5;
      symmetry->grp_pm[31] = 4;

      symmetry->grp_pm[32] = 4;
      symmetry->grp_pm[33] = 5;
      symmetry->grp_pm[34] = 6;
      symmetry->grp_pm[35] = 7;
      symmetry->grp_pm[36] = 0;
      symmetry->grp_pm[37] = 1;
      symmetry->grp_pm[38] = 2;
      symmetry->grp_pm[39] = 3;

      symmetry->grp_pm[40] = 5;
      symmetry->grp_pm[41] = 4;
      symmetry->grp_pm[42] = 7;
      symmetry->grp_pm[43] = 6;
      symmetry->grp_pm[44] = 2;
      symmetry->grp_pm[45] = 3;
      symmetry->grp_pm[46] = 0;
      symmetry->grp_pm[47] = 1;

      symmetry->grp_pm[48] = 6;
      symmetry->grp_pm[49] = 7;
      symmetry->grp_pm[50] = 4;
      symmetry->grp_pm[51] = 5;
      symmetry->grp_pm[52] = 1;
      symmetry->grp_pm[53] = 0;
      symmetry->grp_pm[54] = 3;
      symmetry->grp_pm[55] = 2;

      symmetry->grp_pm[56] = 7;
      symmetry->grp_pm[57] = 6;
      symmetry->grp_pm[58] = 5;
      symmetry->grp_pm[59] = 4;
      symmetry->grp_pm[60] = 3;
      symmetry->grp_pm[61] = 2;
      symmetry->grp_pm[62] = 1;
      symmetry->grp_pm[63] = 0;

}

void generate_group_multiplication_table(SYMMETRY *symmetry, int print, FILES file)

{
    
int i, j, k;
double irr1[9];

    // generate group multiplication table

    if (print) fprintf(file.out,"GROUP MULTIPLICATION TABLE\n\n");
    if (print) fprintf(file.out,"      ");
    for (i = 0; i < symmetry->number_of_operators * symmetry->number_of_operators; i++) 
    //for (i = 0; i < symmetry->number_of_operators * 9; i++) 
    symmetry->grp_k[i] = -1;
    for (i = 0; i < symmetry->number_of_operators; i++) 
    if (print) fprintf(file.out,"%2d ",i);
    if (print) fprintf(file.out,"\n\n");
    for (i = 0; i < symmetry->number_of_operators; i++) {
     if (print) fprintf(file.out,"%3d   ",i);
      for (j = 0; j < symmetry->number_of_operators; j++) {
        double_mat_dot(symmetry->irr + 9 * i, symmetry->irr + 9 * j, irr1);
        //double_mat_dot(symmetry->irr + 9 * symmetry->inverse[i], symmetry->irr + 9 * symmetry->inverse[j], irr1);
        for (k = 0; k < symmetry->number_of_operators; k++) {
          if (check_mat(&symmetry->irr[9 * k], irr1) == 1) {
          symmetry->grp_k[i * symmetry->number_of_operators + j] = k;
          if (print) fprintf(file.out,"%2d ",k);
        } // close if (check_mat ...
      } // close loop on k
    } // close loop on j
   if (print) fprintf(file.out,"\n");
  } // close loop on i
   if (print) fprintf(file.out,"\n");

}

void generate_group_conjugacy_classes(SYMMETRY *symmetry, int print, FILES file)

{
    // generate group conjugacy classes

int i, j, k;

    int taken[symmetry->number_of_operators];
    int total_taken = 0;
    for (i = 0; i < symmetry->number_of_operators; i++) taken[i] = -1;
    for (i = 0; i < symmetry->number_of_operators; i++) symmetry->cls_pos_k[0] = 0;
    symmetry->grp_dim = symmetry->number_of_operators;

    if (print) fprintf(file.out,"GROUP CONJUGACY CLASSES\n\n");
    symmetry->number_of_classes = 0;
    for (i = 0; i < symmetry->number_of_operators; i++) {
      if (taken[i] >= 0) continue;
        symmetry->cls_pos_k[symmetry->number_of_classes] = total_taken; 
        symmetry->cls_num_k[symmetry->number_of_classes] = 0;
        for (j = 0; j < symmetry->number_of_operators; j++) {
          k = symmetry->grp_k[symmetry->inverse[j] * symmetry->number_of_operators + \
          symmetry->grp_k[i * symmetry->number_of_operators + j]];
          if (taken[k] == -1) { 
          (symmetry->cls_num_k[symmetry->number_of_classes])++; 
           symmetry->cls_k[total_taken] = k; 
          (taken[k])++; 
           total_taken++;
         }
           //fprintf(file.out,"%2d %2d %2d %2d\n",i, j, k,total_taken);
         } // close loop on j
          if (symmetry->cls_num_k[symmetry->number_of_classes] > 0) symmetry->number_of_classes++;
        } // close loop on i
          if (print) {
        fprintf(file.out,"CLASS | NO. IN CLASS | POSITION | CLASS MEMBERS\n");
        for (i = 0; i < symmetry->number_of_classes; i++) {
          fprintf(file.out,"%3d        %3d          %3d      ",i,symmetry->cls_num_k[i],symmetry->cls_pos_k[i]);
          for (j = 0; j < symmetry->cls_num_k[i]; j++)  {
            fprintf(file.out," %2d ",symmetry->cls_k[symmetry->cls_pos_k[i] + j]);
	   }
          fprintf(file.out,"\n");
         }
          fprintf(file.out,"\n\n");
       }

}
  
void generate_operator_inverses(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, k;
double *p_irr, *q_irr;

  // Compute inverses of each symmetry operator

  if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "INVERSES OF SYMMETRY OPERATORS\n");
  if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "OPERATOR INVERSE \n");
  for (i = 0; i < symmetry->number_of_operators; i++) {
    for (j = 0; j < symmetry->number_of_operators; j++) {
      p_irr = symmetry->irr + i * 9;
      q_irr = symmetry->irr + j * 9;
      k = check_inverse1(p_irr, q_irr);
      //if (k == 1 && debug == 1) {
      if (job->taskid == 0 && k == 1 && job->verbosity > 1) {
        fprintf(file.out, "%5d %6d \n", i + 1, j + 1);
      }
      if (k == 1) {
        symmetry->inverse[i] = j;
        break;
      }
    }
  }
   //fprintf(file.out,"\n\n");

}

void generate_little_k_group(int k, SYMMETRY *symmetry_little_k_group, FERMI *fermi, KPOINT_TRAN *knet, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int print = 0;

  generate_little_k_group_operators(k,symmetry_little_k_group,symmetry,crystal,knet,fermi,job,file);
  //CHP March2019 generate_little_k_group_operators(k,symmetry_little_k_group,symmetry,knet,fermi,job,file);
  generate_operator_inverses(symmetry_little_k_group, job, file);
  generate_group_multiplication_table(symmetry_little_k_group, print, file);
  generate_group_conjugacy_classes(symmetry_little_k_group, print, file);
  reorder_symmetry_operators_by_class(symmetry_little_k_group,job,file);
  generate_Dirac_characters(symmetry_little_k_group,job,file);
  generate_permutation_group_table(symmetry_little_k_group,job,file);
  print_symmetry_operators(symmetry_little_k_group,job,file);

}
