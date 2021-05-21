/*
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include "mylogical.h"
#include "mycomplex.h"
#include "conversion_factors.h"
#include "myconstants.h"
*/
#include "USER_DATA.h"
#include "ERRORS.h"

void knet_check(int is[3], char *jobname, CRYSTAL *crystal, FILES file)

{

      switch (crystal->type[0]) {

        case 'C':

      if (is[0] < 1 || is[0] > 50 || is[1] < 1 || is[1] > 50 || is[2] < 1 || is[2] > 50) {
      fprintf(file.out,"ERROR: K POINT NET FOR %s INCORRECT. VALUES ENTERED (%d %d %d) MUST BE IN RANGE 1 to 50.",
      jobname,is[0],is[1],is[2]);
      exit(1);
     }

        break;

        case 'S':

      if (is[0] < 1 || is[0] > 50 || is[1] < 1 || is[1] > 50 || is[2] != 1) {
      fprintf(file.out,"ERROR: K POINT NET FOR %s INCORRECT. \
      VALUES ENTERED (%d %d %d) MUST BE IN RANGE 1 to 50 FOR FIRST TWO ENTRIES AND 1 FOR THE THIRD ENTRY.",
      jobname,is[0],is[1],is[2]);
      exit(1);
     }

        break;

        case 'P':

      if (is[0] != 1 || is[1] != 1 || is[2] < 1 || is[2] > 50) {
      fprintf(file.out,"ERROR: K POINT NET FOR %s INCORRECT. \
      VALUES ENTERED (%d %d %d) MUST BE IN RANGE 1 to 50 FOR THIRD ENTRY AND 1 FOR THE FIRST TWO ENTRIES.",
      jobname,is[0],is[1],is[2]);
      exit(1);
     }

        break;

      } // close switch

}

void int_range_check(int range[2], int value[2], char *range_type, char *jobname, FILES file)

{

     if (range[0] < value[0] || range[0] > range[1] || range[1] > value[1]) {
     fprintf(file.out,"ERROR: %s FOR %s INCORRECT. VALUES ENTERED (%d %d) MUST BE IN RANGE %d to %d.",range_type,jobname,range[0],range[1],
     value[0],value[1]);
     exit(1);
    }

}
void double_range_check(double range[2], double value[2], char *range_type, char *jobname, FILES file)

{

     if (range[0] < value[0] || range[0] >= range[1] || range[1] > value[1]) {
     fprintf(file.out,"ERROR: %s FOR %s INCORRECT. VALUES ENTERED (%4.1e %4.1e) MUST BE IN RANGE (%4.1e to %4.1e).",
     range_type,jobname,range[0],range[1],value[0],value[1]);
     exit(1);
    }

}
