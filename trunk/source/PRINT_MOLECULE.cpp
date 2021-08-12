

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
#include "mycomplex.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "PRINT_MOLECULE.h"

using namespace std;

void print_double_matrix(DoubleMatrix *a, FILES file){

  int i, j, k, blocks, remndr;
  blocks = a->iCols / 12;
  remndr = a->iCols - blocks * 12;
  for (k=0;k<blocks;k++){
    for (j=0;j<12;j++){
      fprintf(file.out,"    %3d   ",k * 12 + j + 1);
     }
    fprintf(file.out,"\n\n");
    for (i=0;i<a->iRows;i++){
      //for (j=0;j<a->iCols;j++){
      for (j=0;j<12;j++){
      //fprintf(file.out,"%11.6e  ",a->a[i][k * 12 + j]);
      fprintf(file.out,"%9.5lf  ",a->a[i][k * 12 + j]);
     }
    fprintf(file.out,"\n");
   }
      fprintf(file.out,"\n\n\n");
    }
   if(remndr > 0) {
    for (j=0;j<remndr;j++){
      fprintf(file.out,"    %3d   ",k * 12 + j + 1);
     }
    fprintf(file.out,"\n\n");
     for (i=0;i<a->iRows;i++){
       for(j = 0; j < remndr; j++) {
         //fprintf(file.out,"%11.6e  ",a->a[i][k * 12 + j]);
         fprintf(file.out,"%9.5lf  ",a->a[i][k * 12 + j]);
        }
       fprintf(file.out,"\n");
      }
     fprintf(file.out,"\n\n");
    }

}

void print_complex_eigenvector_matrix3(ComplexMatrix *a, double *e, FILES file){

  int i, j, k, blocks, remndr;
  blocks = a->iCols / 6;
  remndr = a->iCols - blocks * 6;
  for (k=0;k<blocks;k++){
    for (j=0;j<6;j++){
      fprintf(file.out,"     %5d        ",k * 6 + j + 1);
     }
       fprintf(file.out,"\n\n");
    for (j=0;j<6;j++){
      fprintf(file.out," %14.6e   ",*(e + k * 6 + j));
     }
    fprintf(file.out,"\n\n");
    for (i=0;i<a->iRows;i++){
      for (j=0;j<6;j++){
      fprintf(file.out,"%8.3lf%8.3lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
     }
    fprintf(file.out,"\n");
   }
      fprintf(file.out,"\n\n\n");
    }
  if (remndr > 0) {
    for (j=0;j<remndr;j++){
      fprintf(file.out,"     %5d        ",k * 6 + j + 1);
     }
       fprintf(file.out,"\n\n");
    for (j=0;j<remndr;j++){
      fprintf(file.out,"  %14.6e   ",*(e + k * 6 + j));
     }
    fprintf(file.out,"\n\n");
     for (i=0;i<a->iRows;i++){
       for(j = 0; j < remndr; j++) {
         fprintf(file.out,"%8.3lf%8.3lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
        }
       fprintf(file.out,"\n");
      }
     fprintf(file.out,"\n\n");
    }

}

void print_complex_matrix(ComplexMatrix *a, FILES file)

{

  int i, j, k, blocks, remndr;

  blocks = a->iCols / 6;
  remndr = a->iCols - blocks * 6;
  for (k=0;k<blocks;k++){
    for (j=0;j<6;j++){
      fprintf(file.out,"     %5d        ",k * 6 + j + 1);
     }
    fprintf(file.out,"\n\n");
    for (i=0;i<a->iRows;i++){
      for (j=0;j<6;j++){
      fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
     }
    fprintf(file.out,"\n");
   }
      fprintf(file.out,"\n\n\n");
    }
  if (remndr > 0) {
    for (j=0;j<remndr;j++){
      fprintf(file.out,"     %5d        ",k * 6 + j + 1);
     }
    fprintf(file.out,"\n\n");
     for (i=0;i<a->iRows;i++){
       for(j = 0; j < remndr; j++) {
         //fprintf(file.out,"%15.10lf  ",a->a[i][k * 6 + j].real());
         fprintf(file.out,"%11.7lf%9.5lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
         //fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
        }
       fprintf(file.out,"\n");
      }
     fprintf(file.out,"\n\n");
    }

}

void print_double_matrix2(DoubleMatrix *a, int transpose, int cols, double scale_factor, FILES file)

{

  int i, j, k, blocks, remndr, count;

  // print formatted matrix with cols columns. transpose == 0 => not transposed, transpose == 1 => transposed
  // matrix elements are scaled by scale_factor, e.g. au_to_eV

  if (transpose == 0) {
  blocks = a->iCols / cols;
  remndr = a->iCols - blocks * cols;
 }

  if (transpose == 1) {
  blocks = a->iRows / cols;
  remndr = a->iRows - blocks * cols;
  fprintf(file.out,"Matrix Transpose requested\n\n");
 }

  for (k=0;k<blocks;k++){
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|  ");
    for (j=0;j<cols;j++){
      fprintf(file.out," %5d     ",k * cols + j + 1);
     }
    fprintf(file.out,"   |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  if (transpose == 0) {
    count = 1;
    for (i=0;i<a->iRows;i++){
      fprintf(file.out,"|%3d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%9.5lf  ",a->a[i][k * cols + j] * scale_factor);
     }
    fprintf(file.out,"  |\n");
   }
  }
  if (transpose == 1) {
      count = 1;
    for (i=0;i<a->iCols;i++){
      fprintf(file.out,"|%3d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%9.5lf  ",a->a[k * cols + j][i] * scale_factor);
     }
    fprintf(file.out,"  |\n");
   }
  }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

  if (remndr > 0) {
    fprintf(file.out,"|  ");
    for (j=0;j<remndr;j++)
      fprintf(file.out,"  %5d    ",k * cols + j + 1);
    for (j=0;j<cols - remndr;j++)
      fprintf(file.out,"                    ");
    fprintf(file.out,"   |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  if (transpose == 0) {
      count = 1;
     for (i=0;i<a->iRows;i++){
      fprintf(file.out,"|%3d",count);
      count++;
       for (j = 0; j < remndr; j++)
         fprintf(file.out,"%9.5lf  ",a->a[i][k * cols + j] * scale_factor);
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"                    ");
       fprintf(file.out,"  |\n");
      }
     }
  if (transpose == 1) {
      count = 1;
     for (i=0;i<a->iCols;i++){
      fprintf(file.out,"|%3d",count);
      count++;
       fprintf(file.out,"|   ");
       for(j = 0; j < remndr; j++) 
         fprintf(file.out,"%9.5lf  ",a->a[k * cols + j][i] * scale_factor);
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"                    ");
       fprintf(file.out,"  |\n");
      }
     }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

}
