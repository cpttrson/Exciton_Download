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
#include "PRINT_UTIL.h"

using namespace std;

void PrintDoubleMatrix(DoubleMatrix *a, FILES file){
  int i, j;
  for (i=0;i<a->iRows;i++){
    for (j=0;j<a->iCols;j++){
      //fprintf(file.out,"%15.12lf",a->a[i][j]);
      fprintf(file.out,"%6.2lf  ",a->a[i][j]);
    }
    fprintf(file.out,"\n");
  }
    fprintf(file.out,"\n");
}

/*
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
      fprintf(file.out,"%11.6e  ",a->a[i][k * 12 + j]);
      //fprintf(file.out,"%8.3lf  ",a->a[i][k * 12 + j]);
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
         fprintf(file.out,"%11.6e  ",a->a[i][k * 12 + j]);
         //fprintf(file.out,"%8.3lf  ",a->a[i][k * 12 + j]);
        }
       fprintf(file.out,"\n");
      }
     fprintf(file.out,"\n\n");
    }

}
*/

void print_real_matrix2(DoubleMatrix *a, int transpose, int cols, double scale_factor, FILES file)

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
 }

  for (k=0;k<blocks;k++){
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|  ");
    for (j=0;j<cols;j++){
      fprintf(file.out,"%12d",k * cols + j + 1);
     }
    fprintf(file.out,"       |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  if (transpose == 0) {
    count = 1;
    for (i=0;i<a->iRows;i++){
      fprintf(file.out,"|%4d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%12.5lf",a->a[i][k * cols + j] * scale_factor);
     }
    fprintf(file.out,"     |\n");
   }
  }
  if (transpose == 1) {
      count = 1;
    for (i=0;i<a->iCols;i++){
      fprintf(file.out,"|%4d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%12.5lf",a->a[k * cols + j][i] * scale_factor);
     }
    fprintf(file.out,"     |\n");
   }
  }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

  if (remndr > 0) {
    fprintf(file.out,"|  ");
    for (j=0;j<remndr;j++)
      fprintf(file.out,"%12d",k * cols + j + 1);
    for (j=0;j<cols - remndr;j++)
      fprintf(file.out,"            ");
    fprintf(file.out,"       |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  if (transpose == 0) {
      count = 1;
     for (i=0;i<a->iRows;i++){
      fprintf(file.out,"|%4d",count);
      count++;
       for (j = 0; j < remndr; j++)
         fprintf(file.out,"%12.5lf",a->a[i][k * cols + j] * scale_factor);
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"            ");
       fprintf(file.out,"     |\n");
      }
     }
  if (transpose == 1) {
      count = 1;
     for (i=0;i<a->iCols;i++){
      fprintf(file.out,"|%4d",count);
      count++;
       //fprintf(file.out,"|   ");
       for(j = 0; j < remndr; j++) 
         fprintf(file.out,"%12.5lf",a->a[k * cols + j][i] * scale_factor);
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"            ");
       fprintf(file.out,"     |\n");
      }
     }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

}

void print_complex_eigenvector_matrix(ComplexMatrix *a, double *e, FILES file){
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

void print_complex_eigenvector_matrix1(ComplexMatrix *a, double *e, int cols, double scale_factor, FILES file)

{

  int i, j, k, blocks, remndr, count;

  // print formatted matrix with cols columns. transpose == 0 => not transposed, transpose == 1 => transposed
  // matrix elements are scaled by scale_factor, e.g. au_to_eV

  blocks = a->iCols / cols;
  remndr = a->iCols - blocks * cols;

  for (k=0;k<blocks;k++){
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| ");
    for (j=0;j<cols;j++){
      fprintf(file.out,"   %2d  %8.3lf     ",k * cols + j + 1,*(e + k * cols + j) * scale_factor);
     }
    fprintf(file.out,"    |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

    count = 1;
    for (i=0;i<a->iRows;i++){
      fprintf(file.out,"| %2d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * cols + j].real(),a->a[i][k * cols + j].imag());
     }
    fprintf(file.out,"  |\n");
   }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

  if (remndr > 0) {
    fprintf(file.out,"| ");
    for (j=0;j<remndr;j++)
      fprintf(file.out,"   %2d  %8.3lf     ",k * cols + j + 1,*(e + k * cols + j) * scale_factor);
    for (j=0;j<cols - remndr;j++)
      fprintf(file.out,"                    ");
    fprintf(file.out,"    |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

      count = 1;
     for (i=0;i<a->iRows;i++){
      fprintf(file.out,"| %2d",count);
      count++;
       for (j = 0; j < remndr; j++)
         fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * cols + j].real(),a->a[i][k * cols + j].imag());
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"                    ");
       fprintf(file.out,"  |\n");
      }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

}

void print_complex_eigenvector_matrix2(ComplexMatrix *a, double *e, int cols, int last_vector, double scale_factor, FILES file)

{

  int i, j, k, blocks, remndr, count;

  // print formatted matrix with cols columns. Print eigenvectors up to last_vector.

  blocks = a->iRows / cols;
  remndr = a->iRows - blocks * cols;
  if (last_vector < a->iRows) {
  blocks = last_vector / cols;
  remndr = last_vector - blocks * cols;
 }

  for (k=0;k<blocks;k++){
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| ");
    for (j=0;j<cols;j++){
      fprintf(file.out,"  %3d  %8.3lf     ",k * cols + j + 1,*(e + k * cols + j) * scale_factor);
     }
    fprintf(file.out,"    |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

      count = 1;
    for (i=0;i<a->iCols;i++){
      fprintf(file.out,"|%3d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%9.5lf%9.5lf  ",a->a[k * cols + j][i].real(),a->a[k * cols + j][i].imag());
     }
    fprintf(file.out,"  |\n");
   }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

  if (remndr > 0) {
    fprintf(file.out,"| ");
    for (j=0;j<remndr;j++)
      fprintf(file.out,"  %3d  %8.3lf     ",k * cols + j + 1,*(e + k * cols + j) * scale_factor);
    for (j=0;j<cols - remndr;j++)
      fprintf(file.out,"                    ");
    fprintf(file.out,"    |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

      count = 1;
     for (i=0;i<a->iCols;i++){
      fprintf(file.out,"|%3d",count);
      count++;
       //fprintf(file.out,"|   ");
       for(j = 0; j < remndr; j++) 
         fprintf(file.out,"%9.5lf%9.5lf  ",a->a[k * cols + j][i].real(),a->a[k * cols + j][i].imag());
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"                    ");
       fprintf(file.out,"  |\n");
      }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

}

void print_real_eigenvector_matrix(DoubleMatrix *a, double *e, FILES file){

  int i, j, k, blocks, remndr;
  blocks = a->iCols / 8;
  remndr = a->iCols - blocks * 8;
  for (k=0;k<blocks;k++){
    for (j=0;j<8;j++){
      fprintf(file.out,"       %5d  ",k * 8 + j + 1);
     }
    fprintf(file.out,"\n");
       fprintf(file.out,"\n");
    for (j=0;j<8;j++){
      fprintf(file.out,"    %7.3e",*(e + k * 8 + j));
     }
    fprintf(file.out,"\n\n");
    for (i=0;i<a->iRows;i++){
      for (j=0;j<8;j++){
      fprintf(file.out,"       %6.3lf ",a->a[k * 8 + j][i]);
     }
    fprintf(file.out,"\n");
   }
      fprintf(file.out,"\n\n\n");
    }
    for (j=0;j<remndr;j++){
      fprintf(file.out,"       %5d  ",k * 8 + j + 1);
     }
       fprintf(file.out,"\n");
    for (j=0;j<remndr;j++){
      fprintf(file.out,"    %7.3e",*(e + k * 8 + j));
     }
    fprintf(file.out,"\n\n");
     for (i=0;i<a->iRows;i++){
       for(j = 0; j < remndr; j++) {
      fprintf(file.out,"       %6.3lf ",a->a[k * 8 + j][i]);
        }
       fprintf(file.out,"\n");
      }
     fprintf(file.out,"\n\n");

}

void print_real_eigenvector_matrix2(DoubleMatrix *a, double *e, int cols, int last_vector, double scale_factor, FILES file)

{

  int i, j, k, blocks, remndr, count;

  // print formatted matrix with cols columns. Print eigenvectors up to last_vector.

  blocks = a->iRows / cols;
  remndr = a->iRows - blocks * cols;
  if (last_vector < a->iRows) {
  blocks = last_vector / cols;
  remndr = last_vector - blocks * cols;
 }

  for (k=0;k<blocks;k++){
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|    ");
    for (j=0;j<cols;j++){
      fprintf(file.out," %3d%8.3lf",k * cols + j + 1,*(e + k * cols + j) * scale_factor);
     }
    fprintf(file.out,"     |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

      count = 1;
    for (i=0;i<a->iCols;i++){
      fprintf(file.out,"| %3d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%12.5lf",a->a[k * cols + j][i]);
     }
    fprintf(file.out,"     |\n");
   }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

  if (remndr > 0) {
    fprintf(file.out,"|    ");
    for (j=0;j<remndr;j++)
      fprintf(file.out," %3d%8.3lf",k * cols + j + 1,*(e + k * cols + j) * scale_factor);
    for (j=0;j<cols - remndr;j++)
      fprintf(file.out,"            ");
    fprintf(file.out,"     |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

      count = 1;
     for (i=0;i<a->iCols;i++){
      fprintf(file.out,"| %3d",count);
      count++;
       for(j = 0; j < remndr; j++) 
         fprintf(file.out,"%12.5lf",a->a[k * cols + j][i]);
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"            ");
       fprintf(file.out,"     |\n");
      }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

}

/*
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
      //fprintf(file.out,"%15.10lf  ",a->a[i][k * 6 + j].real());
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
         fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
        }
       fprintf(file.out,"\n");
      }
     fprintf(file.out,"\n\n");
    }

}
*/

void print_complex_matrix1(ComplexMatrix *a, int transpose, FILES file)

{

  int i, j, k, blocks, remndr;

  if (transpose == 0) {
  blocks = a->iCols / 6;
  remndr = a->iCols - blocks * 6;
 }

  if (transpose == 1) {
  blocks = a->iRows / 6;
  remndr = a->iRows - blocks * 6;
  fprintf(file.out,"Matrix Transpose requested\n\n");
 }

  for (k=0;k<blocks;k++){
    for (j=0;j<6;j++){
      fprintf(file.out,"     %5d        ",k * 6 + j + 1);
     }
  if (transpose == 0) {
    fprintf(file.out,"\n\n");
    for (i=0;i<a->iRows;i++){
      for (j=0;j<6;j++){
      fprintf(file.out,"%8.3lf%8.3lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
     }
    fprintf(file.out,"\n");
   }
  }
  if (transpose == 1) {
    fprintf(file.out,"\n\n");
    for (i=0;i<a->iCols;i++){
      for (j=0;j<6;j++){
      fprintf(file.out,"%8.3lf%8.3lf  ",a->a[k * 6 + j][i].real(),a->a[k * 6 + j][i].imag());
     }
    fprintf(file.out,"\n");
   }
  }
      fprintf(file.out,"\n\n\n");
    }

  if (remndr > 0) {
    for (j=0;j<remndr;j++){
      fprintf(file.out,"     %5d        ",k * 6 + j + 1);
     }
    fprintf(file.out,"\n\n");
  if (transpose == 0) {
     for (i=0;i<a->iRows;i++){
       for(j = 0; j < remndr; j++) {
         fprintf(file.out,"%8.3lf%8.3lf  ",a->a[i][k * 6 + j].real(),a->a[i][k * 6 + j].imag());
        }
       fprintf(file.out,"\n");
      }
     }
  if (transpose == 1) {
     for (i=0;i<a->iCols;i++){
       for(j = 0; j < remndr; j++) {
         fprintf(file.out,"%8.3lf%8.3lf  ",a->a[k * 6 + j][i].real(),a->a[k * 6 + j][i].imag());
        }
       fprintf(file.out,"\n");
      }
     }

     fprintf(file.out,"\n\n");
    }

}

void print_complex_matrix2(ComplexMatrix *a, int transpose, int cols, double scale_factor, FILES file)

{

  int i, j, k, blocks, remndr, count;

  // print formatted matrix with cols columns. transpose == 0 => not transposed, transpose == 1 => transposed
  // matrix elements are scaled by scale_factor, e.g. au_to_eV

  if (transpose == 0) {
  blocks = a->iCols / cols;
  remndr = a->iCols - blocks * cols;
fprintf(file.out,"%3d %3d\n",blocks,remndr);
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
      fprintf(file.out,"       %5d        ",k * cols + j + 1);
     }
    fprintf(file.out,"   |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  if (transpose == 0) {
    count = 1;
    for (i=0;i<a->iRows;i++){
      fprintf(file.out,"|%3d",count);
      count++;
      for (j=0;j<cols;j++){
      fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * cols + j].real() * scale_factor,a->a[i][k * cols + j].imag() * scale_factor);
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
      fprintf(file.out,"%9.5lf%9.5lf  ",a->a[k * cols + j][i].real() * scale_factor,a->a[k * cols + j][i].imag() * scale_factor);
     }
    fprintf(file.out,"  |\n");
   }
  }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    }

  if (remndr > 0) {
    fprintf(file.out,"|  ");
    for (j=0;j<remndr;j++)
      fprintf(file.out,"       %5d        ",k * cols + j + 1);
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
         fprintf(file.out,"%9.5lf%9.5lf  ",a->a[i][k * cols + j].real() * scale_factor,a->a[i][k * cols + j].imag() * scale_factor);
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
         fprintf(file.out,"%9.5lf%9.5lf  ",a->a[k * cols + j][i].real() * scale_factor,a->a[k * cols + j][i].imag() * scale_factor);
       for (j = 0; j < cols - remndr; j++)
         fprintf(file.out,"                    ");
       fprintf(file.out,"  |\n");
      }
     }

    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
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

void PrintComplexMatrix(ComplexMatrix *a, FILES file){
  int i,j;
  for (i=0;i<a->iRows;i++){
    for (j=0;j<a->iCols;j++){
      //std::printf("(%16.8f,%16.8f)",a->a[i][j].real(),a->a[i][j].imag());
      fprintf(file.out,"%5.2f %5.2f  ",a->a[i][j].real(),a->a[i][j].imag());
    }
    fprintf(file.out,"\n");
    //std::printf("\n");
  }
    fprintf(file.out,"\n");
}

