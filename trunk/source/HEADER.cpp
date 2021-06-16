/*
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <fstream>
*/
#include "USER_DATA.h"
#include "HEADER.h"

using namespace std;

void print_header(FILES file)

{

  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                                                                                                         |\n");
  fprintf(file.out,"|        E X C I T O N 21                                                                                 |\n");
  fprintf(file.out,"|        Revision 232                                                                                     |\n");
  fprintf(file.out,"|                                                                                                         |\n");
  fprintf(file.out,"|        AUTHORS                                                                                          |\n");
  fprintf(file.out,"|        Charles H. Patterson                                                                             |\n");
  fprintf(file.out,"|                                                                                                         |\n");
  fprintf(file.out,"|        School of Physics                                                                                |\n");
  fprintf(file.out,"|        University of Dublin, Trinity College, Dublin 2, Ireland                                         |\n");
  fprintf(file.out,"|        Last modified April 30th 2021                                                                    |\n");
  fprintf(file.out,"|                                                                                                         |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");

  return;
}
