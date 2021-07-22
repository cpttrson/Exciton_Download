

  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

/*
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <fstream>
*/
#include <cstdlib>
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
  fprintf(file.out,"|        CONTRIBUTORS                                                                                     |\n");
  fprintf(file.out,"|        Alin-Marin Elena                                                                                 |\n");
  fprintf(file.out,"|        Svjetlana Gamaic-Mulaomerovic                                                                    |\n");
  fprintf(file.out,"|                                                                                                         |\n");
  fprintf(file.out,"|        School of Physics                                                                                |\n");
  fprintf(file.out,"|        University of Dublin, Trinity College, Dublin 2, Ireland                                         |\n");
  fprintf(file.out,"|        Last modified July 22nd 2021                                                                    |\n");
  fprintf(file.out,"|                                                                                                         |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");

  return;
}
