
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

extern void print_error(const char *str, BOOLEAN label);

extern void print_time(char printing, char *a);

extern BOOLEAN search_file(FILE *fp, const char *keyword);

extern BOOLEAN read_line(FILE * ffile, char * line, int MAXIMUM_LEN_FILE_LINE);

