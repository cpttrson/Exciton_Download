
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void getangle(double x, double y, double z, double a, double b);

int Lebedev_Laikov_npoint(int, FILES);

int Lebedev_Laikov_Oh (int n, double a, double b, double v, double *x, double *y, double *z, double *w);

int Lebedev_Laikov_sphere (int N, double *X, double *Y, double *Z, double *W);

int Lebedev_grid(double *, double *, double *, double *, double *, double *, int *, int, FILES);
