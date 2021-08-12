
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

double e(int, int, int, double, double, double);

double ftuvn(int, int, int, int, double *, VECTOR_DOUBLE);

double cosfactor(int, double);

Complex cosfactor_complex(int, double);

void non_recursive_ftuvn(int, int, double f[][13][13][13], double en[][55], VECTOR_DOUBLE *);

void erf_derivative(double*, int, double, double, RECIPROCAL_LATTICE*);

void erf_exp_derivative(double*, double*, int, double, double, RECIPROCAL_LATTICE*);

void erfc_derivative(double*, int, double, double, RECIPROCAL_LATTICE*);

