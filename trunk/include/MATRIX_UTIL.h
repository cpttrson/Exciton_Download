
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void double_mat_dot(double*, double*, double*);

void double_vec_cross(VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

double double_vec_dot(VECTOR_DOUBLE*, VECTOR_DOUBLE*);

void double_vec_diff(VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

int check_inverse1(double*, double*);

void map_to_wigner(CRYSTAL*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

void map_to_brillouin(CRYSTAL*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_INT*, FILES);

void getcart(VECTOR_INT*, VECTOR_DOUBLE*, int*, CRYSTAL*);

void rotate_vector3(double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

void rotate_vector_int(int*, VECTOR_INT*, VECTOR_INT*);

void rotate_vector_int_latt(int*, int[3], int, VECTOR_INT*, VECTOR_INT*, CRYSTAL*, JOB_PARAM*, FILES);

int check_vec(VECTOR_DOUBLE*, VECTOR_DOUBLE*);

int check_mat(double*, double*) ;

