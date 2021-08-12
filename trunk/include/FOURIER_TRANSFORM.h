
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void fourier_transform(double*, Complex*, KPOINT_TRAN*, int*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_3(double*, Complex*, KPOINT_TRAN*, int*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*,FILES);

void fourier_transform_molecule(double*, double*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_tensor_3(double*, double*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_salc(double*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_salc1(double*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_salc_crystal(double*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, REAL_LATTICE_TABLES*, ATOM*, SHELL*, SYMMETRY*, SYMMETRY*, JOB_PARAM*, FILES);

void count_atom_cluster(int, int, SALC*, int*, int*, int*, SYMMETRY*, ATOM_TRAN*, REAL_LATTICE_TABLES*);

void generate_atom_cluster(IntMatrix*, SALC*, int*, int*, int*, SYMMETRY*);

