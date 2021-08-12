
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void allocate_JOB_PARAM(JOB_PARAM*, FILES);

void free_JOB_PARAM(JOB_PARAM*, FILES);

void allocate_SYMMETRY(SYMMETRY*, JOB_PARAM*, FILES);

void free_SYMMETRY(SYMMETRY*, JOB_PARAM*);

void allocate_REAL_LATTICE(REAL_LATTICE*, JOB_PARAM*, FILES);

void free_REAL_LATTICE(REAL_LATTICE*, JOB_PARAM*);

void allocate_ATOM(ATOM*, JOB_PARAM*, FILES);

void free_ATOM(ATOM*, JOB_PARAM*);

void allocate_ATOM_TRAN(ATOM_TRAN*, ATOM*, SYMMETRY*, JOB_PARAM*, FILES);

void free_ATOM_TRAN(ATOM_TRAN*, JOB_PARAM*);

void allocate_REAL_LATTICE_TABLES(REAL_LATTICE_TABLES*, SYMMETRY*, JOB_PARAM*, FILES);

void free_REAL_LATTICE_TABLES(REAL_LATTICE_TABLES*, JOB_PARAM*);

void allocate_SHELL_GAUSSIAN(SHELL*, GAUSSIAN*, ATOM*, JOB_PARAM*, FILES);

void free_SHELL_GAUSSIAN(SHELL*, GAUSSIAN*, JOB_PARAM*);

void allocate_SALC(SALC*, SYMMETRY*, JOB_PARAM*, FILES);

void free_SALC(SALC*, JOB_PARAM*);

void allocate_PAIR_TRAN(PAIR_TRAN*, ATOM*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void free_PAIR_TRAN(PAIR_TRAN*, JOB_PARAM*);

void allocate_TRIPLE_TRAN(TRIPLE_TRAN*, JOB_PARAM*, FILES);

void free_TRIPLE_TRAN(TRIPLE_TRAN*, JOB_PARAM*);

void allocate_QUAD_TRAN(QUAD_TRAN*, JOB_PARAM*, FILES);

void free_QUAD_TRAN(QUAD_TRAN*, JOB_PARAM*);

void allocate_k_points(KPOINT_TRAN*, CRYSTAL*, JOB_PARAM*, FILES);

void free_k_points(KPOINT_TRAN*, JOB_PARAM*);

/*
void allocate_kq_pairs(KQPOINT_TRAN*, FERMI*, CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void free_kq_pairs(KQPOINT_TRAN*, JOB_PARAM*);

void allocate_k_pairs(KPOINT_PAIR_TRAN*, int[3], CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void free_k_pairs(KPOINT_PAIR_TRAN*, JOB_PARAM*);
*/

void allocate_fermi(FERMI*, ATOM*, JOB_PARAM*, FILES);

void free_fermi(FERMI*, JOB_PARAM*);

void allocate_INT_1E(INT_1E*, int, int[], JOB_PARAM *, FILES);

void free_INT_1E(INT_1E*, int[], JOB_PARAM*, FILES);

void allocate_integral_list(INTEGRAL_LIST*, int, JOB_PARAM*, FILES);

void allocate_integral_list_complex(INTEGRAL_LIST_COMPLEX*, int, JOB_PARAM*, FILES);

void free_integral_list(INTEGRAL_LIST*, JOB_PARAM*);

void free_integral_list_complex(INTEGRAL_LIST_COMPLEX*, JOB_PARAM*);

void allocate_dft_grid(DFT_GRID*, int, ATOM*, JOB_PARAM*, FILES);

void free_dft_grid(DFT_GRID*, JOB_PARAM*);

void DestroyIntArray(int**, int*, JOB_PARAM*);

void AllocateIntArray(int**, int*, JOB_PARAM*);

void ResetIntArray(int*, int*);

void AllocateDoubleArray(double**, int*, JOB_PARAM*);

void ResetDoubleArray(double*, int*);

void DestroyDoubleArray(double**, int*, JOB_PARAM*);

void AllocateComplexArray(Complex**, int*, JOB_PARAM*);

void ResetComplexArray(Complex*, int*);

void DestroyComplexArray(Complex**, int*, JOB_PARAM*);

void AllocateIntMatrix(IntMatrix**, int*, int*, JOB_PARAM*);

void ResetIntMatrix(IntMatrix*);

void DestroyIntMatrix(IntMatrix**, JOB_PARAM*);

void AllocateDoubleMatrix(DoubleMatrix**, int*, int*, JOB_PARAM*);

void ResetDoubleMatrix(DoubleMatrix*);

void DestroyDoubleMatrix(DoubleMatrix**, JOB_PARAM*);

void AllocateComplexMatrix(ComplexMatrix**, int*, int*, JOB_PARAM*);

void ResetComplexMatrix(ComplexMatrix*);

void DestroyComplexMatrix(ComplexMatrix**, JOB_PARAM*);

void array_dimensions(int*, int*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void sh_array_dimensions(int*, int*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void triple_array_dimensions(int*, int*, TRIPLE_TRAN*, ATOM*, JOB_PARAM*, FILES);

void sh_triple_array_dimensions(int*, int*, TRIPLE_TRAN*, ATOM*, JOB_PARAM*, FILES);

