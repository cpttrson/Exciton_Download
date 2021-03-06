
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void atom_scf(ATOM *, int, DoubleMatrix *, double *, double *, int *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);

void scf_atom(ATOM *, int, DoubleMatrix *, double *, double *, int *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);

void init_mix(DoubleMatrix *, DoubleMatrix *, DoubleMatrix *, DoubleMatrix *, DoubleMatrix *);

void mix_density(DoubleMatrix *,DoubleMatrix *,DoubleMatrix *,DoubleMatrix *,double *,double *,double *,int,int*,JOB_PARAM*, FILES);

void initial_density_matrix_atom(DoubleMatrix *, ATOM *, int *, double *, int *, SHELL *, JOB_PARAM *, FILES);

void fock_element_1e_atm(DoubleMatrix*, DoubleMatrix*, int*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void integrals_2e_atom(INTEGRAL_LIST *, int, ATOM *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);

void fock_2e_matrix(DoubleMatrix*, INTEGRAL_LIST*, DoubleMatrix*, int, ATOM*, JOB_PARAM*, FILES);

void total_energy_atom(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, int, double*, double* , ATOM*, JOB_PARAM*, FILES);

void final_total_energy_atom(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, int, double*, double* , ATOM*, JOB_PARAM*, FILES);

double e_atom(int , int , int , double );

double ftuvn_atom(int , int , int , int , double *);

void g000m(int , double *, double , int );
