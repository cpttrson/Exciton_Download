
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void integrals_molecule_ij(double*, PAIR_TRAN*, ATOM*, SHELL*, GAUSSIAN*, REAL_LATTICE*, JOB_PARAM*, FILES);

void nuclear_repulsion_energy(double*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void fock_element_1e1(INT_1E*, int, PAIR_TRAN*, int, int[], REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void fock_element_1e2(INT_1E*, PAIR_TRAN*, int[], REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void fock_element_elecnuc(double*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void Kinetic(double*, int, int, int, int, int, int, int, int, int,  double*, double*, double*, double*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void ElecNuc(double*,int,int,int,int,int,int,int,double*,double*,double*,VECTOR_DOUBLE*,double*,double*,REAL_LATTICE*,RECIPROCAL_LATTICE*,CRYSTAL*,ATOM*,SHELL*,JOB_PARAM*,FILES);

void Momentum(double*, int, int, int, int, int, int, int, int, int, int,  double*, double*, double*, double*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void Overlap(double*, int, int, int, int, int, int, int, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

void Dipole(double*,int,int,int,int,int,int,int,int,int,int,int,double*,double*,double*,double*,VECTOR_DOUBLE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);

void two_centre_coulomb(double*,int,int,int,int,int,int,int,int,int,double*,double*,double*,VECTOR_DOUBLE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);

