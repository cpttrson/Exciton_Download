
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

int init_Complex_array(int, Complex*);

void eigenvec_plot(int*, int*, int, KPOINT_TRAN*, ComplexMatrix*, VECTOR_INT*, int*, VECTOR_DOUBLE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void eigenvec_isosurface(int*, int*, int, KPOINT_TRAN*, VECTOR_INT*, int*, VECTOR_DOUBLE*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, SYMMETRY*, JOB_PARAM*,FILES);

void density_isosurface(double*, int*, int*, int*, VECTOR_DOUBLE*, KPOINT_TRAN*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, REAL_LATTICE_TABLES*, SYMMETRY*, JOB_PARAM*,FILES);
void density_isosurface1(double*, int*, int*, int*, VECTOR_DOUBLE*, KPOINT_TRAN*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, REAL_LATTICE_TABLES*, SYMMETRY*, JOB_PARAM*,FILES);

void band_plot(int*, int*, int, int, double[][5], KPOINT_TRAN*, VECTOR_INT*, VECTOR_INT*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*,  REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void bulk_band_projection(int*, int*, int*, KPOINT_TRAN*, double*, ATOM_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void plot_points(int, int*, KPOINT_TRAN*, int*, VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, Complex*, ComplexMatrix*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void plot_points1(double*, VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, int, double*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void plot_correlated_electron_hole(int[3], VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, int[2], FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void plot_bse_transition_density(int[3], VECTOR_DOUBLE*, VECTOR_DOUBLE*, int[2], FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void plot_bse_transition_density_crystal(double*, int[2],  int*, VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void plot_correlated_electron_hole_crystal(double*, int[2],  int*, VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void plot_electron_hole_transition_density_crystal(double*, int[2],  int*, VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void plot_electron_hole_molecule(int[3], VECTOR_DOUBLE*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void wavefunction_gridpoint(double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, double*, int, int, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void wavefunction_gridpoint_crystal(Complex*, VECTOR_DOUBLE*, int, VECTOR_DOUBLE*, int, KPOINT_TRAN*, Complex*, int, int, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void write_isosurface_xsf(char*, double*, int, int[3], VECTOR_DOUBLE*, int, int, CRYSTAL*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void calc_grid(int *, VECTOR_DOUBLE *, int, int *, VECTOR_DOUBLE *, JOB_PARAM *, FILES);

void calc_grid1(int *, VECTOR_DOUBLE *, int, int *, VECTOR_DOUBLE *, JOB_PARAM *, FILES);

void Rvec_grid_size(int *, int, CRYSTAL *, FILES);

void generate_Rvec_grid(int, int, VECTOR_DOUBLE *, CRYSTAL *, FILES);

