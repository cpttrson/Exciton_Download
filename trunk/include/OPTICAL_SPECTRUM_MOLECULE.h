
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void optical_spectrum_molecule(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void electroabsorption_spectrum_molecule(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void dipole_matrix_elements_molecule(DoubleMatrix*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES* , RECIPROCAL_LATTICE* , JOB_PARAM* , FILES);

void fragment_transition_dipole_moment_molecule(FERMI*, int*, int*, int[][2], DoubleMatrix*, double*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

