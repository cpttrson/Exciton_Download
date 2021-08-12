
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void read_symmetry_group(CRYSTAL *crystal, SYMMETRY *symmetry, int, JOB_PARAM *job, FILES file);

void conventional_unit_cell(CRYSTAL *crystal, FILES file);

void primitive_unit_cell(CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_cartesian_symmetry_operators(CRYSTAL *crystal, SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

void count_real_lattice_vectors(CRYSTAL crystal, REAL_LATTICE *R, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES file);

void generate_real_lattice(CRYSTAL crystal, SYMMETRY s, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

void generate_real_lattice1(CRYSTAL crystal, SYMMETRY s, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

void generate_real_lattice2(CRYSTAL crystal, SYMMETRY s, REAL_LATTICE *R, RECIPROCAL_LATTICE*, JOB_PARAM *job, FILES file);

void count_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, JOB_PARAM*, FILES);

void select_pairs4(PAIR_TRAN*, ATOM*, REAL_LATTICE*, JOB_PARAM*, FILES);

void generate_real_lattice_shells(CRYSTAL crystal, SYMMETRY s, REAL_LATTICE *R, JOB_PARAM *job, FILES file);

void generate_real_lattice_tables(CRYSTAL crystal, SYMMETRY s, REAL_LATTICE *, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

void count_reciprocal_lattice_vectors(CRYSTAL crystal, RECIPROCAL_LATTICE *G, JOB_PARAM*, FILES file);

void generate_reciprocal_lattice(CRYSTAL crystal, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file);

void generate_reciprocal_lattice(CRYSTAL crystal, RECIPROCAL_LATTICE *G, SYMMETRY symmetry, FILES file);

void count_unique_atoms(ATOM *atoms, JOB_PARAM *job, FILES file);

void read_unique_atoms(ATOM *atoms, ATOMTYPE *iratom, CRYSTAL crystal, JOB_PARAM *job, FILES file);

void count_basis_sets(CRYSTAL crystal, ATOM *atoms, JOB_PARAM *job, FILES file);

void read_basis_sets(CRYSTAL crystal, ATOM *atoms, ATOMTYPE *iratom, JOB_PARAM *job, FILES file);

void count_all_atoms(CRYSTAL crystal, ATOM *atoms, ATOMTYPE *iratom, SYMMETRY symmetry, JOB_PARAM *job, FILES file);

void count_all_crystal_atoms(CRYSTAL, ATOM*, ATOMTYPE*, SYMMETRY, JOB_PARAM*, FILES);

void generate_all_crystal_atoms(CRYSTAL , REAL_LATTICE , SYMMETRY *, ATOMTYPE *, ATOMTYPE *, ATOM *, ATOM_TRAN *, ATOM_TRAN *, JOB_PARAM *, FILES );

void generate_all_atoms(CRYSTAL , REAL_LATTICE , SYMMETRY , ATOMTYPE *, ATOMTYPE *, ATOM *, ATOM_TRAN *, ATOM_TRAN *, JOB_PARAM *, FILES );

void generate_runtime_arrays(ATOM *, ATOM_TRAN *, ATOMTYPE *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES );

void generate_rotation_operators(SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

void map_to_wigner1(CRYSTAL *crystal, VECTOR_DOUBLE *vec, VECTOR_DOUBLE *mapped_vec, VECTOR_DOUBLE *map_vec);

void open_files(int argc, char *argv[], int rank, FILES file);

void allocate_CRYSTAL(CRYSTAL crystal);

void allocate_SYMMETRY(SYMMETRY *symmetry, FILES file);

void free_SYMMETRY(SYMMETRY *symmetry);

void allocate_ATOM(ATOM *atoms, FILES file);

void free_ATOM(ATOM *atoms);

void allocate_ATOM_TRAN(ATOM_TRAN *atom_p, ATOM atoms, SYMMETRY symmetry, FILES file);

void free_ATOM_TRAN(ATOM_TRAN *atom_p);

void allocate_SHELL_GAUSSIAN(SHELL *shell, GAUSSIAN *gaussian, ATOM atoms, FILES file);

void free_SHELL_GAUSSIAN(SHELL *shell, GAUSSIAN *gaussian);

void allocate_REAL_LATTICE(REAL_LATTICE *R, REAL_LATTICE_TABLES R_tables, FILES file);

void allocate_REAL_LATTICE1(REAL_LATTICE *R, FILES file);

void allocate_REAL_LATTICE2(REAL_LATTICE *R, FILES file);

void free_REAL_LATTICE(REAL_LATTICE *R);

void allocate_RECIPROCAL_LATTICE(RECIPROCAL_LATTICE *G, FILES file);

void free_RECIPROCAL_LATTICE(RECIPROCAL_LATTICE *G);

void allocate_REAL_LATTICE_TABLES(REAL_LATTICE_TABLES *R_tables, SYMMETRY symmetry, FILES file);

void allocate_REAL_LATTICE_TABLES1(REAL_LATTICE_TABLES *, SYMMETRY, FILES);

void free_REAL_LATTICE_TABLES(REAL_LATTICE_TABLES *R_tables);

void allocate_pairs(PAIR_TRAN *pair_p, int pair_count, int max_pairs, SYMMETRY *symmetry, FILES file);

void allocate_pairs1(PAIR_TRAN*, ATOM*, SYMMETRY*, REAL_LATTICE_TABLES*, FILES);

void allocate_density_pairs(PAIR_TRAN *pair_p, int pair_count, int max_pairs, SYMMETRY *symmetry, ATOM *atoms, REAL_LATTICE_TABLES *, FILES file);

void allocate_density_pairs1(PAIR_TRAN*, ATOM*, SYMMETRY *symmetry, REAL_LATTICE_TABLES*, FILES);

void allocate_quads(PAIR_TRAN*, FILES);

void free_quads(PAIR_TRAN*);

void allocate_triples1(PAIR_TRAN*, ATOM*, FILES);

void allocate_t_quads1(PAIR_TRAN*, ATOM*, FILES);

void free_triples1(PAIR_TRAN*);

void free_pairs(PAIR_TRAN *pair_p);

void free_pairs1(PAIR_TRAN *pair_p);

void free_t_quads1(PAIR_TRAN*);

void free_density_pairs(PAIR_TRAN *pair_p);

void free_density_pairs1(PAIR_TRAN *pair_p);

void allocate_integral_list(INTEGRAL_LIST *integral_list, int num, FILES file);

void free_integral_list(INTEGRAL_LIST *integral_list);

void allocate_JOB_PARAM(JOB_PARAM *, FILES);

void free_JOB_PARAM(JOB_PARAM *, FILES);

void allocate_fermi(FERMI *, int *, int *, ATOM *, JOB_PARAM *, FILES);

void free_fermi(FERMI *);

