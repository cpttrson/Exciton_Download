#ifndef SETUP_RECIPROCAL_LATTICEH
#define SETUP_RECIPROCAL_LATTICEH

void count_reciprocal_lattice_vectors(CRYSTAL, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void generate_reciprocal_lattice(CRYSTAL, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

//void generate_reciprocal_lattice(CRYSTAL, RECIPROCAL_LATTICE*, SYMMETRY, FILES);

void generate_q_lattice(VECTOR_INT*, Q_LATTICE*, FERMI*, RECIPROCAL_LATTICE*, CRYSTAL*, JOB_PARAM*, FILES);

#endif
