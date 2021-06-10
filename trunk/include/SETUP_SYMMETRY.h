#ifndef SETUP_SYMMETRYH
#define SETUP_SYMMETRYH

void read_symmetry_group(CRYSTAL*, SYMMETRY*, int, JOB_PARAM*, FILES);

void generate_cartesian_symmetry_operators(CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void reorder_symmetry_operators_by_class(SYMMETRY*, JOB_PARAM*, FILES);

void generate_operator_inverses(SYMMETRY*, JOB_PARAM*, FILES);

void generate_group_multiplication_table(SYMMETRY*, int, FILES);

void generate_group_conjugacy_classes(SYMMETRY*, int, FILES);

void generate_Dirac_characters(SYMMETRY*, JOB_PARAM*, FILES);

void generate_permutation_group_table(SYMMETRY*, JOB_PARAM*, FILES);

void print_symmetry_operators(SYMMETRY*, JOB_PARAM*, FILES);

void count_little_k_group_operators(int, SYMMETRY*, SYMMETRY*, CRYSTAL*, KPOINT_TRAN*, FERMI*, JOB_PARAM*, FILES);
//void count_little_k_group_operators(int, SYMMETRY*, SYMMETRY*, KPOINT_TRAN*, FERMI*, JOB_PARAM*, FILES);

void generate_little_k_group_operstors(int, SYMMETRY*, SYMMETRY*, CRYSTAL*, KPOINT_TRAN*, FERMI*, JOB_PARAM*, FILES);
//void generate_little_k_group_operstors(int, SYMMETRY*, SYMMETRY*, KPOINT_TRAN*, FERMI*, JOB_PARAM*, FILES);

void generate_little_k_group(int, SYMMETRY*, FERMI*, KPOINT_TRAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);
//void generate_little_k_group(int, SYMMETRY*, FERMI*, KPOINT_TRAN*, SYMMETRY*, JOB_PARAM*, FILES);

#endif
