#ifndef KPOINTSH
#define KPOINTSH

void knet_size(int *, int is[3], CRYSTAL *);

void count_k_points(KPOINT_TRAN*, int is[3], CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_k_points(KPOINT_TRAN*, int[3],  CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_k_point_pairs(KPOINT_TRAN*, int[3], CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void k_point_equivalent_pairs(int, KPOINT_PAIR_TRAN*, KPOINT_TRAN*, int[3], CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void q_point_equivalent_pairs(int, KPOINT_PAIR_TRAN*, KPOINT_TRAN*, int[3], CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void q_point_rotate_pairs(int, int*, int*, int*, int*, FERMI*, KPOINT_TRAN*, KPOINT_TRAN*, int[3], CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void print_knet(KPOINT_TRAN*, int is[3], CRYSTAL*, JOB_PARAM*, FILES);

void tetrahedron_vertices2(int [3], int, int k[10][4], KPOINT_TRAN *, CRYSTAL *, JOB_PARAM *, FILES);

void polygons_vertices_vol(int [], int *, int *, double *, CRYSTAL *, JOB_PARAM *, FILES file);

void polygons_vertices_state_density(int [], int *, int *, double *, CRYSTAL *, JOB_PARAM *, FILES file);

void polygons_vertices_susceptibility(int [], int *, int *, double *, CRYSTAL *, JOB_PARAM *, FILES file);

VECTOR_INT decompose_k_point(int [3], int, CRYSTAL *, JOB_PARAM *, FILES);

int compose_k_point(int [3], VECTOR_INT, VECTOR_INT, CRYSTAL *, JOB_PARAM *, FILES);

VECTOR_INT add_k_points(VECTOR_INT, VECTOR_INT);

VECTOR_INT subtract_k_points(VECTOR_INT, VECTOR_INT);

void occupation_fermi(FERMI *, KPOINT_TRAN*, double *, int, int, JOB_PARAM *, FILES);
//CHANGES2014void occupation_fermi(FERMI *, double *, int, int, JOB_PARAM *, FILES);

void calculate_fermi_level(FERMI *, double *, int *, KPOINT_TRAN *, int, int, ATOM *, JOB_PARAM *, FILES);

void calculate_fermi_level_metal(FERMI *, double *, int *, KPOINT_TRAN *, int, int, ATOM *, JOB_PARAM *, FILES);

void calculate_fermi_level_metal_old(FERMI *, double *, int *, KPOINT_TRAN *, int, int, ATOM *, JOB_PARAM *, FILES);

#endif
