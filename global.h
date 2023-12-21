#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <stdbool.h>
#include <lapacke.h> // Includes <complex.h>. LAPACKE is the C wrapper for the standard F90 LAPACK library.
#include <time.h>
#include <unistd.h>  // For getpid()
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include "mt19937ar.h"

#define EPSILON 1e-10

/* global variables */

extern int d; /* dimension */
extern int L; /* no. of states in the ensemble */
extern int N0; /* sample complexity/no. of measurements */
extern const double zero_epsilon;


// Pauli X
lapack_complex_double X[4];

// Pauli Y
lapack_complex_double Y[4];

// Pauli Z
lapack_complex_double Z[4];

// 2-design elements
lapack_complex_double U_0_0[4];
lapack_complex_double U_pi_0[4];
lapack_complex_double U_pi_2_0[4];
lapack_complex_double U_pi_2_pi[4];
lapack_complex_double U_pi_2_pi_2[4];
lapack_complex_double U_pi_2_3pi_2[4];

// Identity
lapack_complex_double Id[4];

// H
lapack_complex_double H[4];

// S
lapack_complex_double S[4];

// HS
lapack_complex_double HS[4];

// SH
lapack_complex_double SH[4];

// SS(Z)
lapack_complex_double SS[4];

// HSH
lapack_complex_double HSH[4];

// HSS (HZ)
lapack_complex_double HSS[4];

// SHS
lapack_complex_double SHS[4];

// SSH
lapack_complex_double SSH[4];

// SSS (Sd)
lapack_complex_double SSS[4];

// HSHS
lapack_complex_double HSHS[4];

// HSSH
lapack_complex_double HSSH[4];

// HSSS
lapack_complex_double HSSS[4];

// SHSS
lapack_complex_double SHSS[4];

// SSHS
lapack_complex_double SSHS[4];

// HSHSS
lapack_complex_double HSHSS[4];

// HSSHS
lapack_complex_double HSSHS[4];

// SHSSH
lapack_complex_double SHSSH[4];

// SHSSS
lapack_complex_double SHSSS[4];

// SSHSS
lapack_complex_double SSHSS[4];

// HSHSSH
lapack_complex_double HSHSSH[4];

// HSHSSS
lapack_complex_double HSHSSS[4];

// HSSHSS
lapack_complex_double HSSHSS[4];
/*=============================================================================*/


double prior_constant, average_risk_prev;
lapack_complex_double one, zero, i1;

int *outcome_vector_empirical;
double *log_prior, *log_prior_Ginibre;
lapack_complex_double *id, *basis, ***rho, ***rho_mixed;

//lapack_complex_double ***Bayes_est_empirical, **Bayes_est_empirical_sum;

struct risk_types{
    double risk, average_risk;
};

/* functions */
double min(double A[],int size);


double max(double A[],int size);


void convert_pmf_to_cmf(double *pmf, int count);


int is_normalized(double *prob_dist, int size);


void print_matrix(lapack_complex_double* matrix, int size);


lapack_complex_double* random_clifford_qubit();


lapack_complex_double* random_2Design_qubit();


lapack_complex_double* random_1Design_qubit();


int* from_dAlphabet_to_dBase(int *res, int base, int exponent, int inputNum);


int func_inverse_sampling_quantum_states(double *prior);


double func_rel_ent(lapack_complex_double *rho, lapack_complex_double *Bayes_est_copy);


lapack_complex_double* func_Haar_random_B();


lapack_complex_double* func_Ginibre_random_rho();


lapack_complex_double* func_Haar_random_state();


int check_psd_trace_one(lapack_complex_double* mat, int n);


lapack_complex_double* sqrtm_LAPACK(int d, lapack_complex_double* A);


lapack_complex_double* sqrtm_LAPACK_Schur(int d, lapack_complex_double* A);


lapack_complex_double* multiply_three_matrices_LAPACK(int d, lapack_complex_double* A, lapack_complex_double* B, lapack_complex_double* C);


lapack_complex_double trace_LAPACK(int d, lapack_complex_double* A);


double trace_product_LAPACK(int n, lapack_complex_double* A, lapack_complex_double* B);



bool is_normal(lapack_complex_double* A, int d);


double fidelity_mixed(int d, lapack_complex_double* rho1, lapack_complex_double* rho2);


double fidelity_mixed_simple(int d, lapack_complex_double* rho1, lapack_complex_double* rho2);


double fidelity_pure(int d, lapack_complex_double* psi1, lapack_complex_double* psi2);


lapack_complex_double* func_MUBs_prime_powers();


struct risk_types func_Bayes_update_blind(int n);


struct risk_types func_Bayes_update_blind_Ginibre(int n);


double log_sum_exp(double *log_probs, int n);


void normalize_log_probs(double *log_probs, int n);


double fidelity_pure_state_density_matrix(int n, lapack_complex_double* psi, lapack_complex_double* rho);


struct risk_types func_Bayes_update_empirical_average(int x, int n);


double complex* sqrtm(double complex* A);


double complex* pseudoinv(size_t n, double complex* A);


double complex* multiply_two_matrices(double complex* A, double complex* B);


double complex* multiply_three_matrices(int n, double complex* A, double complex* B, double complex* C);


double complex**  generate_random_unitary_matrices(int num);


double complex* generate_random_unitary_matrix();


double complex* generate_random_complex_matrix();


double complex* apply_unital_channel(int L, int d, double complex** U, double complex* A);


double trace_product(int n, double complex* A, double complex* B);


double complex* compute_Bayes_mean_sum(int x, double complex* M_x, double complex** rho, double* p);


void generate_random_density_matrices(int L, double complex** rho);
double complex trace(double complex* A);


double fidelity(int d, double complex* rho1, double complex* rho2);


double* compute_fidelities(double complex* A, double complex** B_x);


double complex** compute_p_rho(int num_matrices, int n, double* p, double complex** rho);


double complex * compute_sum_p_rho(int num_matrices, int n, double* p, double complex **rho);
