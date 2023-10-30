//=============================================================
// Copyright 2023 Maria Quadeer
// Distributed under the GNU General Public License, Version 3.
// (See accompanying file LICENSE.txt or copy at
// https://www.gnu.org/licenses/gpl-3.0.en.html)
//=============================================================

#include "global.h"

/* Define extern variables. */
int d; /* dimension */
int L; /* no. of states in the ensemble */
const double zero_epsilon = 1e-6; /* Defining 0 for comparison with float types. */;

int main(int argc, char* argv[])
{

  if (argc < 3) {
        fprintf(stderr, "Usage: %s <d> <L>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    d = atoi(argv[1]); /* dimension */
    L = atoi(argv[2]); /* no. of states in the ensemble */

  // Seed the random number generator
  unsigned long seed;
  seed = 39845;
  init_genrand(seed);

  //Define global constants
  one = lapack_make_complex_double(1.0,0.0);
  zero = lapack_make_complex_double(0.0,0.0);

  // Declare and define a set of random density matrices, prior_constant.
  prior_constant = (double*) calloc(L, sizeof(double));
  double complex** rho = (double complex**) calloc(L, sizeof(double complex*));

  for(int l = 0; l < L; l++){
    rho[l] = (double complex*) calloc(d*d, sizeof(double complex));
    prior_constant[l]=1.0/L;
  }

  generate_random_density_matrices(L, rho);

  // Testing output of generate_random_density_matrices
    printf("Printing L random density matrices\n");
    for (int i = 0; i < L; i++) {
      printf("i=%d\n", i);
      for (int j = 0; j < d*d; j++) {
        printf("(%lf, %lf)\n", creal(rho[i][j]), cimag(rho[i][j]));
      }
    }

  double complex** E_pgm = (double complex**) calloc(L, sizeof(double complex*));
  double complex** U_E_pgm = (double complex**) calloc(L, sizeof(double complex*));
  double complex** bme = (double complex**) calloc(L, sizeof(double complex*));
  double complex** bme_U = (double complex**) calloc(L, sizeof(double complex*));

  // Compute \sum_j p_j \rho_j
  double complex* sum;
  sum = compute_sum_p_rho(L, d, prior_constant, rho);

  // Compute pseudo-inverse of sum and then square root.
  double complex *sum_pinv, *sum_pinv_sqrt;
  sum_pinv = pseudoinv(d, sum);
  sum_pinv_sqrt = sqrtm(sum);

  // Evaluate the set of scaled density matrices \{p_j\rho_j\}.
  double complex** p_rho;
  p_rho = compute_p_rho(L, d, prior_constant, rho);


  // Generate a set of L random unitary matrices.
  double complex** U = (double complex**)malloc(L * sizeof(double complex*));
  for(int l = 0; l < L; l++)
  {
    U[l] = (double complex*)malloc(d * d * sizeof(double complex));
  }
  U = generate_random_unitary_matrices(L);

  // Testing output of generate_random_unitary_matrices
  printf("Printing L random unitary matrices\n");
    for (int i = 0; i < L; i++) {
      printf("i=%d\n", i);
      for (int j = 0; j < d*d; j++) {
        printf("(%lf, %lf)\n", creal(U[i][j]), cimag(U[i][j]));
      }
    }

  // Evaluate PGM, U(PGM), and bme
  for (size_t x = 0; x < L; x++)
  {
    E_pgm[x] = (double complex*) calloc(d*d, sizeof(double complex));
    E_pgm[x] = multiply_three_matrices(d, sum_pinv_sqrt, p_rho[x], sum_pinv_sqrt);

    U_E_pgm[x] = (double complex*) calloc(d*d, sizeof(double complex));
    U_E_pgm[x] = apply_unital_channel(L, d, U, E_pgm[x]);

    bme[x] = (double complex*) calloc(d*d, sizeof(double complex));
    bme_U[x] = (double complex*) calloc(d*d, sizeof(double complex));
    bme[x] = compute_Bayes_mean_sum(x, E_pgm[x], rho, prior_constant);
    bme_U[x] = compute_Bayes_mean_sum(x, U_E_pgm[x], rho, prior_constant);
  }

  // Testing output of multiply_three_matrices
  printf("Printing PGM POVM for\n");
    for (int x = 0; x < L; x++) {
      printf("x=%d\n", x);
      for (int j = 0; j < d*d; j++) {
        printf("(%lf, %lf)\n", creal(E_pgm[x][j]), cimag(E_pgm[x][j]));
      }
    }

  // Testing output of compute_Bayes_mean_sum
  printf("Printing Bayes estimator and rho for\n");
    for (int x = 0; x < L; x++) {
      printf("x=%d\n", x);
      for (int j = 0; j < d*d; j++) {
        printf("BME: (%lf, %lf), rho:(%lf, %lf)\n", creal(bme[x][j]), cimag(bme[x][j]), creal(rho[x][j]), cimag(rho[x][j]));
      }
    }

  double* fidelities_bme = (double*) calloc(L, sizeof(double));
  double* fidelities_bme_U = (double*) calloc(L, sizeof(double));
  double* fidelities = (double*) calloc(L, sizeof(double));

  fidelities_bme = compute_fidelities(rho[3], bme);
  fidelities_bme_U = compute_fidelities(rho[3], bme_U);
  fidelities = compute_fidelities(rho[3], rho);

    // Open a file for writing
    char filename[256];
    sprintf(filename, "fidelities_BME_Unital_PGM_d=%d_L=%d.txt",d,L);
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("Error opening file\n");
        return 1;
    }

    // Write the fidelities, fidelities_bme to the file
    for(int i = 0; i < L; i++) {
        fprintf(file, "%lf %lf %lf\n", fidelities[i], fidelities_bme[i], fidelities_bme_U[i]);
    }

    // double fid = fidelity(d, rho[20], rho[20]);
    // printf("Manual test for fidelity(rho[20], rho[20]):%f\n", fid);

    // Close the file
    fclose(file);

  // Free dynamic variables

  for(int l = 0; l < L; l++)
  {
    free(rho[l]);
    free(U[l]);
    free(E_pgm[l]);
    free(U_E_pgm[l]);
    free(bme[l]);
  }

  free(rho);
  free(U);
  free(E_pgm);
  free(U_E_pgm);
  free(bme);
  free(bme_U);

  free(fidelities_bme);
  free(fidelities);
  free(prior_constant);

  return 0;
}

