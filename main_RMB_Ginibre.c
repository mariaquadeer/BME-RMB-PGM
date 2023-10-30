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
int N0; /* sample complexity/no. of measurements */
int I0; /* The data_I0 to put the output files in*/
const double zero_epsilon = 1e-6; /* Defining 0 for comparison with float types. */;

int main(int argc, char* argv[])
{

  struct rlimit rl;

  // Get the current stack size limit
  if (getrlimit(RLIMIT_STACK, &rl) != 0) {
    perror("getrlimit");
    return 1;
  }

  // Set the current (soft) limit to the maximum (hard) limit
  rl.rlim_cur = rl.rlim_max;

  // Apply the new limit
  if (setrlimit(RLIMIT_STACK, &rl) != 0) {
    perror("setrlimit");
    fprintf(stderr, "Failed to set the stack size to the maximum allowed value\n");
    return 1;
  }

  if (argc < 5) {
    fprintf(stderr, "Usage: %s <d> <L> <N> <I0>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  d = atoi(argv[1]); /* dimension */
  L = atoi(argv[2]); /* no. of states in the ensemble */
  N0 = atoi(argv[3]); /* no. of measurements */
  I0 = atoi(argv[4]); /* output files in data_I0 folder */

  //printf("STARTING MAIN():\n");
  // getchar();

  /****** Define global variables as declared outside main().******/

  // Define one, zero, i1 and identity matrix id.
  one = lapack_make_complex_double(1.0,0.0);
  zero = lapack_make_complex_double(0.0,0.0);
  i1 = lapack_make_complex_double(0.0,1.0);

  // Calloc initalizes array elements as zero while malloc does not. Malloc is faster than calloc.
  id = (lapack_complex_double *) calloc(d*d, sizeof(lapack_complex_double));
  if (!id) {
    perror("calloc failed for id");
    exit(EXIT_FAILURE);
  }

  // printf("printing identity matrix\n");
  for(int j=0;j<d;j++) // defining id.
  {
    for(int i=0;i<d;i++)
    {
      if (i==j) id[i+j*d] = one;
        else id[i+j*d]=zero;
      // printf("%lf, %lf\n", creal(id[i+j*d]), cimag(id[i+j*d]));
    }
  }


  // Initialise average_risk and risk
  double average_risk = 0.0;
  double risk = 0.0;

  // Allocate memory for pure states, prior and posterior, etc.
  lapack_complex_double* mixed_states[L];
  double prior[L];
  double posterior[L][d];
  double conditional_distribution[L][d];
  double p_x[d];

  // Initialize pure states and uniform prior
  for (int a = 0; a < L; a++) {
    mixed_states[a] = func_Ginibre_random_rho();
    // for (int i = 0; i < d*d; i++) {
    //   printf("mixed_states[%d][%d] = (%lf, %lf)\n", a, i, creal(mixed_states[a][i]), cimag(mixed_states[a][i]));
    //   fflush(stdout);
    // }
    prior[a] = 1.0 / L;
  }

  // Declare basis B
  lapack_complex_double* B;

  // Perform Bayesian updates
  for (int n = 0; n < N0; n++) {

    // Free B from the previous iteration
    if (n > 0) {
      free(B);
    }

    B = func_Haar_random_B();
    // for (int i = 0; i < d*d; i++) {
    //   printf("B[%d] = (%lf, %lf)\n", i, creal(B[i]), cimag(B[i]));
    //   fflush(stdout);
    // }

    // Evaluate conditional distribution and p(x)
    for (int x = 0; x < d; x++) {
      p_x[x] = 0;

      for (int a = 0; a < L; a++) {
        conditional_distribution[a][x] = fidelity_pure_state_density_matrix(d, &B[x * d], mixed_states[a]);
        // printf("conditional_distribution[%d][%d] = %lf\n", a, x, conditional_distribution[a][x]);
        p_x[x] += conditional_distribution[a][x] * prior[a];
        // printf("p_x[%d] = %lf\n", x, p_x[x]);
      }

      // Evaluate posterior for the N-basis measurement
      if(n == N0-1){
        for (int a = 0; a < L; a++){
          posterior[a][x] = conditional_distribution[a][x] * prior[a] / p_x[x];
          // printf("posterior[%d][%d] = %lf\n", a, x, posterior[a][x]);
        }
      }
    }

    // Perform inverse sampling
    int x_n = func_inverse_sampling_quantum_states(p_x);
    // printf("x_n = %d\n", x_n);

    // Update prior only N0-1 times.
    if(n < N0-1){
      double normalization_factor = p_x[x_n];
      for (int a = 0; a < L; a++) {
        prior[a] = conditional_distribution[a][x_n] * prior[a] / normalization_factor;
        // printf("updated prior = %lf\n", prior[a]);
      }
    }

    // Check for normalization
    //  if (is_normalized(prior, L)) {
    //   printf("The probability distribution is normalized.\n");
    // } else {
    //   printf("The probability distribution is not normalized.\n");
    // }

  } // end of N0 updates


  // Compute Bayes estimator for the last measurement n = N0 - 1. Evaluate risk and average_risk
  lapack_complex_double* Bayes_est = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
  for (int x = 0; x < d; x++) {
    // Reset Bayes_est for each x
    memset(Bayes_est, 0, d * d * sizeof(lapack_complex_double));

    for (int a = 0; a < L; a++) {
      // Multiply pure_states[a] with its conjugate transpose and scale by prior[a] and add to Bayes_est; GEMM takes care of sum.
      lapack_complex_double temp = lapack_make_complex_double(posterior[a][x], 0.0);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp, mixed_states[a], d, id, d, &one, Bayes_est, d);
    }

    // test print Bayes_est
    // printf("x = %d\n", x);
    // for (int i = 0; i < d*d; i++) {
    //   printf("Bayes_est[%d] = (%lf, %lf)\n", i, creal(Bayes_est[i]), cimag(Bayes_est[i]));
    //   fflush(stdout);
    // }

    // Evaluate risk of Bayes_est for pure_states[0] as a running sum over x
    double fid_0 = fidelity_mixed_simple(d, mixed_states[0], Bayes_est);
    risk += conditional_distribution[0][x] * (fid_0);

    // Evaluate average_risk as running sum
    for (int a = 0; a < L; a++) {
      double fid = fidelity_mixed_simple(d, mixed_states[a], Bayes_est);
      // printf("fidelity_mixed(%d, %d) = %lf\n", a, x, fid);
      average_risk += prior[a] * conditional_distribution[a][x] * (fid);
    }
  }

  // Open a file for writing
  char filename[256];
  snprintf(filename, sizeof(filename), "plots_RMB/Ginibre/data_%d/Ginibre_average_risk_d=%d_N=%d_L=%d.txt", I0, d, N0, L);

  // Create the directory structure if it doesn't exist
  char temp_path[256];
  char *ptr;
  strncpy(temp_path, filename, sizeof(temp_path));
  ptr = strchr(temp_path, '/');
  while (ptr != NULL) {
    *ptr = '\0';
    // Attempt to create the directory
    mkdir(temp_path, 0755);
    *ptr = '/';
    ptr = strchr(ptr + 1, '/');
  }

  FILE* file = fopen(filename, "w");
  fprintf(file, "Risk: %lf Average Risk: %lf\n", risk, average_risk);
  fclose(file);

  // Free allocated memory for Bayes estimator and pure states
  free(Bayes_est);
  for (int a = 0; a < L; a++) {
    free(mixed_states[a]);
  }

  // Free allocated memory for B and id
  free(B);
  free(id);

  return 0;
}

