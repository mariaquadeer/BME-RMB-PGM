//=============================================================
// Copyright 2023 Maria Quadeer
// Distributed under the GNU General Public License, Version 3.
// (See accompanying file LICENSE.txt or copy at
// https://www.gnu.org/licenses/gpl-3.0.en.html)
//=============================================================

#include "global.h"

/* Original declaration for extern variables. */
int d; /* dimension */
int L; /* no. of states in the ensemble */
int N0; /* sample complexity/no. of measurements */
const double zero_epsilon = 1e-6; /* Defining 0 for comparison with float types. */;

// Define Clifford group elements as static constant arrays
lapack_complex_double Id[4] = {{1, 0}, {0, 0}, {0, 0}, {1, 0}};
lapack_complex_double X[4] = {{0, 0}, {1, 0}, {1, 0}, {0, 0}};
lapack_complex_double Y[4] = {{0, 0}, {0, -1}, {0, 1}, {0, 0}};
lapack_complex_double Z[4] = {{1, 0}, {0, 0}, {0, 0}, {-1, 0}};
lapack_complex_double negI[4] = {{-1, 0}, {0, 0}, {0, 0}, {-1, 0}};
lapack_complex_double negX[4] = {{0, 0}, {-1, 0}, {-1, 0}, {0, 0}};
lapack_complex_double negY[4] = {{0, 0}, {0, 1}, {0, -1}, {0, 0}};
lapack_complex_double negZ[4] = {{-1, 0}, {0, 0}, {0, 0}, {1, 0}};
lapack_complex_double iI[4] = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
lapack_complex_double iX[4] = {{0, 0}, {0, 1}, {0, 1}, {0, 0}};
lapack_complex_double iY[4] = {{0, 0}, {0, 1}, {-1, 0}, {0, 0}};
lapack_complex_double iZ[4] = {{0, 1}, {0, 0}, {0, 0}, {0, -1}};
lapack_complex_double neg_iI[4] = {{0, -1}, {0, 0}, {0, 0}, {0, -1}};
lapack_complex_double neg_iX[4] = {{0, 0}, {0, -1}, {0, -1}, {0, 0}};
lapack_complex_double neg_iY[4] = {{0, 0}, {0, -1}, {0, 1}, {0, 0}};
lapack_complex_double neg_iZ[4] = {{0, -1}, {0, 0}, {0, 0}, {0, 1}};
lapack_complex_double Hadamard_X[4] = {{0.707106781, 0.707106781}, {0, 0}, {0, 0}, {0.707106781, -0.707106781}};
lapack_complex_double Hadamard_X_neg[4] = {{0.707106781, -0.707106781}, {0, 0}, {0, 0}, {0.707106781, 0.707106781}};
lapack_complex_double Hadamard_Y[4] = {{0.707106781, 0.707106781}, {0, 0}, {0, 0}, {-0.707106781, 0.707106781}};
lapack_complex_double Hadamard_Y_neg[4] = {{0.707106781, -0.707106781}, {0, 0}, {0, 0}, {0.707106781, 0.707106781}};
lapack_complex_double Hadamard_Z[4] = {{0.707106781, 0.707106781}, {0, 0}, {0, 0}, {0.707106781, 0.707106781}};
lapack_complex_double Hadamard_Z_neg[4] = {{0.707106781, -0.707106781}, {0, 0}, {0, 0}, {0.707106781, 0.707106781}};
lapack_complex_double XY_combo[4] = {{0.707106781, 0.707106781}, {0, 0}, {0, 0}, {0.707106781, 0.707106781}};
lapack_complex_double XY_combo_neg[4] = {{0.707106781, -0.707106781}, {0, 0}, {0, 0}, {0.707106781, 0.707106781}};

int main(int argc, char* argv[])
  {

    if (argc < 4) {
          fprintf(stderr, "Usage: %s <d> <L> <N>\n", argv[0]);
          exit(EXIT_FAILURE);
      }

      d = atoi(argv[1]); /* dimension */
      L = atoi(argv[2]); /* no. of states in the ensemble */
      N0 = atoi(argv[3]); /* no. of measurements */


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

      //for (int i = 0; i < d*d; i++) {
        // printf("iX[%d] = (%lf, %lf)\n", i, creal(iX[i]), cimag(iX[i]));
        // fflush(stdout);
      //}

      // Seed random number generator
      srand(time(NULL));

      // Initialise average_risk and risk
      double average_risk_clifford = 0.0;
      double average_risk_rmb = 0.0;
      double risk_clifford = 0.0;
      double risk_rmb = 0.0;

     // Allocate memory for pure states, prior and posterior, etc.
     lapack_complex_double* mixed_states[L];
     double prior_clifford[L];
     double prior_rmb[L];
     double posterior_clifford[L][d];
     double posterior_rmb[L][d];
     double conditional_distribution_clifford[L][d];
     double conditional_distribution_rmb[L][d];
     double p_x_clifford[d];
     double p_x_rmb[d];

     // Initialize pure states and uniform prior
     for (int a = 0; a < L; a++) {
       mixed_states[a] = func_Ginibre_random_rho();
       //for (int i = 0; i < d*d; i++) {
         // printf("mixed_states[%d][%d] = (%lf, %lf)\n", a, i, creal(mixed_states[a][i]), cimag(mixed_states[a][i]));
         // fflush(stdout);
       //}
       prior_clifford[a] = 1.0 / L;
       prior_rmb[a] = 1.0 / L;
     }

     // Declare basis B
     lapack_complex_double* B_clifford;
     lapack_complex_double* B_rmb;

     // Perform Bayesian updates
     for (int n = 0; n < N0; n++) {

       // Free B from the previous iteration
       if (n > 0) {
         free(B_clifford);
         free(B_rmb);
         }

      B_clifford = random_clifford_qubit();
      B_rmb = func_Haar_random_B();

      // Evaluate conditional distribution and p(x)
      for (int x = 0; x < d; x++) {
        p_x_clifford[x] = 0;
        p_x_rmb[x] = 0;

        for (int a = 0; a < L; a++) {
          conditional_distribution_clifford[a][x] = fidelity_pure_state_density_matrix(d, &B_clifford[x * d], mixed_states[a]);
          conditional_distribution_rmb[a][x] = fidelity_pure_state_density_matrix(d, &B_rmb[x * d], mixed_states[a]);
          p_x_clifford[x] += conditional_distribution_clifford[a][x] * prior_clifford[a];
          p_x_rmb[x] += conditional_distribution_rmb[a][x] * prior_rmb[a];
          // printf("p_x[%d] = %lf\n", x, p_x[x]);
        }
        // Evaluate posterior for the N-basis measurement
        if(n == N0-1){
          for (int a = 0; a < L; a++){
            posterior_clifford[a][x] = conditional_distribution_clifford[a][x] * prior_clifford[a] / p_x_clifford[x];
            posterior_rmb[a][x] = conditional_distribution_rmb[a][x] * prior_rmb[a] / p_x_rmb[x];
            // printf("posterior[%d][%d] = %lf\n", a, x, posterior[a][x]);
          }
        }
      }

      // Perform inverse sampling
      int x_n_clifford = func_inverse_sampling_quantum_states(p_x_clifford);
      int x_n_rmb = func_inverse_sampling_quantum_states(p_x_rmb);

       // Update prior only N0-1 times.
       if(n < N0-1){
         double normalization_factor_clifford = p_x_clifford[x_n_clifford];
         double normalization_factor_rmb = p_x_rmb[x_n_rmb];
         for (int a = 0; a < L; a++) {
           prior_clifford[a] = conditional_distribution_clifford[a][x_n_clifford] * prior_clifford[a] / normalization_factor_clifford;
           prior_rmb[a] = conditional_distribution_rmb[a][x_n_rmb] * prior_rmb[a] / normalization_factor_rmb;
           }
       }
     } // end of N0 updates


    // Compute Bayes estimator for the last measurement n = N0 - 1. Evaluate risk and average_risk
    lapack_complex_double* Bayes_est_clifford = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* Bayes_est_rmb = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    for (int x = 0; x < d; x++) {
    // Reset Bayes_est for each x
    memset(Bayes_est_clifford, 0, d * d * sizeof(lapack_complex_double));
    memset(Bayes_est_rmb, 0, d * d * sizeof(lapack_complex_double));

    for (int a = 0; a < L; a++) {
      // Multiply pure_states[a] with its conjugate transpose and scale by prior[a] and add to Bayes_est; GEMM takes care of sum.

      lapack_complex_double temp_clifford = lapack_make_complex_double(posterior_clifford[a][x], 0.0);
      lapack_complex_double temp_rmb = lapack_make_complex_double(posterior_rmb[a][x], 0.0);

      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp_clifford, mixed_states[a], d, id, d, &one, Bayes_est_clifford, d);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp_rmb, mixed_states[a], d, id, d, &one, Bayes_est_rmb, d);
    }

    // Evaluate risk of Bayes_est for pure_states[0] as a running sum over x
      double fid_0_clifford = fidelity_mixed_simple(d, mixed_states[0], Bayes_est_clifford);
      double fid_0_rmb = fidelity_mixed_simple(d, mixed_states[0], Bayes_est_rmb);
      risk_clifford += conditional_distribution_clifford[0][x] * (fid_0_clifford);
      risk_rmb += conditional_distribution_rmb[0][x] * (fid_0_rmb);

    // Evaluate average_risk as running sum
      for (int a = 0; a < L; a++) {
      double fid_clifford = fidelity_mixed_simple(d, mixed_states[a], Bayes_est_clifford);
      double fid_rmb = fidelity_mixed_simple(d, mixed_states[a], Bayes_est_rmb);
      // printf("fidelity_mixed(%d, %d) = %lf\n", a, x, fid);
      average_risk_clifford += prior_clifford[a] * conditional_distribution_clifford[a][x] * (fid_clifford);
      average_risk_rmb += prior_rmb[a] * conditional_distribution_rmb[a][x] * (fid_rmb);
     }
    }

    // Open a file for writing
    char filename_clifford[256];
    char filename_rmb[256];
    snprintf(filename_clifford, sizeof(filename_clifford), "Compare_Clifford_Ginibre_average_risk_d=%d_N=%d.txt", d, N0);
    snprintf(filename_rmb, sizeof(filename_rmb), "Compare_RMB_Ginibre_average_risk_d=%d_N=%d.txt", d, N0);
    FILE* file_clifford = fopen(filename_clifford, "w");
    FILE* file_rmb = fopen(filename_rmb, "w");
    fprintf(file_clifford, "Risk: %lf Average Risk: %lf\n", risk_clifford, average_risk_clifford);
    fprintf(file_rmb, "Risk: %lf Average Risk: %lf\n", risk_rmb, average_risk_rmb);
    fclose(file_clifford);
    fclose(file_rmb);

    // Free allocated memory for Bayes estimator and pure states
    free(Bayes_est_clifford);
    free(Bayes_est_rmb);
    for (int a = 0; a < L; a++) {
      free(mixed_states[a]);
    }

    // Free allocated memory for B and id
    free(B_clifford);
    free(B_rmb);
    free(id);

  return 0;
}

