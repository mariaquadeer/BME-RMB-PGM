#include "global.h"

/* Original declaration for extern variables. */
int d; /* dimension */
int L; /* no. of states in the ensemble */
int N0; /* sample complexity/no. of measurements */
const double zero_epsilon = 1e-6; /* Defining 0 for comparison with float types. */;

// Pauli group modulo Identity
/*=============================================================================*/
// Pauli X
lapack_complex_double X[4] = {{0, 0}, {1, 0}, {1, 0}, {0, 0}};

// Pauli Y
lapack_complex_double Y[4] = {{0, 0}, {0, -1}, {0, 1}, {0, 0}};

// Pauli Z
lapack_complex_double Z[4] = {{1, 0}, {0, 0}, {0, 0}, {-1, 0}};

/*=============================================================================*/

// UNITARY 2-DESIGN IN d=2 // Ref. My first paper
/*=============================================================================*/

    // U_0,0 (Identity)
    lapack_complex_double U_0_0[4] = {{1, 0}, {0, 0}, {0, 0}, {1, 0}};

    // U_π,0 (Pauli X)
    lapack_complex_double U_pi_0[4] = {{0, 0}, {1, 0}, {1, 0}, {0, 0}};

    // U_π/2,0 (Hadamard)
    lapack_complex_double U_pi_2_0[4] = {{0.707106781, 0}, {0.707106781, 0}, {-0.707106781, 0}, {0.707106781, 0}};

    // U_π/2,π (Hadamard followed by Z, or S²)
    lapack_complex_double U_pi_2_pi[4] = {{0.707106781, 0}, {-0.707106781, 0}, {0.707106781, 0}, {0.707106781, 0}};

    // U_π/2,π/2 (Hadamard followed by S)
    lapack_complex_double U_pi_2_pi_2[4] = {{0.707106781, 0}, {0, 0.707106781}, {0, 0.707106781}, {0.707106781, 0}};

    // U_π/2,3π/2 (Hadamard followed by S†)
    lapack_complex_double U_pi_2_3pi_2[4] = {{0.707106781, 0}, {0, -0.707106781}, {0, -0.707106781}, {0.707106781, 0}};

/*=============================================================================*/
// Define Clifford group elements in d = 2 as static constant arrays
/*=============================================================================*/
// Identity
lapack_complex_double Id[4] = {{1, 0}, {0, 0}, {0, 0}, {1, 0}};

// H
lapack_complex_double H[4] = {{0.707107, 0}, {0.707107, 0}, {0.707107, 0}, {-0.707107, 0}};

// S
lapack_complex_double S[4] = {{1, 0}, {0, 0}, {0, 0}, {0, 1}};

// HS
lapack_complex_double HS[4] = {{0.707107, 0}, {0, 0.707107}, {0.707107, 0}, {0, -0.707107}};

// SH
lapack_complex_double SH[4] = {{0.707107, 0}, {0, 0.707107}, {0.707107, 0}, {0, -0.707107}};

// SS(Z)
lapack_complex_double SS[4] = {{1, 0}, {0, 0}, {0, 0}, {-1, 0}};

// HSH
lapack_complex_double HSH[4] = {{0.5, 0.5}, {0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}};

// HSS (HZ)
lapack_complex_double HSS[4] = {{0.707107, 0}, {-0.707107, 0}, {0.707107, 0}, {0.707107, 0}};

// SHS
lapack_complex_double SHS[4] = {{0.707107, 0}, {0, 0.707107}, {0, 0.707107}, {0.707107, 0}};

// SSH
lapack_complex_double SSH[4] = {{0.707107, 0}, {0.707107, 0}, {-0.707107, 0}, {0.707107, 0}};

// SSS (Sd)
lapack_complex_double SSS[4] = {{1, 0}, {0, 0}, {0, 0}, {0, -1}};

// HSHS
lapack_complex_double HSHS[4] = {{0.5, 0.5}, {0.5, 0.5}, {0.5, -0.5}, {-0.5, 0.5}};

// HSSH
lapack_complex_double HSSH[4] = {{0, 0}, {1, 0}, {1, 0}, {0, 0}};

// HSSS
lapack_complex_double HSSS[4] = {{0.707107, 0}, {0, -0.707107}, {0.707107, 0}, {0, 0.707107}};

// SHSS
lapack_complex_double SHSS[4] = {{0.707107, 0}, {-0.707107, 0}, {0, 0.707107}, {0, 0.707107}};

// SSHS
lapack_complex_double SSHS[4] = {{0.707107, 0}, {0, 0.707107}, {-0.707107, 0}, {0, 0.707107}};

// HSHSS
lapack_complex_double HSHSS[4] = {{0.5, 0.5}, {-0.5, 0.5}, {0.5, -0.5}, {-0.5, -0.5}};

// HSSHS
lapack_complex_double HSSHS[4] = {{0, 0}, {0, 1}, {1, 0}, {0, 0}};

// SHSSH
lapack_complex_double SHSSH[4] = {{0, 0}, {1, 0}, {0, 1}, {0, 0}};

// SHSSS
lapack_complex_double SHSSS[4] = {{0.707107, 0}, {0, -0.707107}, {0, 0.707107}, {-0.707107, 0}};

// SSHSS
lapack_complex_double SSHSS[4] = {{0.707107, 0}, {-0.707107, 0}, {-0.707107, 0}, {-0.707107, 0}};

// HSHSSH
lapack_complex_double HSHSSH[4] = {{0, 0.707107}, {0.707107, 0}, {0, -0.707107}, {0.707107, 0}};

// HSHSSS
lapack_complex_double HSHSSS[4] = {{0.5, 0.5}, {-0.5, -0.5}, {0.5, -0.5}, {0.5, -0.5}};

// HSSHSS
lapack_complex_double HSSHSS[4] = {{0, 0}, {-1, 0}, {1, 0}, {0, 0}};
/*=============================================================================*/

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
      double average_risk_2Design = 0.0;
      double average_risk_1Design = 0.0;
      double average_risk_rmb = 0.0;
      double risk_clifford = 0.0;
      double risk_2Design = 0.0;
      double risk_1Design = 0.0;
      double risk_rmb = 0.0;

     // Allocate memory for pure states, prior and posterior, etc.
     lapack_complex_double* mixed_states[L];
     double prior_clifford[L];
     double prior_2Design[L];
     double prior_1Design[L];
     double prior_rmb[L];
     double posterior_clifford[L][d];
     double posterior_2Design[L][d];
     double posterior_1Design[L][d];
     double posterior_rmb[L][d];
     double conditional_distribution_clifford[L][d];
     double conditional_distribution_2Design[L][d];
     double conditional_distribution_1Design[L][d];
     double conditional_distribution_rmb[L][d];
     double p_x_clifford[d];
     double p_x_2Design[d];
     double p_x_1Design[d];
     double p_x_rmb[d];

     // Initialize pure states and uniform prior
     for (int a = 0; a < L; a++) {
       mixed_states[a] = func_Ginibre_random_rho();
       //for (int i = 0; i < d*d; i++) {
         // printf("mixed_states[%d][%d] = (%lf, %lf)\n", a, i, creal(mixed_states[a][i]), cimag(mixed_states[a][i]));
         // fflush(stdout);
       //}
       prior_clifford[a] = 1.0 / L;
       prior_2Design[a] = 1.0 / L;
       prior_1Design[a] = 1.0 / L;
       prior_rmb[a] = 1.0 / L;
     }

     // Declare basis B
     lapack_complex_double* B_clifford;
     lapack_complex_double* B_2Design;
     lapack_complex_double* B_1Design;
     lapack_complex_double* B_rmb;

     // Perform Bayesian updates
     for (int n = 0; n < N0; n++) {

       // Free B from the previous iteration
       if (n > 0) {
         free(B_clifford);
         free(B_2Design);
         free(B_1Design);
         free(B_rmb);
         }

      B_clifford = random_clifford_qubit();
      B_2Design = random_2Design_qubit();
      B_1Design = random_1Design_qubit();
      B_rmb = func_Haar_random_B();

      // Evaluate conditional distribution and p(x)
      for (int x = 0; x < d; x++) {
        p_x_clifford[x] = 0;
        p_x_2Design[x] = 0;
        p_x_1Design[x] = 0;
        p_x_rmb[x] = 0;

        for (int a = 0; a < L; a++) {
          conditional_distribution_clifford[a][x] = fidelity_pure_state_density_matrix(d, &B_clifford[x * d], mixed_states[a]);
          conditional_distribution_2Design[a][x] = fidelity_pure_state_density_matrix(d, &B_2Design[x * d], mixed_states[a]);
          conditional_distribution_1Design[a][x] = fidelity_pure_state_density_matrix(d, &B_1Design[x * d], mixed_states[a]);
          conditional_distribution_rmb[a][x] = fidelity_pure_state_density_matrix(d, &B_rmb[x * d], mixed_states[a]);
          p_x_clifford[x] += conditional_distribution_clifford[a][x] * prior_clifford[a];
          p_x_2Design[x] += conditional_distribution_2Design[a][x] * prior_2Design[a];
          p_x_1Design[x] += conditional_distribution_1Design[a][x] * prior_1Design[a];
          p_x_rmb[x] += conditional_distribution_rmb[a][x] * prior_rmb[a];
          // printf("p_x[%d] = %lf\n", x, p_x[x]);
        }
        // Evaluate posterior for the N-basis measurement
        if(n == N0-1){
          for (int a = 0; a < L; a++){
            posterior_clifford[a][x] = conditional_distribution_clifford[a][x] * prior_clifford[a] / p_x_clifford[x];
            posterior_2Design[a][x] = conditional_distribution_2Design[a][x] * prior_2Design[a] / p_x_2Design[x];
            posterior_1Design[a][x] = conditional_distribution_1Design[a][x] * prior_1Design[a] / p_x_1Design[x];
            posterior_rmb[a][x] = conditional_distribution_rmb[a][x] * prior_rmb[a] / p_x_rmb[x];
            // printf("posterior[%d][%d] = %lf\n", a, x, posterior[a][x]);
          }
        }
      }

      // Perform inverse sampling
      int x_n_clifford = func_inverse_sampling_quantum_states(p_x_clifford);
      int x_n_2Design = func_inverse_sampling_quantum_states(p_x_2Design);
      int x_n_1Design = func_inverse_sampling_quantum_states(p_x_1Design);
      int x_n_rmb = func_inverse_sampling_quantum_states(p_x_rmb);

       // Update prior only N0-1 times.
       if(n < N0-1){
         double normalization_factor_clifford = p_x_clifford[x_n_clifford];
         double normalization_factor_2Design = p_x_2Design[x_n_2Design];
         double normalization_factor_1Design = p_x_1Design[x_n_1Design];
         double normalization_factor_rmb = p_x_rmb[x_n_rmb];
         for (int a = 0; a < L; a++) {
           prior_clifford[a] = conditional_distribution_clifford[a][x_n_clifford] * prior_clifford[a] / normalization_factor_clifford;
           prior_2Design[a] = conditional_distribution_2Design[a][x_n_2Design] * prior_2Design[a] / normalization_factor_2Design;
           prior_1Design[a] = conditional_distribution_1Design[a][x_n_1Design] * prior_1Design[a] / normalization_factor_1Design;
           prior_rmb[a] = conditional_distribution_rmb[a][x_n_rmb] * prior_rmb[a] / normalization_factor_rmb;
           }
       }
     } // end of N0 updates


    // Compute Bayes estimator for the last measurement n = N0 - 1. Evaluate risk and average_risk
    lapack_complex_double* Bayes_est_clifford = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* Bayes_est_2Design = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* Bayes_est_1Design = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* Bayes_est_rmb = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    for (int x = 0; x < d; x++) {
    // Reset Bayes_est for each x
    memset(Bayes_est_clifford, 0, d * d * sizeof(lapack_complex_double));
    memset(Bayes_est_2Design, 0, d * d * sizeof(lapack_complex_double));
    memset(Bayes_est_1Design, 0, d * d * sizeof(lapack_complex_double));
    memset(Bayes_est_rmb, 0, d * d * sizeof(lapack_complex_double));

    for (int a = 0; a < L; a++) {
      // Multiply pure_states[a] with its conjugate transpose and scale by prior[a] and add to Bayes_est; GEMM takes care of sum.

      lapack_complex_double temp_clifford = lapack_make_complex_double(posterior_clifford[a][x], 0.0);
      lapack_complex_double temp_2Design = lapack_make_complex_double(posterior_2Design[a][x], 0.0);
      lapack_complex_double temp_1Design = lapack_make_complex_double(posterior_1Design[a][x], 0.0);
      lapack_complex_double temp_rmb = lapack_make_complex_double(posterior_rmb[a][x], 0.0);

      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp_clifford, mixed_states[a], d, id, d, &one, Bayes_est_clifford, d);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp_2Design, mixed_states[a], d, id, d, &one, Bayes_est_2Design, d);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp_1Design, mixed_states[a], d, id, d, &one, Bayes_est_1Design, d);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &temp_rmb, mixed_states[a], d, id, d, &one, Bayes_est_rmb, d);
    }

    // Evaluate risk of Bayes_est for pure_states[0] as a running sum over x
      double fid_0_clifford = fidelity_mixed_simple(d, mixed_states[0], Bayes_est_clifford);
      double fid_0_2Design = fidelity_mixed_simple(d, mixed_states[0], Bayes_est_2Design);
      double fid_0_1Design = fidelity_mixed_simple(d, mixed_states[0], Bayes_est_1Design);
      double fid_0_rmb = fidelity_mixed_simple(d, mixed_states[0], Bayes_est_rmb);
      risk_clifford += conditional_distribution_clifford[0][x] * (fid_0_clifford);
      risk_2Design += conditional_distribution_2Design[0][x] * (fid_0_2Design);
      risk_1Design += conditional_distribution_1Design[0][x] * (fid_0_1Design);
      risk_rmb += conditional_distribution_rmb[0][x] * (fid_0_rmb);

    // Evaluate average_risk as running sum
      for (int a = 0; a < L; a++) {
      double fid_clifford = fidelity_mixed_simple(d, mixed_states[a], Bayes_est_clifford);
      double fid_2Design = fidelity_mixed_simple(d, mixed_states[a], Bayes_est_2Design);
      double fid_1Design = fidelity_mixed_simple(d, mixed_states[a], Bayes_est_1Design);
      double fid_rmb = fidelity_mixed_simple(d, mixed_states[a], Bayes_est_rmb);
      // printf("fidelity_mixed(%d, %d) = %lf\n", a, x, fid);
      average_risk_clifford += prior_clifford[a] * conditional_distribution_clifford[a][x] * (fid_clifford);
      average_risk_2Design += prior_2Design[a] * conditional_distribution_2Design[a][x] * (fid_2Design);
      average_risk_1Design += prior_1Design[a] * conditional_distribution_1Design[a][x] * (fid_1Design);
      average_risk_rmb += prior_rmb[a] * conditional_distribution_rmb[a][x] * (fid_rmb);
     }
    }

    // Open a file for writing
    char filename_clifford[256];
    char filename_2Design[256];
    char filename_1Design[256];
    char filename_rmb[256];
    snprintf(filename_clifford, sizeof(filename_clifford), "Compare_Clifford_Ginibre_average_risk_d=%d_N=%d.txt", d, N0);
    snprintf(filename_2Design, sizeof(filename_2Design), "Compare_2Design_Ginibre_average_risk_d=%d_N=%d.txt", d, N0);
    snprintf(filename_1Design, sizeof(filename_1Design), "Compare_1Design_Ginibre_average_risk_d=%d_N=%d.txt", d, N0);
    snprintf(filename_rmb, sizeof(filename_rmb), "Compare_RMB_Ginibre_average_risk_d=%d_N=%d.txt", d, N0);
    FILE* file_clifford = fopen(filename_clifford, "w");
    FILE* file_2Design = fopen(filename_2Design, "w");
    FILE* file_1Design = fopen(filename_1Design, "w");
    FILE* file_rmb = fopen(filename_rmb, "w");
    fprintf(file_clifford, "Risk: %lf Average Risk: %lf\n", risk_clifford, average_risk_clifford);
    fprintf(file_2Design, "Risk: %lf Average Risk: %lf\n", risk_2Design, average_risk_2Design);
    fprintf(file_1Design, "Risk: %lf Average Risk: %lf\n", risk_1Design, average_risk_1Design);
    fprintf(file_rmb, "Risk: %lf Average Risk: %lf\n", risk_rmb, average_risk_rmb);
    fclose(file_clifford);
    fclose(file_2Design);
    fclose(file_1Design);
    fclose(file_rmb);

    // Free allocated memory for Bayes estimator and pure states
    free(Bayes_est_clifford);
    free(Bayes_est_2Design);
    free(Bayes_est_1Design);
    free(Bayes_est_rmb);
    for (int a = 0; a < L; a++) {
      free(mixed_states[a]);
    }

    // Free allocated memory for B and id
    free(B_clifford);
    free(B_2Design);
    free(B_1Design);
    free(B_rmb);
    free(id);

  return 0;
}
