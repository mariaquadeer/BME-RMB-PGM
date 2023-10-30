//=============================================================
// Copyright 2023 Maria Quadeer
// Distributed under the GNU General Public License, Version 3.
// (See accompanying file LICENSE.txt or copy at
// https://www.gnu.org/licenses/gpl-3.0.en.html)
//=============================================================

/* custom function library */
#include "global.h"

void strev(int* str)
{
    int len = 2;

    for (int i = 0; i < len/2; i++)
    {
        int temp = str[i];
        str[i] = str[len-i-1];
        str[len-i-1] = temp;
    }
}


double min(double A[],int size)
{
	double A_min=A[0];
	for(int i=1;i<size;i++)
	{
	if (A[i]<A_min) A_min=A[i];
	}
return A_min;
}


double max(double A[],int size)
{
	double A_max=A[0];
	for(int i=1;i<size;i++)
	{
	if(A[i]>A_max) A_max=A[i];
	}
return A_max;
}

/* The following function converts a sorted probability mass function (pmf) to a cumulative mass function and re-writes on pmf. */
void convert_pmf_to_cmf(double *pmf, int count)
{
  for(int i = 1; i < count; i++)
  {
    pmf[i]+= pmf[i-1];
  }
}


// function to check if a distribution is normalized
int is_normalized(double *prob_dist, int size)
{
  double sum = 0.0;
  for (int i = 0; i < size; i++) {
    sum += prob_dist[i];
  }

  // Comparing with a small tolerance to account for floating-point inaccuracies
  if (fabs(sum - 1.0) < 1e-6) {
    return 1; // Normalized
  } else {
    return 0; // Not normalized
  }
}


void print_matrix(lapack_complex_double* matrix, int size)
{
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("(%f, %f) ", creal(matrix[i*size + j]), cimag(matrix[i*size + j]));
        }
        printf("\n");
    }
}


lapack_complex_double* random_clifford_qubit()
{

    // List of all Clifford elements for a single qubit
    lapack_complex_double* clifford_elements[24] = {Id, X, Y, Z, negI, negX, negY, negZ, iI, iX, iY, iZ, neg_iI, neg_iX, neg_iY, neg_iZ, Hadamard_X, Hadamard_X_neg, Hadamard_Y, Hadamard_Y_neg, Hadamard_Z, Hadamard_Z_neg, XY_combo, XY_combo_neg};

    // Randomly select a Clifford element
    int choice = rand() % 24;
    // printf("choice = %d\n", choice);
    lapack_complex_double* selected_matrix = clifford_elements[choice];

    // Allocate space for result and copy the selected matrix
    lapack_complex_double* result = malloc(4 * sizeof(lapack_complex_double));
    for (int i = 0; i < 4; ++i) {
        result[i] = selected_matrix[i];
    }

    // for (int i = 0; i < d*d; i++) {
    //   printf("clifford_elements[14][%d] = (%lf, %lf)\n", i, creal(clifford_elements[10][i]), cimag(clifford_elements[10][i]));
    //   fflush(stdout);
    // }

    return result;
}


int func_inverse_sampling_quantum_states(double *prior)
{
	//printf("******Entered func_inverse_sampling_quantum_states.******\n");

  int x = 0;
	double u;
	double* pmf;
	pmf = (double *) calloc(d, sizeof(double));
  if (!pmf) {
    perror("calloc failed for pmf");
    exit(EXIT_FAILURE);
  }
	for (int k = 0; k < d; k++)
	{
		pmf[k] = prior[k];
	//	printf("printing pmf[%d]: %lf\n", k, pmf[k]);
	}

  LAPACKE_dlasrt('I', d, pmf);
	for(int i=0; i<d; i++) //printf("printing sorted pmf[%d]: %lf\n", i, pmf[i]);

  convert_pmf_to_cmf(pmf, d);
	for(int i=0; i<d; i++) //printf("pmf to cmf convert step completed, printing pmf_updated[%d]:%lf\n", i, pmf[i]);

  // Sample 'u' unif. from [0,1].
  u = genrand_real1();

  // Compute cmf_inverse at 'u' : cmf_inverse(u) = inf {i ; u <= pmf[i]}.
  for (int i= 0; i < d; i++)
  {
    if(u <= pmf[i])
    {
     x = i;
     break;
    }
  }

  free(pmf);
  pmf = NULL;

//  printf("******Leaving func_inverse_sampling_quantum_states().******\n");

  return x;
}


static inline bool isEqual_toZero(double x)
{
	//printf("******Entered isEqual().******\n");
	// getchar();

	double epsilon = 1e-2; /* some small number such as 1e-5 for comparing float types. */;
	// printf("epsilon: %lf\n", epsilon);
	// getchar();

	bool flag;
	flag = fabs(x - zero_epsilon) <= epsilon*fabs(x); // see Knuth section 4.2.2 pages 217-218

	//printf("flag: %d\n", flag);
	// getchar();

	return flag;
}


double func_rel_ent(lapack_complex_double* rho, lapack_complex_double* Bayes_est)
{
	//printf("******Entered func_rel_ent().******\n");
	// getchar();

	lapack_int lda = d;
	lapack_int info1, info2;		// lda: leading dimension of matrix
	double rel_ent = 0.0;
	double rel_ent_temp = 0.0;
	double *eigenval_rho_copy, *eigenval_Bayes_est_copy;
	eigenval_rho_copy = (double*) calloc(d, sizeof(double));
  if (!eigenval_rho_copy) {
    perror("calloc failed for eigenval_rho_copy");
    exit(EXIT_FAILURE);
  }
	eigenval_Bayes_est_copy = (double*) calloc(d, sizeof(double));
  if (!eigenval_Bayes_est_copy) {
    perror("calloc failed for eigenval_Bayes_est_copy");
    exit(EXIT_FAILURE);
  }
  lapack_complex_double *rho_copy, *Bayes_est_copy;

	rho_copy = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));
  if (!rho_copy) {
    perror("calloc failed for rho_copy");
    exit(EXIT_FAILURE);
  }
	Bayes_est_copy = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));
  if (!Bayes_est_copy) {
    perror("calloc failed for Bayes_est_copy");
    exit(EXIT_FAILURE);
  }

	cblas_zcopy(d*d, rho, 1, rho_copy, 1);
	cblas_zcopy(d*d, Bayes_est, 1, Bayes_est_copy, 1);


	info1 = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', d, rho_copy, lda, eigenval_rho_copy);

	info2 = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', d, Bayes_est_copy, lda, eigenval_Bayes_est_copy);

// **********TESTING***************** //
	//{
	// if( info2 > 0 ) {
	// 							printf( "The algorithm failed to compute eigenvalues.\n" );
	// 							exit( 1 );
	// 			}
	// else{
	//
	// 	printf("The algorithm computed eigenvalues.\n" );
  //}

		printf("Smallest eigenvalue of Bayes est: %lf\n",eigenval_Bayes_est_copy[0]);
		printf("Smallest eigenvalue of rho: %lf\n",eigenval_rho_copy[0]);
		// getchar();
// **********TESTING COMPLETE ***************** //

	 if(!isEqual_toZero(eigenval_Bayes_est_copy[0]))
		{
			for(int j=0; j<d; j++)
		  	{
			  	for(int i=0; i<d; i++)
			 			{
							for(int k=0;k<d;k++)
							{
								rel_ent_temp += conj(rho_copy[k+j*d])*Bayes_est_copy[k+i*d];
			 				}
							rel_ent += - (eigenval_rho_copy[j]*log(eigenval_Bayes_est_copy[i])*pow(fabs(rel_ent_temp), 2));
							rel_ent_temp = 0.0;
		 				}
	 			}
		}
	 else // when the smallest eigenvalues of both rho and Bayes_est_copy are zero
	 {
	 	rel_ent = INFINITY;
	 }
	printf("rel_ent: %lf.\n", rel_ent);
	return rel_ent;

}


lapack_complex_double* func_Haar_random_B()
{
  lapack_int lda = d; lapack_int info;		// lda: leading dimension of matrix
	lapack_complex_double *U, *R, *tau, *B; /* B is the Haar distributed unitary matrix that defines our basis. */
	U = (lapack_complex_double *) calloc(d*d, sizeof(lapack_complex_double)); // calloc returns void* pointer type; casting isn't really needed.
  if (!U) {
    perror("calloc failed for U");
    exit(EXIT_FAILURE);
  }
	R = (lapack_complex_double *) calloc(d*d, sizeof(lapack_complex_double));
  if (!R) {
    perror("calloc failed for R");
    exit(EXIT_FAILURE);
  }
	B = (lapack_complex_double *) calloc(d*d, sizeof(lapack_complex_double));
  if (!B) {
    perror("calloc failed for B");
    exit(EXIT_FAILURE);
  }
	tau = (lapack_complex_double*) calloc(d, sizeof(lapack_complex_double));
  if (!tau) {
    perror("calloc failed for tau");
    exit(EXIT_FAILURE);
  }

	double x1,x2;
  double y1,y2;
  double w;

  //unsigned long seed = (unsigned long)time(NULL) ^ (unsigned long)getpid() ^ (unsigned long)clock();

  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  unsigned long seed = (unsigned long)ts.tv_nsec ^ (unsigned long)getpid() ^ (unsigned long)clock();
  init_genrand(seed);

	/* Box- Mueller to generate normally distributed random numbers for defining the random matrix elements of U. */
  for(int j=0;j<d;j++)
	{
  	for(int i=0;i<d;i++)
		{
    	do{
      		x1=genrand_real1();
      		x2=genrand_real1();

					x1=2.0*x1-1.0;
      		x2=2.0*x2-1.0;

					w=x1*x1+x2*x2;

				}while(w>=1.0);

			w=sqrt((-2.0*log(w))/w);

			y1=x1*w;
	    y2=x2*w;
			U[i+j*d]=lapack_make_complex_double(y1,y2);  //printf("%lf %lf\n",creal(U1[i+j*n]),cimag(U1[i+j*n]));
		}
	}
// WHY ARE WE DOING QR DECOMPOSITION? The matrix U is a matrix with random entries. Q part of QR decomposition is Unitary!


	info = LAPACKE_zgeqrf(LAPACK_COL_MAJOR, d, d, U, lda, tau); // complex*16 General matrix QR Factorization

	for(int j=0;j<d;j++)
	{
	 	for(int i=0;i<d;i++)
		{
			if (i==j) R[i+j*d]=U[i+j*d]/(cabs(U[i+j*d])); else R[i+j*d] = zero;	//printf("%lf %lf\n",creal(r[i+j*n]),cimag(r[i+j*n]));
		}
	}

	info = LAPACKE_zungqr(LAPACK_COL_MAJOR, d, d, d, U, lda, tau);
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, U, lda, R, lda, &zero, B, d); // this step ensures we have a Haar-distributed U written on B.

  free(U);
  free(R);
  free(tau);

  return B;

}

struct risk_types func_Bayes_update_blind(int n)
{
  printf("**********************************************\n");
  fflush(stdout);
	printf("******Entered func_Bayes_update_blind().******\n");
  fflush(stdout);
  printf("**********************************************\n");
  fflush(stdout);

	struct risk_types r;
	lapack_complex_double **trace, **Bayes_est;
  double **log_likelihood;

	r.risk = 0.0, r.average_risk = 0.0;

	trace = (lapack_complex_double **) calloc(d, sizeof(lapack_complex_double*));
  if (!trace) {
    perror("calloc failed for trace");
    exit(EXIT_FAILURE);
  }

  log_likelihood = (double **) calloc(d, sizeof(double*));
  if (!log_likelihood) {
    perror("calloc failed for likelihood");
    exit(EXIT_FAILURE);
  }

	Bayes_est = (lapack_complex_double **) calloc(d, sizeof(lapack_complex_double*));
  if (!Bayes_est) {
    perror("calloc failed for Bayes_est");
    exit(EXIT_FAILURE);
  }

	// Loop over basis vector indexed by j, {|b_j>} and copy each vector into basis_vec_temp.
		for(int j=0; j < d; j++)
		{
      // Evaluate Bayes_est for j-th basis vector.
		  // printf("******Evaluating Bayes estimator for j=%d.******\n",j);
      // printf("\n");
      // getchar();

      lapack_complex_double* basis_vec_temp = (lapack_complex_double *) calloc(d, sizeof(lapack_complex_double));
      if (!basis_vec_temp) {
        perror("calloc failed for basis_vec_temp");
        exit(EXIT_FAILURE);
      }

      // loop over elements of j-th vector and copy corresponding d-entries of B into basis_vec_temp.
			for(int k=0; k<d; k++)
			{
			     basis_vec_temp[k] = basis[k+j*d]; // measurement basis
			     // printf("printing (re, imag) of measurement basis_vec_temp[%d]: (%lf, %lf).\n", k, creal(basis_vec_temp[k]), cimag(basis_vec_temp[k]));
			     // getchar();
			}

      // evaluate trace[], log_likelihood[], and posterior[]
      trace[j] = (lapack_complex_double*) calloc(L*d, sizeof(lapack_complex_double));
      log_likelihood[j] = (double*) calloc(L*d, sizeof(double));
      if (!log_likelihood[j] || !trace[j]) {
        perror("calloc failed for log_likelihood[j] or trace[j] or posterior[j]");
        exit(EXIT_FAILURE);
      }

      // start loop over states
      for (int l = 0; l < L; l++)
      {
        for(int i = 0; i < d; i++)
  		  {
  			  // for computing tr(P_j \rho_i) = trace[i+j*d+]: ZDOTC <b_j|rho[l][i]> and then mod square it.
  			  trace[j][i+l*d] = cblas_zdotc(d, basis_vec_temp, 1, rho[l][i], 1);
  		  	trace[j][i+l*d] = pow(cabs(trace[j][i+l*d]), 2)/((double) d);

  	  		// printf("(re, imag) of trace[%d + d*%d + %d*d]: (%lf, %lf).\n", i, j, l, creal(trace[i+j*d+l*d]), cimag(trace[i+j*d+l*d]));
  		  	// getchar();

  		  	// evaluate likelihood
          if (n == 0){
            log_likelihood[j][i+l*d] = log(trace[j][i+l*d]) + log(prior_constant);
          }
          else {
            log_likelihood[j][i+l*d] = log(trace[j][i+l*d]) + log_prior[i+l*d];
          }

  	  		// printf("(re, imag) of likelihood[%d + %d*d + %d*d]: (%lf, %lf).\n", i, j, l, creal(exp(log_likelihood[i+j*d+l*d])), cimag(exp(log_likelihood[i+j*d+l*d])));
  		  	// getchar();

        } // end of i-loop
      } // end of l-loop

      // Normalize log-likelihood to ensure it sums up to 1 in the log-space
      normalize_log_probs(log_likelihood[j], L*d);

     // When n = N0-1 obtain Bayes_est[j]
     // Bayes_est = \sum_{l,i} posterior[i+j*d+l*d]*rho[l][i]. ZGEMM takes care of the sum.
     Bayes_est[j] = (lapack_complex_double *) calloc(d*d, sizeof(lapack_complex_double));
     if (!Bayes_est[j]) {
       perror("calloc failed for Bayes_est[j]");
       exit(EXIT_FAILURE);
     }
     if (n == N0-1) // evaluate the Bayes estimator
     {
       lapack_complex_double* density = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));
       if (!density) {
         perror("calloc failed for density");
         exit(EXIT_FAILURE);
       }
       for (int l = 0; l < L; l++)
       {
         for(int i=0; i < d; i++)
         {
           lapack_complex_double posterior;
           posterior = lapack_make_complex_double(exp(log_likelihood[j][i+l*d]), 0.0);
           // printf("printing posterior[%d][%d + %d*d]: (%lf, %lf).\n", j, i, l, creal(posterior), cimag(posterior));
           // getchar();

           // Using gemm for vec --> op.
           cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, 1, &one, rho[l][i], d, rho[l][i], d, &zero, density, d);

           // Multiply density(l, i) with posterior and add to Bayes_est[j].
           cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &posterior, density, d, id, d, &one, Bayes_est[j], d);

           // evaluate average risk for Bayes_est[j]
           r.average_risk += creal(posterior)*creal(trace[j][i+l*d])*fidelity_pure_state_density_matrix(d, rho[l][i], Bayes_est[j]);

         } // end of i-loop
       } // end of l-loop

       // free density
       free(density);
       density = NULL;

      //  printing Bayes_est[j]:
       printf("Printing elements of Bayes_est[%d] for testing.\n", j);
       fflush(stdout);
       getchar();
       for(int k=0; k<d; k++)
       {
         for(int i=0; i<d; i++)
         {
           printf("(re, imag) of Bayes_est[%d][%d][%d]: (%lf, %lf).\n",j, i, k, creal(Bayes_est[j][i+k*d]), cimag(Bayes_est[j][i+k*d]));
           fflush(stdout);
           getchar();
         }
       }

      // evaluate risk for an arbitrary state from the ensemble: compute trace between rho[1][0] and outcome_prev and average over this for all outcomes j
      //r.risk += trace[1+0*d+j*d]*fidelity_pure_state_density_matrix(d, rho[1][0], Bayes_est[j]);

     } // end of if condition

    // free basis_vec_temp
    free(basis_vec_temp);
    basis_vec_temp = NULL;
    } // end of j-loop

  // print r.average_risk
  // printf("Printing average risk: %lf\n", r.average_risk);
  // getchar();

    // measurement simulation by inverse sampling from the total probability of outcome[j]
    double sum_measure_prob = 0.0;
    double *measure_prob;
    measure_prob = (double*) calloc(d, sizeof(double));
    for (int j = 0; j < d; j++) {
      for(int l = 0; l < L; l++) {
        for(int i = 0; i < d; i++){
          measure_prob[j] += exp(log_likelihood[j][i+l*d]);
        }
      }
      sum_measure_prob += measure_prob[j];
    }

    for (int j = 0; j < d; j++) {
      measure_prob[j] /= sum_measure_prob;
    //  printf("measure_prob[%d] = %lf\n", j, measure_prob[j]);
    }

    int outcome_prev = func_inverse_sampling_quantum_states(measure_prob);

    // free measure_prob
    free(measure_prob);
    measure_prob = NULL;

    // update log_prior
    cblas_zcopy(L*d, log_likelihood[outcome_prev], 1, log_prior, 1);

    // printf("\n");
    // printf("Measurement outcome: %d\n", outcome_prev);
    // printf("\n");

    printf("******Leaving func_Bayes_update_blind().******\n");
    fflush(stdout);
	  printf("\n");
    fflush(stdout);

  // Freeing the inner allocated arrays
  for(int j = 0; j < d; j++){
    free(trace[j]);
    trace[j] = NULL;
    free(log_likelihood[j]);
    log_likelihood[j] = NULL;
    free(Bayes_est[j]);
    Bayes_est[j] = NULL;
  }

  // Freeing the outer allocated arrays
   free(trace);
   trace = NULL;
   free(log_likelihood);
   log_likelihood = NULL;
   free(Bayes_est);
   Bayes_est = NULL;

  return r;

}

lapack_complex_double* func_Ginibre_random_rho()
{
    lapack_complex_double *B; /* B is the Ginibre distributed matrix. */
    lapack_complex_double *rho; /* rho will be the resulting density matrix. */

    B = (lapack_complex_double *) calloc(d*d,sizeof(lapack_complex_double));
    rho = (lapack_complex_double *) calloc(d*d,sizeof(lapack_complex_double));

    if (!B || !rho) {
        perror("calloc failed for B or rho");
        exit(EXIT_FAILURE);
    }

    double x1,x2;
    double y1,y2;
    double w;

    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    unsigned long seed = (unsigned long)ts.tv_nsec ^ (unsigned long)getpid() ^ (unsigned long)clock();
    init_genrand(seed);
    // printf("seed = %lu\n", seed);

    /* Box-Muller to generate normally distributed random numbers for defining the random matrix elements of B. */
    for(int j=0;j<d;j++)
    {
        for(int i=0;i<d;i++)
        {
            do{
                x1=genrand_real1();
                x2=genrand_real1();

                x1=2.0*x1-1.0;
                x2=2.0*x2-1.0;

                w=x1*x1+x2*x2;

            }while(w>=1.0);

            w=sqrt((-2.0*log(w))/w);

            y1=x1*w;
            y2=x2*w;
            B[i+j*d]=lapack_make_complex_double(y1,y2);
        }
    }

    // for (int j = 0; j < d*d; j++) {
    //   printf("B[%d]=(%lf, %lf)\n", j, creal(B[j]), cimag(B[j]));
    // }

    /* Compute the Hermitian matrix B * B^H to get a positive semi-definite matrix. */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &one, B, d, B, d, &zero, rho, d);

    free(B);
    B = NULL;

    /* Compute the trace of rho to normalize it to a trace-1 density matrix. */
    lapack_complex_double trace = 0.0;
    for(int i=0;i<d;i++)
    {
        trace += rho[i+i*d];
    }

    double trace_inv = 1.0 / creal(trace); // assume that the trace is a real number.

    /* Normalize rho to make it a density matrix. */
    for(int j=0;j<d;j++)
    {
        for(int i=0;i<d;i++)
        {
            rho[i+j*d] = rho[i+j*d] * trace_inv;
        }
    }

    // for (int j = 0; j < d*d; j++) {
    //   printf("rho[%d]=(%lf, %lf)\n", j, creal(rho[j]), cimag(rho[j]));
    // }

    return rho;
}

lapack_complex_double* func_Haar_random_state()
{
  lapack_int lda = d; lapack_int info;		// lda: leading dimension of matrix
	lapack_complex_double *U, *R, *tau, *B; /* B is the Haar distributed unitary matrix that defines our basis. */
	U = (lapack_complex_double *) calloc(d*d,sizeof(lapack_complex_double)); // calloc returns void* pointer type; casting isn't really needed.
  if (!U) {
    perror("calloc failed for U");
    exit(EXIT_FAILURE);
  }
	R = (lapack_complex_double *) calloc(d*d,sizeof(lapack_complex_double));
  if (!R) {
    perror("calloc failed for R");
    exit(EXIT_FAILURE);
  }
	B = (lapack_complex_double *) calloc(d*d,sizeof(lapack_complex_double));
  if (!B) {
    perror("calloc failed for B");
    exit(EXIT_FAILURE);
  }
	tau = (lapack_complex_double*) calloc(d, sizeof(lapack_complex_double));
  if (!tau) {
    perror("calloc failed for tau");
    exit(EXIT_FAILURE);
  }

	double x1,x2;
  double y1,y2;
  double w;

  //unsigned long seed = (unsigned long)time(NULL) ^ (unsigned long)getpid() ^ (unsigned long)clock();

  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  unsigned long seed = (unsigned long)ts.tv_nsec ^ (unsigned long)getpid() ^ (unsigned long)clock();
  init_genrand(seed);
  //printf("seed = %lu\n", seed);


	/* Box- Mueller to generate normally distributed random numbers for defining the random matrix elements of U. */
  for(int j=0;j<d;j++)
	{
  	for(int i=0;i<d;i++)
		{
    	do{
      		x1=genrand_real1();
      		x2=genrand_real1();

					x1=2.0*x1-1.0;
      		x2=2.0*x2-1.0;

					w=x1*x1+x2*x2;

				}while(w>=1.0);

			w=sqrt((-2.0*log(w))/w);

			y1=x1*w;
	    y2=x2*w;
			U[i+j*d]=lapack_make_complex_double(y1,y2);  //printf("%lf %lf\n",creal(U1[i+j*n]),cimag(U1[i+j*n]));
		}
	}
// WHY ARE WE DOING QR DECOMPOSITION? The matrix U is a matrix with random entries. Q part of QR decomposition is Unitary!


	info = LAPACKE_zgeqrf(LAPACK_COL_MAJOR,d,d,U,lda,tau); // complex*16 General matrix QR Factorization

	for(int j=0;j<d;j++)
	{
	 	for(int i=0;i<d;i++)
		{
			if (i==j) R[i+j*d]=U[i+j*d]/(cabs(U[i+j*d]));else R[i+j*d]=zero;	//printf("%lf %lf\n",creal(r[i+j*n]),cimag(r[i+j*n]));
		}
	}

	info = LAPACKE_zungqr(LAPACK_COL_MAJOR,d,d,d,U,lda,tau);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,d,d,d,&one,U,lda,R,lda,&zero,B,d); // this step ensures we have a Haar-distributed U written on B.

  // Extract first column as random pure state
  lapack_complex_double* state = (lapack_complex_double*) malloc(d * sizeof(lapack_complex_double));
  for (int i = 0; i < d; i++) {
    state[i] = B[i];
  }

 // Free allocated memory
  free(U);
  U = NULL;
  free(R);
  R = NULL;
  free(tau);
  tau = NULL;
  free(B);
  B = NULL;

  return state;

}

struct risk_types func_Bayes_update_blind_Ginibre(int n)
{
  printf("******************************************************\n");
  fflush(stdout);
	printf("******Entered func_Bayes_update_blind_Ginibre().******\n");
  fflush(stdout);
  printf("******************************************************\n");
  fflush(stdout);

	struct risk_types r;
	lapack_complex_double **trace, **Bayes_est;
  double **log_likelihood;

	r.risk = 0.0, r.average_risk = 0.0;

	trace = (lapack_complex_double **) calloc(d, sizeof(lapack_complex_double*));
  if (!trace) {
    perror("calloc failed for trace");
    exit(EXIT_FAILURE);
  }

	Bayes_est = (lapack_complex_double **) calloc(d, sizeof(lapack_complex_double*));
  if (!Bayes_est) {
    perror("calloc failed for Bayes_est");
    exit(EXIT_FAILURE);
  }

  log_likelihood = (double **) calloc(d, sizeof(double*));
  if (!log_likelihood) {
    perror("calloc failed for likelihood");
    exit(EXIT_FAILURE);
  }

	// Loop over basis vector indexed by j, {B_j} and copy each vector into basis_vec_temp.
		for(int j=0; j < d; j++)
		{
      //Evaluate Bayes_est for j-th basis vector.
		  // printf("******Evaluating log_likelihood for j=%d.******\n",j);
      // printf("\n");
      //getchar();

      // loop over elements of j-th vector and copy corresponding d-entries of measurement basis into basis_vec_temp
      lapack_complex_double* basis_vec_temp = (lapack_complex_double *) calloc(d, sizeof(lapack_complex_double));
      if (!basis_vec_temp) {
        perror("calloc failed for basis_vec_temp");
        exit(EXIT_FAILURE);
      }

			for(int k = 0; k < d; k++)
			{
			  basis_vec_temp[k] = basis[k+j*d];
        // printf("basis_vec_temp[%d] = (%lf, %lf).\n", k, creal(basis_vec_temp[k]), cimag(basis_vec_temp[k]));
        //getchar();
		  }

      // evaluate trace[], log_likelihood[], and posterior[]
      trace[j] = (lapack_complex_double*) calloc(L*d, sizeof(lapack_complex_double));

      log_likelihood[j] = (double*) calloc(L*d, sizeof(double));

      if (!log_likelihood[j] || !trace[j]) {
        perror("calloc failed for log_likelihood[j] or trace[j]");
        exit(EXIT_FAILURE);
      }

      // loop over states
      for (int l = 0; l < L; l++)
      {
        for(int i = 0; i < d; i++) // for computing tr(basis_vec_temp \rho_{l,i}) = trace[j][i+j*d]:
  		  {
          // Temporarily storing the result of rho_mixed[l][i] * basis_vec_temp.
          lapack_complex_double* temp_vec = (lapack_complex_double*) calloc(d, sizeof(lapack_complex_double));
          if (!temp_vec) {
            perror("calloc failed for temp_vec");
            exit(EXIT_FAILURE);
          }
          // Note that rho_mixed[l][i] should be a square matrix and basis_vec_temp a vector.
          cblas_zgemv(CblasColMajor, CblasNoTrans, d, d, &one, rho_mixed[l][i], d, basis_vec_temp, 1, &zero, temp_vec, 1);

          // for (int i = 0; i < d; i++) {
          //   printf("temp_vec[%d] = (%lf, %lf).\n", i, creal(temp_vec[i]), cimag(temp_vec[i]));
          //   getchar();
          // }

          // Compute the dot product of basis_vec_temp and temp_vec
          trace[j][i+l*d] = cblas_zdotc(d, basis_vec_temp, 1, temp_vec, 1);

          // free temp_vec
          free(temp_vec);
          temp_vec = NULL;

          // free basis_vec_temp
          free(basis_vec_temp);
          basis_vec_temp = NULL;

          trace[j][i+l*d] = pow(cabs(trace[j][i+l*d]), 2)/((double) d);

          // printf("Printing trace[%d+ %d*d + %d*d] = (%lf, %lf).\n", i, j, l, creal(trace[j][i+l*d]), cimag(trace[j][i+l*d]));
          // getchar();

  		    // evaluate log_likelihood

          if (n == 0){
              // printf("n==0!\n");
              log_likelihood[j][i+l*d] = log(trace[j][i+l*d]) + log(prior_constant);
            }
          else{
              // printf("n!=0!\n");
              log_likelihood[j][i+l*d] = log(trace[j][i+l*d]) + log_prior_Ginibre[i+l*d];
            }

          // printf("log_likelihood[%d][%d+%d*d]: %lf\n", j, i, l, log_likelihood[j][i+l*d]);

        } // end of i-loop
      } // end of l-loop

      // Normalize log-likelihood to ensure it sums up to 1 in the log-space
      normalize_log_probs(log_likelihood[j], L*d);

      // When n = N0-1 obtain Bayes_est[j]
		  // Bayes_est = \sum_{l,i} posterior[j][i+l*d]*rho[l][i]. ZGEMM takes care of the sum.
      Bayes_est[j] = (lapack_complex_double *) calloc(d*d, sizeof(lapack_complex_double));
      if (!Bayes_est[j]) {
        perror("calloc failed for Bayes_est[j]");
        exit(EXIT_FAILURE);
      }
      if (n == N0-1)
      {
       // printf("last update\n");
       // printf("\n");
       for (int l = 0; l < L; l++)
       {
         for(int i = 0; i < d; i++)
         {
           lapack_complex_double posterior;
           posterior = lapack_make_complex_double(exp(log_likelihood[j][i+l*d]), 0.0);
           // printf("printing posterior[%d][%d + %d*d]: (%lf, %lf).\n", j, i, l, creal(posterior), cimag(posterior));
           // getchar();

           // Multiply rho_mixed[l][i] with posterior and add to Bayes_est[j].
           cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &posterior, rho_mixed[l][i], d, id, d, &one, Bayes_est[j], d);

           // evaluate average risk for Bayes_est[j]
           r.average_risk += creal(posterior)*creal(trace[j][i+l*d])*fidelity_mixed(d, rho_mixed[l][i], Bayes_est[j]);
         }
       }

      //printing Bayes_est[j]:
      printf("Printing elements of Bayes_est[%d] for testing.\n", j);
      fflush(stdout);
      getchar();

      for(int k=0; k<d; k++){
         for(int i=0; i<d; i++){
            printf("(re, imag) of Bayes_est[%d][%d][%d]: (%lf, %lf).\n",j,i,k, creal(Bayes_est[j][i+k*d]), cimag(Bayes_est[j][i+k*d]));
            fflush(stdout);
            getchar();
         }
       }

       // print average_risk
       // printf("Average risk = %lf\n", r.average_risk);

       // evaluate risk for an arbitrary state: compute trace between rho_mixed[1][1] and outcome_prev and average over this for all outcomes j
       //r.risk += trace[1+0*d+j*d]*fidelity_mixed(d, rho_mixed[1][0], Bayes_est[j]);

      } // end of if condition


    } // end of j-loop


     // measurement simulation by inverse sampling from the total probability of outcome[j]
     double sum_measure_prob = 0.0;
     double *measure_prob;
     measure_prob = (double*) calloc(d, sizeof(double));
     for (int j = 0; j < d; j++) {
       for(int l = 0; l < L; l++) {
         for(int i = 0; i < d; i++){
           measure_prob[j] += exp(log_likelihood[j][i+l*d]);
         }
       }
      sum_measure_prob += measure_prob[j];
     }

     for (int j = 0; j < d; j++) {
      measure_prob[j] /= sum_measure_prob;
     }

     int outcome_prev = func_inverse_sampling_quantum_states(measure_prob);

     // free measure_prob
     free(measure_prob);
     measure_prob = NULL;

     // update log_prior_Ginibre
     cblas_zcopy(L*d, log_likelihood[outcome_prev], 1, log_prior_Ginibre, 1);

     // printf("\n");
     // printf("Measurement outcome: %d\n", outcome_prev);
     // printf("\n");

    printf("******Leaving func_Bayes_update_blind_Ginibre().******\n");
    fflush(stdout);
	  printf("\n");
    fflush(stdout);

    for (int j = 0; j < d; j++) {
      free(trace[j]);
      trace[j] = NULL;
      free(log_likelihood[j]);
      log_likelihood[j] = NULL;
      free(Bayes_est[j]);
      Bayes_est[j] = NULL;
    }

    free(trace);
    trace = NULL;
    free(log_likelihood);
    log_likelihood = NULL;
    free(Bayes_est);
    Bayes_est = NULL;

    return r;

}


// Function to check if a matrix is PSD and has trace 1
int check_psd_trace_one(lapack_complex_double* mat, int n)
{
    printf("Entered check_psd_trace_one()....\n");
    double sum_eig = 0.0;

    // Allocate space for eigenvalues
    double* eig = (double*) malloc(n * sizeof(double));
    if (!eig) {
        perror("malloc failed for eig");
        exit(EXIT_FAILURE);
    }

    // Compute eigenvalues with LAPACKE zheevd routine
    lapack_int info = LAPACKE_zheevd(LAPACK_COL_MAJOR, 'N', 'U', n, mat, n, eig);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zheevd failed with info = %d\n", info);
        exit(EXIT_FAILURE);
    }

    // Check if all eigenvalues are non-negative (PSD check) and compute their sum (trace check)

    for (int i = 0; i < n; i++) {
        printf("eig[%d] = %lf.\n", i, eig[i]);
        if (eig[i] < 0.0) {
            free(eig);
            eig = NULL;
            return 0;  // Not PSD
        }
        sum_eig += eig[i];
    }

    free(eig);
    eig = NULL;
    printf("PSD! Now checking trace...\n");
    // Check if trace is 1 (within EPSILON tolerance)
    if (fabs(sum_eig - 1.0) < EPSILON) {
        return 1;  // PSD and trace 1
    } else {
        return 0;  // Not trace 1
    }
}


// Function to calculate the logarithm of the sum of exponentials (log-sum-exp)
double log_sum_exp(double *log_probs, int n)
{
    double max_val = log_probs[0];
    for (int i = 1; i < n; i++) {
        if (log_probs[i] > max_val) {
            max_val = log_probs[i];
        }
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += exp(log_probs[i] - max_val);
    }

    return max_val + log(sum);
}

// Function to normalize log-probabilities to ensure they sum up to 1
void normalize_log_probs(double *log_probs, int n)
{
    double log_sum = log_sum_exp(log_probs, n);
    for (int i = 0; i < n; i++) {
        log_probs[i] -= log_sum;
    }
}


/* Function to convert a given number between 0 and d^n
 to a base d number: 'dernary'. */
int* from_dAlphabet_to_dBase(int *res, int base, int exponent, int inputNum)
{
		printf("entering from_dAlphabet_to_dBase\n");
	//	getchar();
    int index = exponent - 1;  // Initialize index of result array res[] as last location of res[] to ensure leading zeroes and a top-to-bottom readout after the base conversion.

    // repeatedly dividing inputNum by base and taking remainder
		if(inputNum == 0)
		{
			printf("inputNum is zero!\n");
		//	getchar();
			printf("printing res[0]: %d\n", res[0]);
		//	getchar();
			return res;
		}
		else{
			while (inputNum > 0)
    	{
        res[index--] = inputNum % base;
        inputNum /= base;
    	}
		}

		return res;
}


double fidelity_pure_state_density_matrix(int n, lapack_complex_double* psi, lapack_complex_double* rho_original)
{
    // Copy rho_original to another array so we can get eigenvectors without altering rho_original
    lapack_complex_double* rho = (lapack_complex_double*) malloc(n * n * sizeof(lapack_complex_double));
    memcpy(rho, rho_original, n * n * sizeof(lapack_complex_double));

    // Step 1: Compute Eigenvalues and Eigenvectors
    double* w = (double*) malloc(n * sizeof(double));
    lapack_int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', n, rho, n, w);
    // printf("info = %d\n", info);

    // Step 2: Threshold small eigenvalues
    double epsilon = 1e-12;  // Threshold for setting eigenvalues to zero
    for (int i = 0; i < n; i++) {
        // printf("w[%d]=%lf\n", i, w[i]);
        if (w[i] < epsilon) {
            w[i] = 0.0;
        }
    }

    // Step 2: Compute Adjusted Fidelity
    lapack_complex_double* psi_prime = (lapack_complex_double*) malloc(n * sizeof(lapack_complex_double));
    cblas_zgemv(CblasColMajor, CblasConjTrans, n, n, &one, rho, n, psi, 1, &zero, psi_prime, 1);

    double fid = 0.0;
    for (int i = 0; i < n; i++) {
        fid += w[i] * (creal(psi_prime[i]) * creal(psi_prime[i]) + cimag(psi_prime[i]) * cimag(psi_prime[i]));
    }

    free(w);
    free(rho);
    free(psi_prime);

    return fid;

}

// Function to compute square root of a matrix
lapack_complex_double* sqrtm_LAPACK(int d, lapack_complex_double* A)
{
    // printf("Entered sqrtm_LAPACK...\n");
    // for (int j = 0; j < d*d; j++) {
    //     printf("A[%d]=(%lf, %lf)\n", j, creal(A[j]), cimag(A[j]));
    // }

    // int psd = check_psd_trace_one(A, d); // this modifies mat so comment after testing!
    // printf("psd test = %d\n", psd);

    lapack_complex_double *sqrtA = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    lapack_complex_double *eigvec = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    double *w = (double*) calloc(d, sizeof(double));

    lapack_complex_double *eigvec_rescale = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    cblas_zcopy(d*d, A, 1, eigvec, 1);

    for (int j = 0; j < d*d; j++) {
        // printf("eigvec[%d]=(%lf, %lf)\n", j, creal(eigvec[j]), cimag(eigvec[j]));
    }

    // Perform the eigenvalue decomposition
    lapack_int info = LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', d, eigvec, d, w);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zheevd failed with info = %d\n", info);
        exit(EXIT_FAILURE);
    }

    cblas_zcopy(d*d, eigvec, 1, eigvec_rescale, 1);

    // Take square root of eigenvalues and scale the eigenvectors
    for(int i = 0; i < d; i++) {
        // printf("w[%d] = %lf\n", i, w[i]);
        w[i] = sqrt(w[i]);
        cblas_zdscal(d, w[i], &eigvec_rescale[i*d], 1);
    }

    // Reconstruct sqrtA = eigvec_rescale * eigvec^H
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &one, eigvec_rescale, d, eigvec, d, &zero, sqrtA, d);

    free(w);
    w = NULL;
    free(eigvec);
    eigvec = NULL;
    free(eigvec_rescale);
    eigvec_rescale = NULL;

    return sqrtA;
}

// Function to compute square root of a matrix using Schur decomposition
lapack_complex_double* sqrtm_LAPACK_Schur(int d, lapack_complex_double* A)
{
    // Allocate memory for Q, T, and sqrtA
    lapack_complex_double* Q = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* T = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* sqrtA = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* w = (lapack_complex_double*) calloc(d, sizeof(lapack_complex_double));
    lapack_complex_double* temp = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));

    // Check if allocation succeeded
    if (!Q || !T || !sqrtA || !temp || !w) {
        perror("calloc failed");
        exit(EXIT_FAILURE);
    }

    cblas_zcopy(d*d, A, 1, T, 1);

    // Compute the Schur decomposition
    // LAPACKE_zgees (matrix_layout, jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs);
    lapack_int sdim;
    // Now, pass &sdim as the parameter
    LAPACKE_zgees(LAPACK_COL_MAJOR, 'V', 'N', NULL, d, T, d, &sdim, w, Q, d);

    // Compute the square root of T by taking square root of each element on the diagonal
    for (int i = 0; i < d; i++) {
        T[i*d+i] = csqrt(T[i*d+i]);
        printf("W[%d] = (%lf, %lf)\n", i, creal(w[i]), cimag(w[i]));
        printf("sqrt:T[%d] = (%lf, %lf)\n", i, creal(T[i*d+i]), cimag(T[i*d+i]));
    }

    // Compute sqrt(A) = Q * sqrt(T) * Q*
    // temp = Q * sqrt(T)
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, Q, d, T, d, &zero, temp, d);
    // sqrtA = temp * Q*
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &one, temp, d, Q, d, &zero, sqrtA, d);

    // Free allocated memory
    free(Q);
    Q = NULL;
    free(T);
    T = NULL;
    free(w);
    w = NULL;
    free(temp);
    temp = NULL;

    return sqrtA;
}


// Function to compute the product of three matrices
lapack_complex_double* multiply_three_matrices_LAPACK(int d, lapack_complex_double* A, lapack_complex_double* B, lapack_complex_double* C)
{
    printf("Entered multiply_three_matrices_LAPACK...\n");
    lapack_complex_double* temp1 = calloc(d*d, sizeof(lapack_complex_double));
    lapack_complex_double* result = calloc(d*d, sizeof(lapack_complex_double));

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, A, d, B, d, &zero, temp1, d);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, temp1, d, C, d, &zero, result, d);

    int psd = check_psd_trace_one(result, d);
    printf("check_psd_trace_one() on output of multiply_three_matrices, psd = %d.\n", psd);

    free(temp1);
    temp1 = NULL;

    return result;
}


// Function to compute trace
lapack_complex_double trace_LAPACK(int d, lapack_complex_double* A)
{
    lapack_complex_double tr = 0.0;
    for (int i = 0; i < d; i++)
        tr += A[i * d + i];
    return tr;
}

// Function to compute the trace of the product of two PSD matrices
double trace_product_LAPACK(int n, lapack_complex_double* A, lapack_complex_double* B)
{
    // Temporary storage for the product
    lapack_complex_double* temp = (lapack_complex_double*) calloc(n*n, sizeof(lapack_complex_double));
    if (!temp) {
      perror("calloc failed for temp");
      exit(EXIT_FAILURE);
    }

    // Compute the product: temp = A * B
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, A, n, B, n, &zero, temp, n);

    // Compute the trace
    double trace = 0.0;
    for(int i = 0; i < n; i++)
    {
        trace += temp[i * n + i];
    }

    // Free the temporary array
    free(temp);

    return trace;
}


// Checks if a matrix is normal.
bool is_normal(lapack_complex_double* A, int d)
{
    // Allocate memory for A*, AA*, and A*A
    lapack_complex_double* A_conj_transpose = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* AA_conj_transpose = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));
    lapack_complex_double* A_conj_transposeA = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));

    // Create A*
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            A_conj_transpose[j * d + i] = conj(A[i * d + j]);
        }
    }

    // Compute AA*
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, A, d, A_conj_transpose, d, &zero, AA_conj_transpose, d);

    // Compute A*A
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, A_conj_transpose, d, A, d, &zero, A_conj_transposeA, d);

    // Check if AA* = A*A
    bool normal = true;
    for (int i = 0; i < d * d; i++) {
        if (cabs(AA_conj_transpose[i] - A_conj_transposeA[i]) > 1e-10) {
            normal = false;
            break;
        }
    }

    // Free allocated memory
    free(A_conj_transpose);
    A_conj_transpose = NULL;
    free(AA_conj_transpose);
    AA_conj_transpose = NULL;
    free(A_conj_transposeA);
    A_conj_transposeA = NULL;

    return normal;
}


double fidelity_mixed(int d, lapack_complex_double* rho1, lapack_complex_double* rho2)
{
     // printf("Entered func fidelity_mixed...\n");

     // for (int j = 0; j < d*d; j++) {
     //     printf("rho1[%d]=(%lf, %lf)\n", j, creal(rho1[j]), cimag(rho1[j]));
     // }


    // Compute the square root of rho1
    lapack_complex_double* sqrt_rho1;
    sqrt_rho1 = sqrtm_LAPACK(d, rho1);

    // Allocate memory for the eigenvalues.
    double* w = (double*) calloc(d, sizeof(double));

    lapack_complex_double* U = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    lapack_complex_double* D_sqrt = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    lapack_complex_double* temp1 = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    lapack_complex_double* temp2 = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    lapack_complex_double* temp = (lapack_complex_double*) calloc(d*d, sizeof(lapack_complex_double));

    // Copy rho2 to U for sending to zheev.
    cblas_zcopy(d*d, rho2, 1, U, 1);

    // Perform the eigenvalue decomposition of rho2.
    int info1 = LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', d, U, d, w);

    // Check for success.
    if (info1 != 0) {
        fprintf(stderr, "LAPACKE_zheevd failed with info = %d\n", info1);
        exit(EXIT_FAILURE);
    }

    // Calculate the square root matrix.
    for (int i = 0; i < d; i++) {
      D_sqrt[i+i*d] = sqrt(w[i]);
    }

    free(w);
    w = NULL;

   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, D_sqrt, d, U, d, &zero, temp1, d);

   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, temp1, d, sqrt_rho1, d, &zero, temp2, d);

   cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, d, d, d, &one, temp2, d, temp2, d, &zero, temp, d);

   // int psd = check_psd_trace_one(temp, d);
   // printf("X^dagger X gives PSD = %d\n", psd);


   // Allocate memory for the eigenvalues.
   double* w1 = (double*) calloc(d, sizeof(double));

   // Perform the eigenvalue decomposition of temp = temp2^\dagger*temp2
   int info = LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', d, temp, d, w1);

   // Check for success.
   if (info != 0) {
       fprintf(stderr, "LAPACKE_zheevd failed with info = %d\n", info);
       exit(EXIT_FAILURE);
   }


   double fid = 0.0;
   for (int i = 0; i < d; i++) {
     if(fabs(w1[i]) < EPSILON) w1[i] = 0.0;
     fid += sqrt(w1[i]);
    // printf("fidelity = %lf\n", fid);
   }

   free(w1);
   w1 = NULL;

   fid *= fid;

  // Free the temporary arrays
    free(sqrt_rho1);
    sqrt_rho1 = NULL;
    free(U);
    U = NULL;
    free(D_sqrt);
    D_sqrt = NULL;
    free(temp1);
    temp1 = NULL;
    free(temp2);
    temp2 = NULL;
    free(temp);
    temp = NULL;

  return fid;
}

// Fidelity for mixed states using a simple method
double fidelity_mixed_simple(int d, lapack_complex_double* rho1, lapack_complex_double* rho2)
{
    // Allocate memory for product matrix
    lapack_complex_double* product = (lapack_complex_double*) calloc(d * d, sizeof(lapack_complex_double));

    // Compute the product rho1 * rho2
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, rho1, d, rho2, d, &zero, product, d);

    // Allocate memory for eigenvalues
    lapack_complex_double* eigenvalues = (lapack_complex_double*) calloc(d, sizeof(lapack_complex_double));

    // Perform eigenvalue decomposition
    // Initialize variables
    lapack_int info;
    lapack_complex_double* vl = NULL;  // Don't compute left eigenvectors
    lapack_complex_double* vr = NULL;  // Don't compute right eigenvectors

    // Compute eigenvalues
    info = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', d, product, d, eigenvalues, vl, d, vr, d);

    // Check for success
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zgeev failed with info = %d\n", info);
        exit(EXIT_FAILURE);
    }

    // lapack_int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', d, product, d, eigenvalues);
    // if (info != 0) {
    //     fprintf(stderr, "LAPACKE_zheev failed with info = %d\n", info);
    //     exit(EXIT_FAILURE);
    // }

    // Sum up square roots of eigenvalues
    double sum_sqrt_eigenvalues = 0.0;
    for (int i = 0; i < d; i++) {
        if (creal(eigenvalues[i]) > 0) {  // Numerical errors might lead to tiny negative eigenvalues
            sum_sqrt_eigenvalues += sqrt(creal(eigenvalues[i]));
        }
    }

    // Free allocated memory
    free(product);
    free(eigenvalues);

    return sum_sqrt_eigenvalues;
}


double fidelity_pure(int d, lapack_complex_double* psi1, lapack_complex_double* psi2)
{
  lapack_complex_double inner_product = 0.0;

  // Compute the inner product <psi1|psi2>
  for (int i = 0; i < d; i++) {
    inner_product += conj(psi1[i]) * psi2[i];
  }

  // Compute the squared magnitude of the inner product
  double fidelity = creal(inner_product) * creal(inner_product) +
                    cimag(inner_product) * cimag(inner_product);

  return fidelity;
}

/************************************************************************************/
// Functions used for analysis with Pretty good measurement. //
/************************************************************************************/

/* Function to compute the square root of a positive semidefinite Hermitian matrix. */
double complex* sqrtm(double complex* A)
{
    double complex* sqrtA = (double complex*) calloc(d*d, sizeof(double complex));
    if (!sqrtA) {
      perror("calloc failed for sqrtA");
      exit(EXIT_FAILURE);
    }

    // Allocate arrays for eigenvalues and eigenvectors
    double* w = (double*) calloc(d, sizeof(double));
    if (!w) {
      perror("calloc failed for w");
      exit(EXIT_FAILURE);
    }
    double complex* eigvec = (double complex*) calloc(d*d, sizeof(double complex));
    if (!eigvec) {
      perror("calloc failed for eigvec");
      exit(EXIT_FAILURE);
    }
    double complex* eigvec_rescale = (double complex*) calloc(d*d, sizeof(double complex));
    if (!eigvec_rescale) {
      perror("calloc failed for eigvec_rescale");
      exit(EXIT_FAILURE);
    }

    // Copy A into eigvec, as LAPACK routine will overwrite input
    for(size_t i = 0; i < d*d; i++) {
        eigvec[i] = A[i];
    }

    // Perform the eigenvalue decomposition
    lapack_int info = LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', d, eigvec, d, w);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zheevd failed with info = %d\n", info);
        exit(EXIT_FAILURE);
    }


    // Take square root of eigenvalues
    for(int i = 0; i < d; i++) {

        printf("w[%d]=%f\n", i, w[i]);
        w[i] = sqrt(w[i]);
        printf("sqrt(w[%d)=%f\n", i, w[i]);

        for (int j = 0; j < d; j++) {
          eigvec_rescale[i*d+j] = w[i]*eigvec[i*d+j];
        }

    }

    // Reconstruct sqrtA = eigvec_rescale * eigvec^H
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d,
                &one, eigvec_rescale, d, eigvec, d, &zero, sqrtA, d);

    // Free allocated arrays
    free(w);
    free(eigvec);
    free(eigvec_rescale);

    return sqrtA;
}

double complex* pseudoinv(size_t n, double complex* A)
{
    double complex* A_loc = (double complex*) calloc(d*d, sizeof(double complex));
    if (!A_loc) {
      perror("calloc failed for A_loc");
      exit(EXIT_FAILURE);
    }

    double complex* pinvA = (double complex*) calloc(n*n, sizeof(double complex));
    if (!pinvA) {
      perror("calloc failed for pinvA");
      exit(EXIT_FAILURE);
    }

    // Allocate arrays for eigenvalues and eigenvectors
    double* w = (double*) calloc(n, sizeof(double));
    if (!w) {
      perror("calloc failed for w");
      exit(EXIT_FAILURE);
    }
    double complex* eigvec = (double complex*) calloc(n*n, sizeof(double complex));
    if (!eigvec) {
      perror("calloc failed for eigvec");
      exit(EXIT_FAILURE);
    }

    // Copy A into eigvec, as LAPACK routine will overwrite input
    for(size_t i = 0; i < n * n; i++)
    {
        eigvec[i] = A[i];
    }

    // Perform the eigenvalue decomposition
    LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', n, eigvec, n, w);

    // Compute reciprocal of eigenvalues for non-zero values
    for(size_t i = 0; i < n; ++i) {
        if(w[i] > 1e-10) { // tolerance to check for zero, adjust as needed
            w[i] = 1.0 / w[i];
        } else {
            w[i] = 0.0;
        }
    }

    // Reconstruct pinvA = eigvec * w * eigvec'
    // Compute pinvA = eigvec * w
    for(size_t j = 0; j < n; ++j) {
        for(size_t i = 0; i < n; ++i) {
            pinvA[j*n + i] = eigvec[j*n + i] * w[i];
        }
    }
    // Compute pinvA = pinvA * eigvec'
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, n, n, n,
                &one, pinvA, n, eigvec, n, &zero, A_loc, n);

    // Copy A back to pinvA
    for(size_t i = 0; i < n * n; ++i) {
        pinvA[i] = A_loc[i];
    }

    // Free allocated arrays
    free(w);
    free(eigvec);
    free(A_loc);

    return pinvA;
}

/* Function to multiply two matrices: C = A * B. */
double complex* multiply_two_matrices(double complex* A, double complex* B)
{
    // Variables for interfacing with CBLAS (note the factor of 2 for complex)
    int lda = d, ldb = d, ldc = d;
    double complex* C = (double complex*) calloc(d*d, sizeof(double complex));
    if (!C) {
      perror("calloc failed for C");
      exit(EXIT_FAILURE);
    }

    // Compute C = alpha*A*B + beta*C
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, A, lda, B, ldb, &zero, C, ldc);

    return C;
}

/* Function to multiply three matrices: D = A * B * C. */
double complex* multiply_three_matrices(int n, double complex* A, double complex* B, double complex* C)
{
    double complex* res = (double complex*) calloc(n*n, sizeof(double complex));
    if (!res) {
      perror("calloc failed for res");
      exit(EXIT_FAILURE);
    }

    // Allocate an array for the intermediate result
    double complex* temp = (double complex*) calloc(n*n, sizeof(double complex));
    if (!temp) {
      perror("calloc failed for temp");
      exit(EXIT_FAILURE);
    }

    // Compute temp = A * B
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, A, n, B, n, &zero, temp, n);

    // Compute res = temp * C
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, temp, n, C, n, &zero, res, n);

    // Free allocated array
    free(temp);

    return res;
}

/* Function to generate a set of L random unitary matrices of dimension d. */
double complex**  generate_random_unitary_matrices(int num)
{
    // Allocate arrays for QR decomposition
    double complex* tau = (double complex*) calloc(d, sizeof(double complex));
    if (!tau) {
      perror("calloc failed for tau");
      exit(EXIT_FAILURE);
    }
    double complex** U = (double complex**) calloc(L, sizeof(double complex *));
    if (!U) {
      perror("calloc failed for U");
      exit(EXIT_FAILURE);
    }

    // For each matrix
    for(int l = 0; l < num; l++)
    {
        U[l] = (double complex*) calloc(d*d, sizeof(double complex));
        if (!U[l]) {
          perror("calloc failed for U[l]");
          exit(EXIT_FAILURE);
        }

        // Generate a random complex matrix
        for(int i = 0; i < d * d; i++) {
            double real_part = genrand_real1() * 2.0 - 1.0;
            double imag_part = genrand_real1() * 2.0 - 1.0;
            U[l][i] = real_part + imag_part * I;
        }

        // Perform QR decomposition
        LAPACKE_zgeqrf(LAPACK_COL_MAJOR, d, d, U[l], d, tau);

        // Generate unitary matrix
        LAPACKE_zungqr(LAPACK_COL_MAJOR, d, d, d, U[l], d, tau);
    }

    // Free allocated arrays
    free(tau);
    return U;

}

double complex* generate_random_unitary_matrix()
{
  // Allocate a matrix and fill it with random complex numbers
    double complex* A = (double complex*) calloc(d * d, sizeof(double complex));
    if (!A) {
      perror("calloc failed for A");
      exit(EXIT_FAILURE);
    }

    for(int i = 0; i < d * d; i++)
    {
        double real_part = genrand_real1() * 2.0 - 1.0;  // Random real part between -1 and 1
        double imag_part = genrand_real1() * 2.0 - 1.0;  // Random imag part between -1 and 1
        A[i] = real_part + imag_part * I;
    }

    // Perform QR decomposition
    double complex* tau = (double complex*)calloc(d, sizeof(double complex));
    if (!tau) {
      perror("calloc failed for tau");
      exit(EXIT_FAILURE);
    }
    LAPACKE_zgeqrf(LAPACK_COL_MAJOR, d, d, A, d, tau);

    // Convert Q to a full matrix
    LAPACKE_zungqr(LAPACK_COL_MAJOR, d, d, d, A, d, tau);

    // Free the tau array
    free(tau);

    // Return the Q part of the QR decomposition, which is a unitary matrix
    return A;
}

/* Function to compute the sum U_k * A * U_k^H for a set of unitary matrices U_k and a matrix A. */
double complex* apply_unital_channel(int L, int d, double complex** U, double complex* A)
{
    double complex* result = (double complex*) calloc(d*d, sizeof(double complex));
    if (!result) {
      perror("calloc failed for result");
      exit(EXIT_FAILURE);
    }

    // Temporary storage for the intermediate results
    double complex* temp1 = (double complex*) calloc(d * d, sizeof(double complex));
    if (!temp1) {
      perror("calloc failed for temp1");
      exit(EXIT_FAILURE);
    }
    double complex* temp2 = (double complex*) calloc(d * d, sizeof(double complex));
    if (!temp2) {
      perror("calloc failed for temp2");
      exit(EXIT_FAILURE);
    }

    // Start with zero result
    for(int i = 0; i < d * d; ++i)
    {
        result[i] = 0;
    }

    // For each matrix U_k
    for(int k = 0; k < L; ++k)
    {
        // Compute temp1 = U_k * A
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &one, U[k], d, A, d, &zero, temp1, d);

        // Compute temp2 = temp1 * U_k^H
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &one, temp1, d, U[k], d, &zero, temp2, d);

        // Add temp2 to the result
        for(int i = 0; i < d * d; ++i)
        {
            result[i] += temp2[i];
        }
    }

    // Free the temporary arrays
    free(temp1);
    free(temp2);

    return result;
}

/* Function to compute the trace of the product of two matrices. */
double trace_product(int n, double complex* A, double complex* B)
{
    // Temporary storage for the product
    double complex* temp = (double complex*) calloc(n*n, sizeof(double complex));
    if (!temp) {
      perror("calloc failed for temp");
      exit(EXIT_FAILURE);
    }

    // Compute the product: temp = A * B
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, A, n, B, n, &zero, temp, n);

    // Compute the trace
    double trace = 0.0;
    for(int i = 0; i < n; i++)
    {
        trace += temp[i * n + i];
    }

    // Free the temporary array
    free(temp);

    return trace;
}

/* Function to compute the Bayes mean: the sum p(j|x) * rho_j for a set of probabilities p_j, density matrices rho_j, and a measurement matrix M_x. */
double complex* compute_Bayes_mean_sum(int x, double complex* M_x, double complex** rho, double* p)
{
    double complex* bme_x = (double complex*) calloc(d*d, sizeof(double complex));
    if (!bme_x) {
      perror("calloc failed for bme_x");
      exit(EXIT_FAILURE);
    }

    // Temporary storage for the intermediate results
    double complex* temp1 = (double complex*) malloc(d * d * sizeof(double complex));
    if (!temp1) {
      perror("calloc failed for temp1");
      exit(EXIT_FAILURE);
    }
    double complex* temp2 = (double complex*) malloc(d* d * sizeof(double complex));
    if (!temp2) {
      perror("calloc failed for temp2");
      exit(EXIT_FAILURE);
    }
    double p_x = 0.0;

    // Compute p(x) = sum_j p_j * tr(M_x * rho_j)
    for(int j = 0; j < L; j++)
    {
        p_x += (p[j] * trace_product(d, M_x, rho[j]));
    }

    // Compute the sum for each j
    for(int j = 0; j < L; j++)
    {
        // Compute p(j|x) = tr(M_x * rho_j) * p_j / p_x
        double p_j_x = trace_product(d, M_x, rho[j]) * p[j] / p_x;

        // Multiply rho_j by p(j|x) and store in temp1
        for(int i = 0; i < d*d; i++)
        {
            temp1[i] = p_j_x * rho[j][i];
        }

        // Add temp1 to the result
        if(j == 0)
        {
            // Copy temp1 to result[x] for the first iteration
            for(int i = 0; i < d * d; i++)
            {
                bme_x[i] = temp1[i];
            }
        }
        else
        {
            // Add temp1 to bme_x for subsequent iterations
            cblas_zaxpy(d * d, &one, temp1, 1, bme_x, 1);
        }
    }

    // Free the temporary arrays
    free(temp1);
    free(temp2);

    return bme_x;
}

/* Function to generate a set of random density matrices.
Assuming the functions generate_random_unitary_matrix, multiply_two_matrices, sqrtm and trace are defined as follows:

void generate_random_unitary_matrix(int d, double complex* U);
void multiply_two_matrices(int d, double complex* A, double complex* B, double complex* result);
void sqrtm(int d, double complex* A, double complex* sqrtA);
double complex trace(int d, double complex* A);

*/
double complex* generate_random_complex_matrix()
{
  double complex *mat = (double complex*) calloc(d*d, sizeof(double complex));
  if (!mat) {
    perror("calloc failed for *mat");
    exit(EXIT_FAILURE);
  }

    for(int j=0; j<d; j++)
    {
      for(int i=0; i<d; i++)
        {
            double real_part = genrand_real1();
            double imag_part = genrand_real1();
            mat[j*d + i] = real_part + imag_part * I;
        }
    }
    return mat;
}


void generate_random_density_matrices(int L, double complex** rho)
{

    // Generate a local random matrix.
    double complex *C_loc;

    // For each density matrix
    for(int l = 0; l < L; l++)
    {
        // Allocate memory for a random matrix
        C_loc = (double complex*) calloc(d*d, sizeof(double complex));
        if (!C_loc) {
          perror("calloc failed for C_loc");
          exit(EXIT_FAILURE);
        }
        C_loc = generate_random_complex_matrix();

        // Compute the product: rho_l = C * C^H
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &one, C_loc, d, C_loc, d, &zero, rho[l], d);

        // Normalize the density matrix by its trace
        double complex tr = trace(rho[l]);
        for(int i = 0; i < d * d; i++) {
            rho[l][i] /= tr;
        }

        // Free the temporary array
        free(C_loc);
    }
}

/* Function to  compute the trace of a matrix. */
double complex trace(double complex* A)
{
    double complex trace = 0.0;
    for(int i = 0; i < d; i++)
    {
        trace += A[i * d + i];
    }
    return trace;
}

/* Function to compute the fidelity between two quantum states.
Assuming the functions sqrtm and multiply_three_matrices are defined as follows:

void sqrtm(int n, double complex* A, double complex* sqrtA);
void multiply_three_matrices(int n, double complex* A, double complex* B, double complex* C, double complex* result);

*/
double fidelity(int d, double complex* rho1, double complex* rho2)
{
    // Compute the square root of rho1
    double complex* sqrt_rho1 = (double complex*)malloc(d * d * sizeof(double complex));
    sqrt_rho1 = sqrtm(rho1);

    // Compute the product: temp = sqrt(rho1) * rho2 * sqrt(rho1)
    double complex* temp = (double complex*)malloc(d * d * sizeof(double complex));
    temp = multiply_three_matrices(d, sqrt_rho1, rho2, sqrt_rho1);

    // Compute the square root of temp
    double complex* sqrt_temp = (double complex*)malloc(d * d * sizeof(double complex));
    sqrt_temp = sqrtm(temp);

    // Compute the trace
    double complex tr = trace(sqrt_temp);

    // Compute the fidelity
    double fid = creal(tr) * creal(tr) + cimag(tr) * cimag(tr);

    // Free the temporary arrays
    free(sqrt_rho1);
    free(temp);
    free(sqrt_temp);

    return fid;
}

/* Function to compute the fidelities between a matrix A and a set of matrices B_x.
Assuming the functions sqrtm and multiply_three_matrices are defined as follows:

void sqrtm(int n, double complex* A, double complex* sqrtA);
void multiply_three_matrices(int n, double complex* A, double complex* B, double complex* C, double complex* result);

*/
double* compute_fidelities(double complex* A, double complex** B_x)
{
    double* fidelities = (double*) calloc(L, sizeof(double));
    if (!fidelities) {
      perror("calloc failed for fidelities");
      exit(EXIT_FAILURE);
    }

    // Compute the square root of A
    double complex* sqrt_A = (double complex*) calloc(d*d, sizeof(double complex));
    if (!sqrt_A) {
      perror("calloc failed for sqrt_A");
      exit(EXIT_FAILURE);
    }
    sqrt_A = sqrtm(A);

    //testing output of sqrtm()
    printf("Func: compute_fidelities. Printing output of sqrt of rho[0]\n");
    for (int j = 0; j < d*d; j++) {
      printf("(%lf, %lf)\n", creal(sqrt_A[j]), cimag(sqrt_A[j]));
    }

    // Temporary storage for the intermediate results
    double complex* temp = (double complex*) calloc(d*d, sizeof(double complex));
    if (!temp) {
      perror("calloc failed for temp");
      exit(EXIT_FAILURE);
    }
    double complex* sqrt_temp = (double complex*) calloc(d*d, sizeof(double complex));
    if (!sqrt_temp) {
      perror("calloc failed for sqrt_temp");
      exit(EXIT_FAILURE);
    }

    // For each matrix B_x
    for(int x = 0; x < L; x++) {
        // Compute the product: temp = sqrt(A) * B_x * sqrt(A)
        temp = multiply_three_matrices(d, sqrt_A, B_x[x], sqrt_A);

        // Compute the square root of temp
        sqrt_temp = sqrtm(temp);

        // Compute the trace
        double complex tr = trace(sqrt_temp);

        // testing output of trace()
        printf("Printing output of trace()\n");
        printf("(%lf, %lf)\n", creal(tr), cimag(tr));

        // Compute the fidelity
        fidelities[x] = creal(tr) * creal(tr) + cimag(tr) * cimag(tr);
    }

    // Free the temporary arrays
    free(sqrt_A);
    free(temp);
    free(sqrt_temp);

    return fidelities;
}


/* Function to compute rescaled matrices.*/
double complex** compute_p_rho(int num_matrices, int n, double* p, double complex** rho)
{
    double complex **p_rho = (double complex**)calloc(num_matrices, sizeof(double complex*));
    if (!p_rho) {
      perror("calloc failed for p_rho");
      exit(EXIT_FAILURE);
    }

    for (int j = 0; j < num_matrices; j++)
    {
        p_rho[j] = (double complex*) calloc(d*d, sizeof(double complex));
        if (!p_rho[j]) {
          perror("calloc failed for p_rho[j]");
          exit(EXIT_FAILURE);
        }
        // Compute p_j * rho_j and store it in p_rho_j
        for (int i = 0; i < n * n; i++)
        {
            p_rho[j][i] = p[j] * rho[j][i];
        }
    }
    return p_rho;
}

/* Function to compute the sum \sum_j p_j \rho_j. */
double complex * compute_sum_p_rho(int num_matrices, int n, double* p, double complex **rho)
{
    double complex *sum_p_rho;
    // Create an intermediate storage for p_j * rho_j
    double complex **p_rho = (double complex**)malloc(num_matrices * sizeof(double complex*));
    for(int i=0; i<num_matrices; i++)
        p_rho[i] = (double complex*)malloc(n * n * sizeof(double complex));

    // Compute p_j * rho_j for all j
    p_rho = compute_p_rho(num_matrices, n, p, rho);

    // Allocate and initialize sum_p_rho with zeros
    sum_p_rho = (double complex *)calloc(n * n, sizeof(double complex));
    if (!sum_p_rho) {
      perror("calloc failed for sum_p_rho");
      exit(EXIT_FAILURE);
    }

    // Compute the sum
    for (int j = 0; j < num_matrices; j++)
    {
        for (int i = 0; i < n * n; i++)
        {
            sum_p_rho[i] += p_rho[j][i];
        }
    }
    // Free the allocated memory
    for(int i=0; i<num_matrices; i++)
        free(p_rho[i]);
    free(p_rho);

    return sum_p_rho;
}
