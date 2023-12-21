#include "global.h"

#define NUM_CLIFFORDS 24

/* The following is copied from main_CLifford_vs_RMB_qubit.c for consistency. */
/*=============================================================================*/
/* Original declaration for extern variables. */
int d; /* dimension */
int L; /* no. of states in the ensemble */
int N0; /* sample complexity/no. of measurements */
const double zero_epsilon = 1e-6; /* Defining 0 for comparison with float types. */;
/*=============================================================================*/

// Hadamard
lapack_complex_double H[4] = {{0.707106781, 0}, {0.707106781, 0}, {0.707106781, 0}, {-0.707106781, 0}};

// Phase
lapack_complex_double S[4] = {{1, 0}, {0, 0}, {0, 0}, {0, 1}};

// Function to multiply two 2x2 matrices and store result in the first matrix
void multiply_matrices(lapack_complex_double* A, lapack_complex_double* B)
{
    lapack_complex_double result[4];

    // Column-major order multiplication
    result[0] = A[0] * B[0] + A[2] * B[1]; // Element at first row, first column
    result[1] = A[1] * B[0] + A[3] * B[1]; // Element at second row, first column
    result[2] = A[0] * B[2] + A[2] * B[3]; // Element at first row, second column
    result[3] = A[1] * B[2] + A[3] * B[3]; // Element at second row, second column

    memcpy(A, result, 4 * sizeof(lapack_complex_double));  // Copy result back into A
}

// Function to generate and store Clifford group matrices based on a sequence of H and S gates
void generate_clifford_matrices(lapack_complex_double clifford_matrices[NUM_CLIFFORDS][4], char sequences[NUM_CLIFFORDS][7])
{
    for (int i = 0; i < NUM_CLIFFORDS; i++) {
        // Start with the identity matrix
        clifford_matrices[i][0] = 1 + 0 * I;
        clifford_matrices[i][1] = 0 + 0 * I;
        clifford_matrices[i][2] = 0 + 0 * I;
        clifford_matrices[i][3] = 1 + 0 * I;

        // Iterate through the sequence of gates
        for (int j = 0; j < strlen(sequences[i]); j++) {
            switch (sequences[i][j]) {
                case 'H':
                          multiply_matrices(clifford_matrices[i], H);
                    break;
                case 'S':
                          multiply_matrices(clifford_matrices[i], S);
                    break;
            }
        }
    }
}

int main() {
    // Array of strings representing the sequences of gates
    char sequences[NUM_CLIFFORDS][7] = {
        "I", "H", "S", "HS", "SH", "SS", "HSH", "HSS", "SHS", "SSH", "SSS",
        "HSHS", "HSSH", "HSSS", "SHSS", "SSHS", "HSHSS", "HSSHS", "SHSSH", "SHSSS", "SSHSS",
        "HSHSSH", "HSHSSS", "HSSHSS"
    };

    // Array to store the generated Clifford matrices
    lapack_complex_double clifford_matrices[NUM_CLIFFORDS][4];

    // Generate the matrices
    generate_clifford_matrices(clifford_matrices, sequences);

    // Print the generated matrices
    for (int i = 0; i < NUM_CLIFFORDS; i++) {
        printf("Matrix for sequence %s:\n", sequences[i]);
        for (int j = 0; j < 4; j++) {
            printf("(%f, %f) ", creal(clifford_matrices[i][j]), cimag(clifford_matrices[i][j]));
            if ((j + 1) % 2 == 0) {
                printf("\n");
            }
        }
        printf("\n");
    }

    return 0;
}
