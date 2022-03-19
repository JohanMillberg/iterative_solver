#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

/*
The purpose of this program is to solve equations systems of N equations using n_threads
amount of threads. Both of these parameters are provided by the user as input when running
the program. The program then first utilizes the Jacobi method to update every odd-indexed element
in the solution vector. Then, the even-indexed elements of the solution vector are updated using
both the old values and the previously updated values. This algorithm is similar to the
Gauss Seidel algorithm, but this version is perfectly parallel.
*/

/**
 * @brief Generates a pseudo random diagonally dominant coefficient matrix
 *
 * @param amount_equations The amount of equations in the system
 * @return The complete generated matrix
 */
double** get_A_matrix(int amount_equations, int n_threads) {
    int i, j;
    double temp_sum;
    double start = omp_get_wtime();

    //Generate a random matrix
    double** A = (double**) malloc(amount_equations * sizeof(double*));
    for (i = 0; i < amount_equations; i++) {
        A[i] = malloc(amount_equations * sizeof(double));
    }
    #pragma omp parallel for private(i, j) num_threads(n_threads)
    for (i = 0; i < amount_equations; i++) {
        for (j = i; j < amount_equations; j++) {
            A[i][j] = (double)rand()/RAND_MAX*10.0;
            A[j][i] = A[i][j];
        }
    }

    //Ensure that matrix is diagonally dominant, necessary for GS.
    for (i = 0; i < amount_equations; i++) {
        temp_sum = 0;

        for (j = 0; j < amount_equations; j++) {
            temp_sum = temp_sum + fabs(A[i][j]);
        }

        A[i][i] += temp_sum;
    }
    printf("Creation of matrix took %lf seconds.\n", omp_get_wtime()-start);
    return A;
}

/**
 * @brief Function used to verify the results of the iterative algorithm
 *
 * @param A The coefficient matrix
 * @param b The right hand side vector of the original equation
 * @param x The result vector calculated using the iterative algorithm
 * @param amount_equations The amount of equations/variables
 */
void verify_result(double** A, double* b, double* x, int amount_equations) {
    int i, j;
    double error_sum = 0;
    double* result = malloc(amount_equations * sizeof(double));

    double temp;

    /*
    Multiplies A with the result vector x, and checks the difference
    between the product and the original right hand side b
    */
    for (i = 0; i  < amount_equations; i++) {
        temp = 0;
        for (j = 0; j < amount_equations; j++) {
            temp += A[i][j] * x[j];
        }

        error_sum += fabs(temp - b[i]);
    }

    printf("Sum of error: %lf\n", error_sum);
    free(result);
}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Syntax to run program is: ./iterative_solver num_equations num_threads num_iterations\n");
        return 0;
    }
    int i, j;
    const int amount_equations = atoi(argv[1]);
    const int n_threads = atoi(argv[2]);
    const int max_iterations = atoi(argv[3]);

    const int eq_per_loop = 8;

    if (amount_equations % eq_per_loop != 0) {
        printf("Loop handles %d equations per loop. Amount of equations need to be divisible with %d.\n", eq_per_loop, eq_per_loop);
        return 0;
    }


    double** A = get_A_matrix(amount_equations, n_threads);

    double* b = (double*) malloc(amount_equations * sizeof(double));

    for (i = 0; i < amount_equations; i++) {
        b[i] = (double)rand()/RAND_MAX*1000.0;
    }

    double* x_new = malloc(amount_equations * sizeof(double));
    double* x_current = malloc(amount_equations*sizeof(double));
    for (i = 0; i < amount_equations; i++) {
        x_current[i] = 0;
        x_new[i] = 0;
    }

    double start = omp_get_wtime();
    int iter;
    for (iter = 0; iter < max_iterations; iter++){

        // The odd-indexed elements in x_current are updated using Jacobi, utilizing loop unrolling
        #pragma omp parallel num_threads(n_threads)
        {
            #pragma omp for private(i, j)
            for (i = 1; i < amount_equations; i += 2) {
                double temp_sum_1 = 0;


                for (j = 0; j < amount_equations; j++) {
                    temp_sum_1 = temp_sum_1 + (A[i][j]*x_current[j]);
                }

                temp_sum_1 = temp_sum_1 - (A[i][i]*x_current[i]);

                double new_val_1 = (1/A[i][i]) * (b[i]-temp_sum_1);

                x_new[i] = new_val_1;
            }

        // The even-indexed elements in x_current are updated, utilizing loop unrolling
            #pragma omp for private(i, j)
            for (i = 0; i < amount_equations; i += 2) {
                double temp_sum_1 = 0;

                for (j = 0; j < amount_equations; j++) {
                    temp_sum_1 = temp_sum_1 + (A[i][j]*x_new[j]);
                }
                temp_sum_1 = temp_sum_1 - (A[i][i]*x_new[i]);

                double new_val_1 = (1/A[i][i]) * (b[i]-temp_sum_1);

                x_current[i] = new_val_1;
            }
        }

        // Synchs the two arrays
        for (i = 1; i < amount_equations; i+=2) {
            x_current[i] = x_new[i];
            x_new[i-1] = x_current[i-1];
        }

    }

    //Uncomment the line below in order to call the verification function
    //verify_result(A, b, x_new, amount_equations);

    for (i = 0; i < amount_equations; i++) {
        free(A[i]);
    }
    free(A);
    free(x_new);
    free(x_current);
    free(b);
    printf("Calculations took %lf seconds\n", omp_get_wtime() - start);

    return 0;
}