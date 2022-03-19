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
double** get_A_matrix(int amount_equations) {
    int i, j;
    double temp_sum;

    //Generate a random matrix
    double** A = (double**) malloc(amount_equations * sizeof(double*));
    for (i = 0; i < amount_equations; i++) {
        A[i] = malloc(amount_equations * sizeof(double));
    }

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

    return A;
}



int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Syntax to run program is: ./iterative_solver num_equations num_iterations\n");
        return 0;
    }
    int i, j;
    const int amount_equations = atoi(argv[1]);
    const int max_iterations = atoi(argv[2]);

    const int eq_per_loop = 4;

    if (amount_equations % eq_per_loop != 0) {
        printf("Loop treats %d equations per loop. Amount of equations need to be divisible with %d.\n", eq_per_loop, eq_per_loop);
        return 0;
    }

    double** A = get_A_matrix(amount_equations);

    /* Equations used for debugging below! Expected result is:
    x_0 = -0.25262
    x_1 = -0.097055
    x_2 = 0.13662
    x_3 = 0.73922
    */

/* 
    double** A = (double**) malloc(amount_equations * sizeof(double*));
    for (i = 0; i < amount_equations; i++) {
        A[i] = (double*) malloc(amount_equations * sizeof(double));
    }
    A[0][0] = 10; //x0
    A[0][1] = 2; //y0
    A[0][2] = -1; //z0
    A[0][3] = 4; //a0
    A[1][0] = -5; //x1
    A[1][1] = 15; //y1
    A[1][2] = -4; //z1
    A[1][3] = 1;//a1
    A[2][0] = 1; //x2
    A[2][1] = 2; //y2
    A[2][2] = 16; //z2
    A[2][3] = -1;// a2
    A[3][0] = 1;//x3
    A[3][1] = 1;//y3
    A[3][2] = -7;//z3
    A[3][3] = 18;//a3

    double* b = (double*) malloc(amount_equations * sizeof(double));
    b[0] = 0.1; //b0
    b[1] = 0; //b1
    b[2] = 1; //b2
    b[3] = 12;//b3
 */

    double* b = (double*) malloc(amount_equations * sizeof(double));
    for (i = 0; i < amount_equations; i++) {
        b[i] = (double)rand()/RAND_MAX*1000.0;
    }

    double* x_new = malloc(amount_equations * sizeof(double));
    double* x_current = malloc(amount_equations*sizeof(double));
    for (int i = 0; i < amount_equations; i++) {
        x_current[i] = 0;
        x_new[i] = 0;
    }

    double start = omp_get_wtime();
    int iter;
    for (iter = 0; iter < max_iterations; iter++){

        /*
        printf("Pre Jacobi: ");
        for (i = 0; i < amount_equations; i++) {
            printf("x_%d = %lf ", i, x_current[i]);
        }
        printf("\n");
        */

        // The odd-indexed elements in x_current are updated using Jacobi
        for (i = 1; i < amount_equations; i += 2) {
            double temp_sum = 0;
            for (j = 0; j < amount_equations; j++) { //Auto
                temp_sum = temp_sum + (A[i][j]*x_current[j]);
            }
            temp_sum = temp_sum - (A[i][i]*x_current[i]);
            double new_val = (1/A[i][i]) * (b[i]-temp_sum);
            x_new[i] = new_val;
        }

        /*
        printf("Pre GS: ");
        for (i = 0; i < amount_equations; i++) {
            printf("x_%d = %lf ", i, x_current[i]);
        }
        printf("\n\n");
        */

        for (i = 0; i < amount_equations; i += 2) {
            double temp_sum = 0;
            for (j = 0; j < amount_equations; j++) {
                temp_sum = temp_sum + (A[i][j]*x_new[j]);
            }
            temp_sum = temp_sum - (A[i][i]*x_new[i]);
            double new_val = (1/A[i][i]) * (b[i]-temp_sum);
            x_current[i] = new_val;
        }

        // Synchs the two arrays
        for (i = 1; i < amount_equations; i+=2) {
            x_current[i] = x_new[i];
            x_new[i-1] = x_current[i-1];
        }

    }

    //Prints all results for debugging purposes

/*
    for (int i = 0; i < amount_equations; i++) {
        printf("x_%d = %lf\n", i, x_current[i]);
    }
*/


    for (int i = 0; i < amount_equations; i++) {
        free(A[i]);
    }
    free(A);
    free(x_new);
    free(x_current);
    free(b);
    printf("Program took %lf seconds\n", omp_get_wtime() - start);
}