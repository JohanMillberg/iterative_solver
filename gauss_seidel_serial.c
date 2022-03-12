#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

double** get_A_matrix(int amount_equations) {
    int i, j;
    double temp_sum;
    double diagonal_value;

    //Generate a random matrix
    double** A = (double**) malloc(amount_equations * sizeof(double*));
    for (i = 0; i < amount_equations; i++) {
        A[i] = malloc(amount_equations * sizeof(double));
        for (j = i; j < amount_equations; j++) {
            A[i][j] = (double)rand()/RAND_MAX*10.0;
            A[j][i] = A[j][i];
        }
    }

    //Ensure that matrix is diagonally dominant, necessary for GS.
    for (i = 0; i < amount_equations; i++) {
        diagonal_value = fabs(A[i][i]);
        temp_sum = 0;

        for (j = 0; j < amount_equations; j++) {
            temp_sum += fabs(A[i][j]);
        }

        A[i][i] += temp_sum;
    }

    return A;
}

int main(int argc, char *argv[]) {
    int i, j;
    int amount_equations = atoi(argv[1]);
    double start = omp_get_wtime();
    double** A = get_A_matrix(amount_equations);
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
        b[i] = (double)rand()/RAND_MAX*10.0;
    }

    double change = 0;
    double threshold = 0.00001;

    double* x_current = malloc(amount_equations*sizeof(double));
    for (int i = 0; i < amount_equations; i++) {
        x_current[i] = 0;
    }
    double* x_new = malloc((amount_equations/2)*sizeof(double));
    double temp_sum;
    double new_val;
    int counter;

    do {
        printf("Pre Jacobi: ");
        for (i = 0; i < amount_equations; i++) {
            printf("x_%d = %lf ", i, x_current[i]);
        }
        printf("\n");
        // The odd-indexed elements in x_current are updated using Jacobi
        change = 0;
        for (i = 1; i < amount_equations; i += 2) {
            temp_sum = 0;
            for (j = 0; j < amount_equations; j++) {
                if (i != j) {
                    temp_sum += (A[i][j]*x_current[j]);
                }
            }
            new_val = (1/A[i][i]) * (b[i]-temp_sum);
            change += (x_current[i] - new_val)*(x_current[i] - new_val);
            x_new[(i-1)/2] = new_val;
        }

        // Update the new values in x_current
        counter = 0;
        for (i = 1; i < amount_equations; i += 2) {
            x_current[i] = x_new[counter];
            counter++;
        }

        printf("Pre GS: ");
        for (i = 0; i < amount_equations; i++) {
            printf("x_%d = %lf ", i, x_current[i]);
        }
        printf("\n\n");
        // The even-indexed elements in x_current are updated using Gauss-Seidel
        for (i = 0; i < amount_equations; i += 2) {
            temp_sum = 0;
            for (j = 0; j < amount_equations; j++) {
                if (i != j) {
                    temp_sum += (A[i][j]*x_current[j]);
                }
            }
            new_val = (1/A[i][i]) * (b[i]-temp_sum);
            change = (x_current[i] - new_val)*(x_current[i] - new_val);
            x_current[i] = new_val;
        }

    } while (change > threshold);
    for (int i = 0; i < amount_equations; i++) {
        printf("x_%d = %lf\n", i, x_current[i]);
    }
    free(A);
    free(x_new);
    free(x_current);
    free(b);
    printf("Program took %lf seconds\n", omp_get_wtime() - start);
}