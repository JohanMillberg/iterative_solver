#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>

typedef struct args {
    int start;
    int eq_per_thread;
    int n_equations;
    double** A;
    double* b;
    double* x_current;
    double* x_new;
} args_t;

void* jacobi_threaded(void* arg) {
    int i, j, counter = 0;
    args_t* args = (args_t*)arg;
    int start = args->start;
    int n_equations = args->n_equations;
    double** A = args->A;
    double* b = args->b;
    double* x_current = args->x_current;
    double* x_new = args->x_new;
    double temp_sum;
    double new_val;

    for (i = start; i < start + (args->eq_per_thread); i+=2) {
        temp_sum = 0;
        for (j = 0; j < n_equations; j++) {
            if (i != j) {
                temp_sum += (A[i][j]*x_current[j]);
            }
        }
        new_val = (1/A[i][i]) * (b[i]-temp_sum);

        x_new[counter] = new_val;
        counter++;
    }
    return NULL;
}

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
    int n_threads = atoi(argv[2]);
    int equations_per_thread = amount_equations / n_threads;

    //double** A = get_A_matrix(amount_equations);

    double** A = (double**) malloc(amount_equations * sizeof(double*));
    for (i = 0; i < amount_equations; i++) {
        A[i] = (double*) malloc(amount_equations * sizeof(double));
    }
    A[0][0] = 16; //x0
    A[0][1] = 3; //y0
    //A[0][2] = 3; //z0
    A[1][0] = 7; //x1
    A[1][1] = -11; //y1
    //A[1][2] = -5; //z1
    //A[2][0] = 2; //x2
    //A[2][1] = -1; //y2
    //A[2][2] = 13; //z2

    double* b = (double*) malloc(amount_equations * sizeof(double));
    b[0] = 11; //b0
    b[1] = 13; //b1
    //b[2] = 8; //b2

    /*
    double* b = (double*) malloc(amount_equations * sizeof(double));
    for (i = 0; i < amount_equations; i++) {
        b[i] = (double)rand()/RAND_MAX*10.0;
    }
    */
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

    args_t** args = malloc(n_threads*sizeof(args_t*));
    for (i = 0; i < n_threads; i++) {
        args[i] = malloc(sizeof(args_t));
        args[i]->start = (i*equations_per_thread)+1;
        args[i]->eq_per_thread = equations_per_thread;
        args[i]->A = A;
        args[i]->b = b;
        args[i]->x_current = x_current;
        args[i]->x_new = malloc(equations_per_thread*sizeof(double));
    }

    pthread_t threads[n_threads];

    do {
        printf("Pre Jacobi: ");
        for (i = 0; i < amount_equations; i++) {
            printf("x_%d = %lf ", i, x_current[i]);
        }
        printf("\n");
        // The odd-indexed elements in x_current are updated using Jacobi
        for (i = 0; i < n_threads; i++) {
            pthread_create(&threads[i], NULL, jacobi_threaded, args[i]);
        }

        for (i = 0; i < n_threads; i++) {
            pthread_join(threads[i], NULL);
            //for (j = 1; j < equations_per_thread; j+=2) {
            //    x_current[i*equations_per_thread+j] =
            //}
        }

        counter = 0;
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
            x_new[counter] = new_val;
            counter++;
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
}