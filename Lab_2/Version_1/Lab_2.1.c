#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 100
#define t 10e-6
#define e 10e-9

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("invalid number of arguments\n");
        return 0;
    }

    int number_of_processes = atoi(argv[1]);

    omp_set_num_threads(number_of_processes);

    double *matrix = (double*)malloc(N*N * sizeof(double));
    double *b = (double*)malloc(N * sizeof(double));
    double *x = (double*)malloc(N * sizeof(double));
    double *tmp_x = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i) {
        b[i] = N + 1;
        x[i] = 0;

        for (int j = 0; j < N; ++j) {
            if (i == j) {
                matrix[i*N + j] = 2.0;
            } else {
                matrix[i*N + j] = 1.0;
            }
        }
    }

    double start_time = omp_get_wtime();

    double b_length = 0;

#pragma omp parallel for reduction(+:b_length)
    for (int i = 0; i < N; ++i) {
        b_length += b[i] * b[i];
    }

    b_length = sqrt(b_length);

    int is_calculation_complete = 1;

    while (is_calculation_complete) {
        double result_length = 0;

#pragma omp parallel for reduction(+:result_length)
        for (int i = 0; i < N; ++i) {
            double process_x = 0;
            for (int j = 0; j < N; ++j) {
                process_x += matrix[i * N + j] * x[j]; //Ax
            }
            process_x = process_x - b[i]; //Ax - b
            tmp_x[i] = x[i] - process_x * t; //x^n=x-t(Ax-b)

            result_length += process_x * process_x;
        }

        for (int i = 0; i < N; ++i) {
            x[i] = tmp_x[i];
        }

        result_length = sqrt(result_length);
        is_calculation_complete = (result_length / b_length >= e);
    }

    double finish_time = omp_get_wtime();

    printf("Processes: %d; Matrix Size: %dx%d; Time: %lf\n", number_of_processes, N, N, (finish_time - start_time));

    printf("Result: \n");
    for (int i = 0; i < N; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");

    free(matrix);
    free(b);
    free(x);
    free(tmp_x);

    return 0;
}