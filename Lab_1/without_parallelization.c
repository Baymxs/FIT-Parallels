#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000
#define METHOD_PARAMETER 10e-5
#define LESS_THAN 10e-8

double vectorLength(double *vector) {
    double result = 0;

    for (int i = 0; i < N; i++) {
        result += vector[i] * vector[i];
    }

    return sqrt(result);
}

double* matrixAndVectorMultiplication(const double *matrix, const double *vector) {
    double *result = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }

    return result;
}

double* vectorSubtraction(const double *vector_1, const double *vector_2) {
    double *result = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        result[i] = vector_1[i] - vector_2[i];
    }

    return result;
}

double* vectorAndScalarMultiplication(const double *vector, double scalar) {
    double *result = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        result[i] = vector[i] * scalar;
    }

    return result;
}

int main() {
    double *matrix = (double*)malloc(N*N * sizeof(double));
    double *b = (double*)malloc(N * sizeof(double));
    double *result = (double*)malloc(N * sizeof(double));

    double *temp = NULL;

    for (int i = 0; i < N; i++) {
        b[i] = N + 1;
        result[i] = 0;

         for (int j = 0; j < N; j++) {
             if (i == j) {
                 matrix[i*N + j] = 2.0;
             } else {
                 matrix[i*N + j] = 1.0;
             }
         }
    }

//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            printf("%f ", matrix[i*N + j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//
//    for (int i = 0; i < N; i++) {
//        printf("%f ", b[i]);
//    }
//    printf("\n");
//    printf("\n");
//
//    for (int i = 0; i < N; i++) {
//        printf("%f ", result[i]);
//    }
//    printf("\n");

    int is_calculation_complete = 1;

    while (is_calculation_complete) {
        temp = vectorSubtraction(matrixAndVectorMultiplication(matrix, result), b);
        if (vectorLength(temp)/vectorLength(b) < LESS_THAN) {
            is_calculation_complete = 0;
        } else {
            result = vectorSubtraction(result, vectorAndScalarMultiplication(temp, METHOD_PARAMETER));
        }
    }

//    for (int i = 0; i < N; i++) {
//        printf("%f ", result[i]);
//    }

    free(matrix);
    free(b);
    free(result);

    return 0;
}