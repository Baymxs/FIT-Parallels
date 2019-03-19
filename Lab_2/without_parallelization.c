#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000
#define t 10e-5
#define e 10e-8

double vectorLength(double *vector) {
    double result = 0;
    int i;

    for (i = 0; i < N; i++) {
        result += vector[i] * vector[i];
    }

    return sqrt(result);
}

double* matrixAndVectorMultiplication(const double *matrix, const double *vector) {
    double *result = (double*)malloc(N * sizeof(double));
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }

    return result;
}

double* vectorSubtraction(const double *vector_1, const double *vector_2) {
    double *result = (double*)malloc(N * sizeof(double));
    int i;

    for (i = 0; i < N; i++) {
        result[i] = vector_1[i] - vector_2[i];
    }

    return result;
}

double* vectorAndScalarMultiplication(const double *vector, double scalar) {
    double *result = (double*)malloc(N * sizeof(double));
    int i;

    for (i = 0; i < N; i++) {
        result[i] = vector[i] * scalar;
    }

    return result;
}

int main() {
    double *matrix = (double*)malloc(N*N * sizeof(double));
    double *b = (double*)malloc(N * sizeof(double));
    double *result = (double*)malloc(N * sizeof(double));

    double *temp = NULL;

    int i, j;

    for (i = 0; i < N; i++) {
        b[i] = N + 1;
        result[i] = 0;

         for (j = 0; j < N; j++) {
             if (i == j) {
                 matrix[i*N + j] = 2.0;
             } else {
                 matrix[i*N + j] = 1.0;
             }
         }
    }

    int is_calculation_complete = 1;

    while (is_calculation_complete) {
        temp = vectorSubtraction(matrixAndVectorMultiplication(matrix, result), b); //Ax-b
        if (vectorLength(temp)/vectorLength(b) < e) {
            is_calculation_complete = 0;
        } else {
            result = vectorSubtraction(result, vectorAndScalarMultiplication(temp, t)); //x-t(Ax-b)
        }
    }

    free(matrix);
    free(b);
    free(result);

    return 0;
}