#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 10000
#define t 10e-6
#define e 10e-9

int main(int argc, char *argv[]){
    int number_of_processes, current_process;

    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    //printf("Number of processes: %d; Current process: %d\n", number_of_processes, current_process);

    int i, j;

    int first_process_line;

    //Массив для хранения кол-ва строк в i-ом процессе
    int *lines_per_process = (int*)malloc(number_of_processes * sizeof(int));

    /*
     temp - количество процессов, в которых будет идти обработка N / numberOfProcesses строк
     Тогда numberOfProcesses - temp - это кол-во процессов, в которых количество строк будет увеличиваться на 1,
     чтобы покрыть все количество строк в матрице
    */
    int temp = number_of_processes - (N % number_of_processes);
    for (i = 0; i < number_of_processes; i++) {
        if (i < temp) {
            lines_per_process[i] = N / number_of_processes;
        } else {
            lines_per_process[i] = N / number_of_processes + 1;
        }
        /* Находим первую строку матрицы для i процесса  */
        if(i < current_process) {
            first_process_line += lines_per_process[i];
        }
    }

//    printf("Lines in this process: %d\n", lines_per_process[current_process]);
//    printf("First process line: %d\n\n", first_process_line);


    double *matrix = (double*)malloc(lines_per_process[current_process]*N * sizeof(double));
    for(i = 0; i < lines_per_process[current_process]; i++) {
        for(j = 0; j < N; j++) {
            if(first_process_line + i == j){
                matrix[i*N + j] = 2;
            } else{
                matrix[i*N + j] = 1;
            }
        }
    }

    double *full_x = (double*)malloc(N * sizeof(double));
    for (i = 0; i < N; i++) {
        full_x[i] = 0;
    }

    double *part_x = (double*)malloc(lines_per_process[current_process] * sizeof(double));
    for (i = 0; i < lines_per_process[current_process]; i++) {
        part_x[i] = 0;
    }

    //Рассылаем часть вектора x остальным процессам для сборки у них полного вектора x
    for (i = 0; i < number_of_processes; i++) {
        if (i != current_process) {
            MPI_Send(part_x, lines_per_process[current_process], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }


//    printf("Matrix:\n");
//    for (i = 0; i < lines_per_process[current_process]; i++) {
//        for (j = 0; j < N; j++) {
//            printf("%f ", matrix[i*N + j]);
//        }
//        printf("\n");
//    }
//    printf("\n");

    double *b = (double*)malloc(lines_per_process[current_process] * sizeof(double));
    for(i = 0; i < lines_per_process[current_process]; i++) {
        b[i] = N + 1;
    }

//    printf("Vector b:\n");
//    for (i = 0; i < lines_per_process[current_process]; i++) {
//        printf("%f ", b[i]);
//    }
//    printf("\n");


    double start_time, end_time;
    double b_length = 0;
    int is_calculation_complete = 1;

    if (current_process == 0) {
        start_time = MPI_Wtime();

        for (i = 0; i < N; i++) {
            b_length += (N + 1) * (N + 1);
        }
        b_length = sqrt(b_length);

        //printf("Vector b length: %f\n\n", b_length);
    }
    while (is_calculation_complete) {
        /*
         * Сборка full_x
         */

        //Записываем в full_x part_x current процесса
        memcpy(full_x + first_process_line, part_x, lines_per_process[current_process] * sizeof(double));

        //Записываем в full_x part_x остальных процессов (!= current_process)
        int current_line = 0;
        for (i = 0; i < number_of_processes; i++) {
            if (i != current_process) {
                MPI_Recv(&full_x[current_line], lines_per_process[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            current_line += lines_per_process[i];
        }

        double process_length = 0;
        for (i = 0; i < lines_per_process[current_process]; i++) {
            double temp_sum = 0;
            for (j = 0; j < N; j++) {
                temp_sum += matrix[i * N + j] * full_x[j]; //Ax
            }
            temp_sum = temp_sum - b[i]; //Ax-b

            part_x[i] = part_x[i] - t * temp_sum; //x-t(Ax-b)
            process_length += temp_sum * temp_sum; //part of result length
        }

        if (current_process!=0) {
            MPI_Send(&process_length, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }

        for (i = 0; i < number_of_processes; i++) {
            if (i != current_process) {
                MPI_Send(part_x, lines_per_process[current_process], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }

        if(current_process == 0) {
            double tmp_length = process_length;
            for (i = 1; i < number_of_processes; i++) {
                double tmp_tmp_length;
                MPI_Recv(&tmp_tmp_length, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                tmp_length += tmp_tmp_length;
            }
            //printf("Tmp length: %lf", tmp_length);
            tmp_length = sqrt(tmp_length);

            is_calculation_complete = (tmp_length / b_length >= e);
        }
        MPI_Bcast(&is_calculation_complete, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if (current_process == 0) {
        end_time = MPI_Wtime();

        memcpy(full_x + first_process_line, part_x, lines_per_process[current_process] * sizeof(double));
        int current_line = 0;
        for (i = 0; i < number_of_processes; i++) {
            if (i != current_process) {
                MPI_Recv(&full_x[current_line], lines_per_process[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            current_line += lines_per_process[i];
        }

        printf("Processes: %d; Time: %f\n", number_of_processes, end_time - start_time);

        printf("Result: ");
        for (i = 0; i < N; i++) {
            printf("%f ", full_x[i]);
        }
        printf("\n");
    }

    free(lines_per_process);
    free(matrix);
    free(b);
    free(full_x);
    free(part_x);

    /* Функция закрывает все MPI-процессы и ликвидирует все области связи */
    MPI_Finalize();
    return 0;
}
