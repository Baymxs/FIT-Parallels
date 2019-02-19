#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 10
#define t 10e-6
#define e 10e-9

int main(int argc, char *argv[]){
    int number_of_processes, current_process;
    /*
     Любая прикладная MPI-программа (приложение) должна начинаться с вызова функции инициализации MPI: функции MPI_Init.
     В результате выполнения этой функции создается группа процессов, в которую помещаются все процессы приложения,
     и создается область связи, описываемая предопределенным коммуникатором MPI_COMM_WORLD. Эта область связи
     объединяет все процессы-приложения. Процессы в группе упорядочены и пронумерованы от 0 до groupsize-1,
     где groupsize равно числу процессов в группе. Кроме этого, создается предопределенный коммуникатор MPI_COMM_SELF,
     описывающий свою область связи для каждого отдельного процесса
     В программах на C каждому процессу при инициализации передаются аргументы функции main, полученные из командной строки
    */
    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    //printf("Number of processes: %d; Current process: %d\n", number_of_processes, current_process);

    int first_process_line;
    int i, j;

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

    //printf("Lines in this process: %d\n", lines_per_process[current_process]);
    //printf("First process line: %d\n\n", first_process_line);


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


    double *result = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        result[i] = 0;
    }

    double start_time, end_time;
    double b_length = 0;
    int is_calculation_complete = 1;

    if (current_process == 0) {
        start_time = MPI_Wtime();

        for (i = 0; i < N; i++) {
            b_length += (N + 1) * (N + 1);
        }
        b_length = sqrt(b_length);

//        printf("Vector b length: %f\n\n", b_length);
    }

    double *process_result = (double*)malloc(lines_per_process[current_process] * sizeof(double));
    for (int i = 0; i < lines_per_process[current_process]; i++) {
        process_result[i] = 0;
    }

    while (is_calculation_complete) {
        double process_length = 0;
        for (i = 0; i < lines_per_process[current_process]; i++) {
            double temp_sum = 0;
            for (j = 0; j < N; j++) {
                temp_sum += matrix[i*N + j] * result[j]; //Ax
            }
            temp_sum = temp_sum - b[i]; //Ax-b

            process_result[i] = result[i + first_process_line] - t * temp_sum; //x-t(Ax-b)
            process_length += temp_sum * temp_sum; //part of result length
        }

        if(current_process != 0) {
            /* int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
               Функция выполняет посылку count элементов типа datatype сообщения с идентификатором tag процессу dest в области связи коммуникатора comm
            */
            MPI_Send(process_result, lines_per_process[current_process], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&process_length, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        } else {
            for (i = 0; i < lines_per_process[current_process]; i++) {
                result[i] = process_result[i];
            }

            double result_length = process_length;

            int current_line = lines_per_process[current_process];

            //Заполнение результируещего вектора посредством получения результатных строк из других процессов
            for (i = 1; i < number_of_processes; i++) {
                /*
                 int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
                 Функция выполняет прием count элементов типа datatype сообщения с идентификатором tag от процесса source в области связи коммуникатора comm.
                 */
                MPI_Recv(&result[current_line], lines_per_process[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                current_line += lines_per_process[i];

                double tmp_length;
                MPI_Recv(&tmp_length, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                result_length += tmp_length;
            }
            result_length = sqrt(result_length);

            is_calculation_complete = (result_length / b_length >= e);
        }
        /*
         int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
         Широковещательная рассылка данных выполняется с помощью функции MPI_Bcast. Процесс с номером root рассылает сообщение из своего буфера передачи всем процессам области связи коммуникатора comm.
        */
        MPI_Bcast(&is_calculation_complete, 1, MPI_INT, 0, MPI_COMM_WORLD);//what,how many, type,root,communicator
        MPI_Bcast(result, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (current_process == 0) {
        end_time = MPI_Wtime();

        printf("Processes: %d; Time: %f\n", number_of_processes, end_time - start_time);

        printf("Result: ");
        for (i = 0; i < N; i++) {
            printf("%f ", result[i]);
        }
        printf("\n");
    }

    free(lines_per_process);
    free(matrix);
    free(b);
    free(result);
    free(process_result);

    /* Функция закрывает все MPI-процессы и ликвидирует все области связи */
    MPI_Finalize();
    return 0;
}
