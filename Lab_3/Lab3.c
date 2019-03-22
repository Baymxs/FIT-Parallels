#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

/* DIMS_NUMS - размер декартовой топологии. */
#define DIMS_NUM 2

#define P0 4
#define P1 4

#define N1 1600
#define N2 1600
#define N3 1600

#define N1 1600
#define N2 1600
#define N3 1600

int main(int argc, char *argv[]) {
    int number_of_processes, current_process; // кол-во процессов, номер нынешнего процесса

    int dims[DIMS_NUM]; // число процессов вдоль измерения
    int periods[DIMS_NUM]; // переодичность граничных условий
    int reorder = 0; // перенумерация процессов

    double n[3]; // размеры матриц: n[0] = N1, n[1] = N2, n[2] = N3
   // размер двумерной решетки: p[0] = dims[0] x p[1] = dims[1]

    double *A, *B, *C; // A,B - исходные матрицы, над которыми выполняется операция умножения
                       // С - результат перемножения A и B

    double *AA, *BB, *CC;  // AA - локальная подматрица матрицы А, состоящая из p[0] горизонтальных полос
                           // BB - локальная подматрица матрицы B, состоящая из p[1] вертикальных полос
                           // СС - локальная подматрица матрицы С

    int nn[2]; // nn[0] - размер подматрицы AA
               // nn[1] - размер подматрицы BB

    int coords[2]; // декартовы координаты данного процесса
    int current_grid_rank;

    MPI_Comm grid_comm; // коммуникатор для 2d решетки

    MPI_Init(&argc,&argv); //
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes); //Функция определения числа процессов в области связи MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process); //Функция определения номера процесса в области связи MPI_COMM_WORLD

    for (int i = 0; i < DIMS_NUM; i++) {
        dims[i] = 0;
        periods[i] = 0;
    }
    /* Функция определения оптимальной конфигурации сетки
     *
     * MPI_Dims_create(int nnodes, int ndims, int *dims)
     *
     * IN nnodes - общее число узлов в сетке;
     * IN ndims	- число измерений;
     * INOUT dims -	массив целого типа размерности ndims, в который помещается рекомендуемое число процессов вдоль каждого измерения.
     *
     * На входе в процедуру в массив dims должны быть занесены целые неотрицательные числа. Если элементу массива dims[i] присвоено положительное число,
     * то для этой размерности вычисление не производится (число процессов вдоль этого направления считается заданным).
     * Вычисляются только те компоненты dims[i], для которых перед обращением к процедуре были присвоены значения 0.
     * Функция стремится создать максимально равномерное распределение процессов вдоль направлений, выстраивая их по убыванию,
     * т.е. для 12-ти процессов она построит трехмерную сетку 4 х 3 х 1.
     * Результат работы этой процедуры может использоваться в качестве входного параметра для процедуры MPI_Cart_create.*/
 //   MPI_Dims_create(number_of_processes, DIMS_NUM, dims); // определение размеров решетки

    /* Функция создания коммуникатора с декартовой топологией
     *
     * MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods, int reorder, MPI_Comm *comm_cart)
     *
     *IN comm_old -	родительский коммуникатор;
     *IN ndims - число измерений;
     *IN dims - массив размера ndims, в котором задается число процессов вдоль каждого измерения;
     *IN periods - логический массив размера ndims для задания граничных условий (true - периодические, false - непериодические);
     *IN reorder - логическая переменная, указывает, производить перенумерацию процессов (true) или нет (false);
     *OUT comm_cart - новый коммуникатор.
     *Функция является коллективной, т.е. должна запускаться на всех процессах, входящих в группу коммуникатора comm_old.
     * При этом, если какие-то процессы не попадают в новую группу, то для них возвращается результат MPI_COMM_NULL.
     * В случае, когда размеры заказываемой сетки больше имеющегося в группе числа процессов, функция завершается аварийно.
     * Значение параметра reorder=false означает, что идентификаторы всех процессов в новой группе будут такими же, как в старой группе.
     * Если reorder=true, то MPI будет пытаться перенумеровать их с целью оптимизации коммуникаций.
     */
    dims[0] = P0;
    dims[1] = P1;

    MPI_Cart_create(MPI_COMM_WORLD, DIMS_NUM, dims, periods, reorder, &grid_comm); // создание коммуникатора
    n[0] = N1;
    n[1] = N2;
    n[2] = N3;
    if (current_process == 0) {

        A = (double *) malloc(N1 * N2 * sizeof(double));
        B = (double *) malloc(N2 * N3 * sizeof(double));
        C = (double *) malloc(N1 * N3 * sizeof(double));

        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < N2; ++j) {
                A[N2 * i + j] = 1;
            }
        }

        // Сразу транспонируем матрицу для ее дальнейшего перемножения
        for (int i = 0; i < N3; ++i) {
            for (int j = 0; j < N2; ++j) {
                B[N2 * i + j] = 2;
            }
        }

        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < N3; ++j) {
                C[N3 * i + j] = 0;
            }
        }


    }

    nn[0] = n[0] / dims[0]; // находим размер подматрицы AA для p[0] процессов
    nn[1] = n[2] / dims[1]; // находим размер подматрицы BB для p[1] процессов

    AA = (double*)malloc(nn[0] * n[1] * sizeof(double));
    BB = (double*)malloc(n[1] * nn[1] * sizeof(double));
    CC = (double*)malloc(nn[0] * nn[1] * sizeof(double));

    double start_time = MPI_Wtime();

    MPI_Comm_rank(grid_comm, &current_grid_rank); //Функция определения номера процесса в области grid_comm

    /* Функция определения координат процесса по его идентификатору
     * MPI_Cart_coords(MPI_Comm comm, int rank, int ndims, int *coords)
     * IN comm - коммуникатор с декартовой топологией;
     * IN rank - идентификатор процесса;
     * IN ndim - число измерений;
     * OUT coords - координаты процесса в декартовой топологии.
     */
    MPI_Cart_coords(grid_comm, current_grid_rank, DIMS_NUM, coords);

    MPI_Comm row_comm, col_comm;

    /* MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
     *IN comm - родительский коммуникатор;
     * IN color - признак подгруппы;
     * IN key - управление упорядочиванием;
     * OUT newcomm - новый коммуникатор.
     *
     * Функция расщепляет группу, связанную с родительским коммуникатором, на непересекающиеся подгруппы по одной
     * на каждое значение признака подгруппы color. Значение color должно быть неотрицательным.
     * Каждая подгруппа содержит процессы с одним и тем же значением color.
     * Параметр key управляет упорядочиванием внутри новых групп: меньшему значению key соответствует меньшее
     * значение идентификатора процесса. В случае равенства параметра key для нескольких процессов упорядочивание
     * выполняется в соответствии с порядком в родительской группе.
     */
    MPI_Comm_split(grid_comm, coords[0], coords[1], &row_comm);
    MPI_Comm_split(grid_comm, coords[1], coords[0], &col_comm);

    //Рассылаем подматрицы AA из процессов с coords[1] = 0 по процессам того же столбцового коммуникатора
    if (coords[1] == 0) {
        /*
         * Функция MPI_Scatter разбивает сообщение из буфера посылки процесса root на равные части размером sendcount
         * и посылает i-ю часть в буфер приема процесса с номером i (в том числе и самому себе).
         * Процесс root использует оба буфера (посылки и приема), поэтому в вызываемой им подпрограмме все параметры
         * являются существенными. Остальные процессы группы с коммуникатором comm являются только получателями,
         * поэтому для них параметры, специфицирующие буфер посылки, не существенны.
         *
         * int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
         *
         * IN sendbuf - адрес начала размещения блоков распределяемых данных (используется только в процессе-отправителе root);
         * IN sendcount - число элементов, посылаемых каждому процессу;
         * IN sendtype - тип посылаемых элементов; OUT recvbuf - адрес начала буфера приема;
         * IN recvcount - число получаемых элементов;
         * IN recvtype - тип получаемых элементов;
         * IN root - номер процесса-отправителя;
         * IN comm - коммуникатор.
         */
        MPI_Scatter(A, nn[0]*n[1], MPI_DOUBLE, AA, nn[0]*n[1], MPI_DOUBLE, 0, col_comm);
    }

    /*
     * Широковещательная рассылка данных выполняется с помощью функции MPI_Bcast. Процесс с номером root рассылает сообщение из своего буфера передачи всем процессам области связи коммуникатора comm.
     *
     * int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
     *
     * INOUT buffer - адрес начала расположения в памяти рассылаемых данных;
     * IN count - число посылаемых элементов;
     * IN datatype - тип посылаемых элементов;
     * IN root - номер процесса-отправителя;
     * IN comm - коммуникатор.
     */

    //Принимаем подматрицы AA из процесса с coords[1] = 0 столбцового коммуникатора процесса с coords[1]=0
    MPI_Bcast(AA, nn[0]*n[1], MPI_DOUBLE, 0, row_comm);


    //Рассылаем подматрицы BB из процессов с coords[0] = 0 по процессам того же строкового коммуникатора
    if (coords[0] == 0) {
        /*
         * Функция MPI_Scatter разбивает сообщение из буфера посылки процесса root на равные части размером sendcount
         * и посылает i-ю часть в буфер приема процесса с номером i (в том числе и самому себе).
         * Процесс root использует оба буфера (посылки и приема), поэтому в вызываемой им подпрограмме все параметры
         * являются существенными. Остальные процессы группы с коммуникатором comm являются только получателями,
         * поэтому для них параметры, специфицирующие буфер посылки, не существенны.
         *
         * int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
         *
         * IN sendbuf - адрес начала размещения блоков распределяемых данных (используется только в процессе-отправителе root);
         * IN sendcount - число элементов, посылаемых каждому процессу;
         * IN sendtype - тип посылаемых элементов; OUT recvbuf - адрес начала буфера приема;
         * IN recvcount - число получаемых элементов;
         * IN recvtype - тип получаемых элементов;
         * IN root - номер процесса-отправителя;
         * IN comm - коммуникатор.
         */
        MPI_Scatter(B, nn[1]*n[1], MPI_DOUBLE, BB, nn[1]*n[1], MPI_DOUBLE, 0, row_comm);
    }
    /*
     * Широковещательная рассылка данных выполняется с помощью функции MPI_Bcast. Процесс с номером root рассылает сообщение из своего буфера передачи всем процессам области связи коммуникатора comm.
     *
     * int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
     *
     * INOUT buffer - адрес начала расположения в памяти рассылаемых данных;
     * IN count - число посылаемых элементов;
     * IN datatype - тип посылаемых элементов;
     * IN root - номер процесса-отправителя;
     * IN comm - коммуникатор.
     */

    //Принимаем подматрицы BB из процесса с coords[0] = 0 строкового коммуникатора процесса с coords[0]=0
    MPI_Bcast(BB, nn[1]*n[1], MPI_DOUBLE, 0, col_comm);

    //Вычисляем подматрицу CC перемножением подматриц AA и BB
    for(int i = 0; i < nn[0]; i++){
        for(int j = 0; j < nn[1]; j++){
            CC[nn[1]*i + j] = 0.0;
            for(int k = 0; k < n[1]; k++){
		int pos1 = n[1]*i+k;
		int pos2 = n[1]*j+k;
                CC[nn[1]*i+j] += AA[pos1] * BB[pos2];
            }
        }
    }

    if (current_process!=0) {
        /*
         * Функция передачи сообщения MPI_Send
         *
         * int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
         *
         * IN buf - адрес начала расположения пересылаемых данных;
         * IN count - число пересылаемых элементов;
         * IN datatype - тип посылаемых элементов;
         * IN dest - номер процесса-получателя в группе, связанной с коммуникатором comm;
         * IN tag - идентификатор сообщения (аналог типа сообщения функций nread и nwrite PSE nCUBE2);
         * IN comm - коммуникатор области связи.
         */
	//
	MPI_Send(CC, nn[0] * nn[1], MPI_DOUBLE, 0, 0, grid_comm);
    }

    //Сборка матрицы С в 0 процессе
    if (current_process == 0) {
        //Идем вдоль оси строк
        for (int i = 0; i < dims[0]; i++) {
            //Идем вдоль столбцов
            for (int j = 0; j < dims[1]; j++) {
                if (i != 0 || j != 0) {
                    int tmp[2];
                    tmp[0] = i;
                    tmp[1] = j;
                    int rank_tmp;
                    MPI_Cart_rank(grid_comm, tmp, &rank_tmp);
                    MPI_Recv(CC, nn[0] * nn[1], MPI_DOUBLE, rank_tmp, 0, grid_comm, MPI_STATUS_IGNORE);
                }
                int base_i1 = nn[0]*i;
                int base_j1 = nn[1]*j;
                for (int i1 = 0; i1 < nn[0]; i1++) {
                    for (int j1 = 0; j1 < nn[1]; j1++) {
                        C[(i1+base_i1)*N3+(j1+base_j1)] = CC[i1*nn[1]+j1];
                    }
                }
            }
        }
    }

    double finish_time = MPI_Wtime();

    free(AA);
    free(BB);
    free(CC);

    if (current_process == 0) {
       printf("Processes: %d\n", number_of_processes);
       printf("Grid: %dx%d\n", P0, P1);
       printf("Time: %lf\n", finish_time - start_time); 
       
       /*for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N3; j++) {
                printf("%lf ", C[N3*i + j]);
            }
            printf("\n");
        }*/
       printf("Test val: %lf\n", C[52*N3+23]);

        free(A);
        free(B);
        free(C);
    }
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&grid_comm);
    MPI_Finalize();

    return 0;
}
