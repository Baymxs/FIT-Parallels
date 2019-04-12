#include <cstring>
#include <algorithm>
#include <iostream>
#include <mpi.h>

using namespace std;

const double eps = 1e-8;
const double a = 1e5;
const int N = 100; //Кол-во точек в каждом измерении
const double min_x = -1, min_y = -1, min_z = -1;
const double max_x = 1,  max_y = 1,  max_z = 1;
const bool print_matrix = true;

//Считаем значение функции фи по переданным координатам 
double phi(double x, double y, double z) {
	return x*x+y*y+z*z;	
}

double ro(double x, double y, double z) {
	return 6-a*phi(x,y,z);
}

double recalc_layer(int base_z, int z, double* matrix, double* tmp_matrix, double delta_x, double delta_y, double delta_z) {
	//Абсолютный номер слоя
	int abs_z = base_z+z;
	//Если это граничный слой
	if (abs_z==0 || abs_z==N-1) {
		//Копируем этот слой в новый массив на старое место, не пересчитывая
		memcpy(tmp_matrix + z*N*N, matrix + z*N*N, N * N * sizeof(double));
		return 0;
	}
	//Иначе пересчитываем каждый элемент слоя с помощью итерационной формулы Якоби
	double max_delta = 0;
	double cur_z = min_z+abs_z*delta_z;
	for (int i=0;i<N;i++) {
		double cur_x = min_x+i*delta_x;
		for (int j=0;j<N;j++) {
			double cur_y = min_y+j*delta_y;
			//Если элемент находится на границе слоя, то не пересчитываем его 
			if (i==0 || i==N-1 || j==0 || j==N-1) {
				tmp_matrix[z*N*N+i*N+j] = matrix[z*N*N+i*N+j];
				continue;
			}
			int cell = z*N*N+i*N+j;
			double tmp = (matrix[z*N*N+(i+1)*N+j]+matrix[z*N*N+(i-1)*N+j])/(delta_x*delta_x);
			tmp += (matrix[z*N*N+i*N+(j+1)]+matrix[z*N*N+i*N+(j-1)])/(delta_y*delta_y);
			tmp += (matrix[(z+1)*N*N+i*N+j]+matrix[(z-1)*N*N+i*N+j])/(delta_z*delta_z);
			tmp -= ro(cur_x, cur_y, cur_z);
			tmp_matrix[cell] = 1/(2/(delta_x*delta_x)+2/(delta_y*delta_y)+2/(delta_z*delta_z)+a);
			tmp_matrix[cell]*=tmp;
			max_delta = max(max_delta, abs(tmp_matrix[cell]-matrix[cell]));
		}
	}
	return max_delta;
}

int main(int argc, char** argv) {
	int sz, id; 

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	double st_time = MPI_Wtime();

	double delta_x = (max_x - min_x) / (N - 1); //Расстояние между двумя соседними точками в измерении Ox
	double delta_y = (max_y - min_y) / (N - 1); //Расстояние между двумя соседними точками в измерении Oy
	double delta_z = (max_z - min_z) / (N - 1); //Расстояние между двумя соседними точками в измерении Oz

	int part_z = N / sz; //Делим область решения на подобласти

	double* matrix = new double[(part_z + 2) * N * N]; //Трехмерный массив для хранения значений функции в точках слоя, +2 - добавляем образ верхнего и нижнего внешних граничных слоев

	double* tmp_matrix = new double[(part_z + 2) * N * N];

	int base_z = id * part_z-1; //Номер по Oz, с которой начинается первый элемент слоя для данного процесса 

	//Инициализация слоя 
	for (int i = 0; i < part_z+2; i++) {
		int real_i = i + base_z; //Позиция текущего эелемента по Oz
		double cur_z = min_z + delta_z*real_i; //Координата по Oz, от которой будет браться функция
		for (int j = 0; j < N; j++) {
			double cur_x = min_x + delta_x*j; //Координата по Ox, от которой будет браться функция
			for (int k = 0; k < N; k++) {
				double cur_y = min_y + delta_y*k; //Координата по Oy, от которой будет браться функция
				int cell = i*N*N + j*N + k; //По координате получаем индекс в одномерном массиве
				//Вычисляем значение на границе
				if (real_i == 0 || real_i == N-1 || j == 0 || j == N-1 || k == 0 || k == N-1) {
					//Считаем значение функции фи
					matrix[cell] = phi(cur_x, cur_y, cur_z);
				} else {
					matrix[cell] = 0;
				}
			}
		}
	}

	//Если размер матрицы не делится нацело на количество процессов, то завершаем работу программы
	if (N % sz) {
		if (id == 0) {
			cout << "Invalid number of processes" << endl;
		}
		MPI_Finalize();
		return 0;
	}
	
	//Бесконечный цикл
	for (int iter = 0;;iter++) {
		double max_delta = 0;
		//Пересчитываем граничный верхний слой(плоскость)
		double tmp_delta = recalc_layer(base_z, 1, matrix, tmp_matrix, delta_x, delta_y, delta_z);
		max_delta = max(max_delta, tmp_delta);
		//Пересчитываем граничный нижний слой(плоскость)
		tmp_delta = recalc_layer(base_z, part_z, matrix, tmp_matrix, delta_x, delta_y, delta_z);
		max_delta = max(max_delta, tmp_delta);
		//Массив для хранения информации о запрсах
		MPI_Request rq[4];
		//Если это ненулевой процесс
		if (id != 0) {
			//Отправляем верхний граничный слой(плоскость) предыдущему процеесу
			/*
				Функция передачи сообщения без блокировки MPI_Isend

				C:
				int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest,
				int tag, MPI_Comm comm, MPI_Request *request)


				IN	buf	- адрес начала расположения передаваемых данных;
				IN	count	- число посылаемых элементов;
				IN	datatype	- тип посылаемых элементов;
				IN	dest	- номер процесса-получателя;
				IN	tag	- идентификатор сообщения;
				IN	comm	- коммуникатор;
				OUT	request	- "запрос обмена"
			*/
			MPI_Isend(tmp_matrix+N*N, N*N, MPI_DOUBLE, id-1, 123, MPI_COMM_WORLD, &rq[0]);		

			//Получаем верхний внешний граничный слой от предыдущего процесса
			/*
				Функция приема сообщения без блокировки MPI_Irecv

				C:
				int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source,
				int tag, MPI_Comm comm, MPI_Request *request)
		
				OUT	buf	- адрес для принимаемых данных;
				IN	count	- максимальное число принимаемых элементов;
				IN	datatype	- тип элементов принимаемого сообщения;
				IN	source	- номер процесса-отправителя;
				IN	tag	- идентификатор сообщения;
				IN	comm	- коммуникатор;
				OUT	request	- "запрос обмена"
			*/
			MPI_Irecv(tmp_matrix, N*N, MPI_DOUBLE, id-1, 123, MPI_COMM_WORLD, &rq[2]);
		}	
		//Если это не последний процесс
		if (id != sz-1) {
			//Отправляем нижний граничный слой(плоскость) следующему процеесу
			MPI_Isend(tmp_matrix+part_z*N*N, N*N, MPI_DOUBLE, id+1, 123, MPI_COMM_WORLD, &rq[1]);
			//Получаем нижний внешний граничный слой от следующего процесса
			MPI_Irecv(tmp_matrix+(part_z+1)*N*N, N*N, MPI_DOUBLE, id+1, 123, MPI_COMM_WORLD, &rq[3]);
		}
		//Если это ненулевой процесс
		if (id != 0) {
			/*
				C:
				int MPI_Wait(MPI_Request *request, MPI_Status *status)
				
				request	- "запрос обмена";
				OUT	status	- атрибуты сообщения.

				Это нелокальная блокирующая операция. Возврат происходит после завершения операции, связанной с запросом request. 
				В параметре status возвращается информация о законченной операции.
			*/
			MPI_Wait(&rq[0], MPI_STATUS_IGNORE);
			MPI_Wait(&rq[2], MPI_STATUS_IGNORE);
		}
		//Если это не последний процесс
		if (id != sz-1) {
			MPI_Wait(&rq[1], MPI_STATUS_IGNORE);
			MPI_Wait(&rq[3], MPI_STATUS_IGNORE);
		}
		//Пересчитываем все элементы не граничных и не граничных внешних слоев
		for (int i = 2; i < part_z; i++) {
			double tmp_delta = recalc_layer(base_z, i, matrix, tmp_matrix, delta_x, delta_y, delta_z);
			max_delta = max(max_delta, tmp_delta);
		}
		//"Архивируем" полностью пересчитанную область для данного процесса  	
		memcpy(matrix, tmp_matrix, (part_z+2)*N*N*sizeof(double));
		double max_delta_shared;
		//Находим максимальную дельта из все процессов
		/*
			Функция MPI_Reduce выполняется следующим образом. Операция глобальной редукции, указанная параметром op, 
			выполняется над первыми элементами входного буфера, и результат посылается в первый элемент буфера приема процесса root. 
			Затем то же самое делается для вторых элементов буфера и т.д.

			С:
			int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
			MPI_Op op, int root, MPI_Comm comm)

			IN	sendbuf	-	адрес начала входного буфера;
			OUT	recvbuf	-	адрес начала буфера результатов (используется только в процессе-получателе root);
			IN	count	-	число элементов во входном буфере;
			IN	datatype	-	тип элементов во входном буфере;
			IN	op	-	операция, по которой выполняется редукция;
			IN	root	-	номер процесса-получателя результата операции;
			IN	comm	-	коммуникатор.
		*/
		MPI_Reduce(&max_delta, &max_delta_shared, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		//Отправляем найденную ранее максимальную дельта всем остальным процессам
		/*
			Широковещательная рассылка данных выполняется с помощью функции MPI_Bcast. Процесс с номером root рассылает сообщение из своего буфера передачи всем процессам области связи коммуникатора comm.

			С:
			int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root,
			MPI_Comm comm)
				
			INOUT	buffer	- адрес начала расположения в памяти рассылаемых данных;
			IN	count	- число посылаемых элементов;
			IN	datatype	- тип посылаемых элементов;
			IN	root	- номер процесса-отправителя;
			IN	comm	- коммуникатор.

			После завершения подпрограммы каждый процесс в области связи коммуникатора comm, включая и самого отправителя, получит копию сообщения от процесса-отправителя root. 
		*/
		MPI_Bcast(&max_delta_shared, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//Условия выхода из бесконечного цикла
		if (max_delta_shared < eps) break;
	}

	double* fullRes;
	if (id == 0) {
		fullRes = new double[N*N*N];
	}

	//Выполняем сборку посчитанных областей 
	MPI_Gather(matrix+N*N, part_z*N*N, MPI_DOUBLE, fullRes, part_z*N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//В нулевом процессе считаем максимальное отклонение от эталонного ответа
	if (id == 0) {
		double max_delta = 0;
		for (int layer = 0;layer < N;layer++){
			double z = min_z + layer*delta_z;
			for (int j = 0;j < N;j++) {
				double x= min_x + j*delta_x;
				for (int k = 0; k < N; k++) {
					double y = min_y + k*delta_y;
					double tmp = fullRes[layer*N*N + j*N + k];
					double val = phi(x, y, z);
					max_delta = max(max_delta, abs(tmp-val));
				}
			}
		}
		cout << "Answer delta: " << max_delta << endl;
		delete[] fullRes;
	}

	double en_time = MPI_Wtime();

	double delta_time = en_time - st_time;

	if (id == 0) {
		printf("Time: %lf\n", delta_time);
	}

	MPI_Finalize();
	return 0;
}
