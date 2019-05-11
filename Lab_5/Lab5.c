#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <mpi.h>

int createThreads();
void* doTasks(void*);
void* sendTasks(void*);

//Кол-во задач
const int TASKS = 20;

//Кол-во задач у текущего процесса 
int processTaskNum;

//Флаг балансировки
int balance;

//Указатель на массив задач
int *tasks;

//Позиция "головы" очереди
int nextTask;

//Два объекта типа "описатель потока"
pthread_t thrs[2];

//Мьютекс
/*
  Для организации взаимного исключения потоков при доступе к разделяемым данным используются мьютексы (mutex = mutual exclusion), 
  объекты типа pthread_mutex_t. Мьютекс должен быть инициализирован перед использованием.

  Он либо заблокирован (занят, залочен, захвачен) каким-то потоком, либо свободен. 
  Поток, который захватил мьютекс, работает с участком кода. Остальные потоки, когда достигают мьютекса, ждут его разблокировки. 
  Разблокировать мьютекс может только тот поток, который его захватил. 
  Обычно освобождение занятого мьютекса происходит после исполнения критичного к совместному доступу участка кода.
*/
pthread_mutex_t mutex;

//Кол-во MPI-процессов
int size;

//Номер текущего MPI-процесса
int rank;

int main(int argc, char* argv[]) {
  //Фактически предоставленный уровень "поддержки потоков"
  int provided;
  int i;
  double t1, t2, mt;

  //Получаем флаг балансировки
  balance = atoi(argv[1]);

  //Инициилизируем многопоточный MPI
  /*
    int MPI_Init_thread(int *argc, char *((*argv)[]), int required, int *provided)

    argc      IN  То же самое, что и в MPI_Init
    argv      IN  То же самое, что и в MPI_Init

    required  IN  
    Запрашиваемый уровень "поддержки потоков":

    MPI_THREAD_SINGLE -- программист обещает, что программа будет иметь лишь один поток.
    MPI_THREAD_FUNNELED -- программист обещает, что только тот поток, который инициализировал MPI будет в дальнейшем вызывать функции MPI.
    MPI_THREAD_SERIALIZED -- программист обещает, что потоки не будут вызывать функции MPI одновременно.
    MPI_THREAD_MULTIPLE -- программист ничего не обещает. Любой поток может вызывать MPI-функции независимо от других потоков. Если реализация MPI поддерживает этот уровень, то это и есть thread-compliant реализация MPI.

    provided  OUT Фактически предоставленный уровень  


    MPI_Init озанчает MPI_Init_thread с уровнем MPI_THREAD_SINGLE.    
  */
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(provided != MPI_THREAD_SERIALIZED) {
    printf("Can't get needed level!\n");
    MPI_Finalize();
    return 1;
  }

  //Инициилизируем mutex перед использованием
  /*
    int pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t *attr);

    Параметр attr (IN) задает атрибуты мьютекса. Можно передать NULL для принятия атрибутов по умолчанию.
  */
  pthread_mutex_init(&mutex, NULL);

  //Вычисляем кол-во задач для MPI-процессса
  processTaskNum = TASKS / size;

  //Выделяем память для массива задач 
  tasks = (int*)malloc(sizeof(int) * processTaskNum);

  //Заполняем массив задач случайными задачами 
  for(i = 0; i < processTaskNum; i++){
    tasks[i] = TASKS / size * rank + rand() % (TASKS / size);
  }

  //Выполняем измерение времени
  t1 = MPI_Wtime();

  //Функция, предназначенная для создания потоков
  createThreads();

  //Выполняем измерение времени
  t2 = MPI_Wtime();

  //Определяем разницу между показаниями таймера в начале и в конце выполнения функции 
  t2 -= t1;

  //Находим максимальное время работы функции createThreads() - это время и будет итоговым временем работы программы
  /*
    int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)


    IN  sendbuf   - адрес начала входного буфера;
    OUT recvbuf   - адрес начала буфера результатов (используется только в процессе-получателе root);
    IN  count     - число элементов во входном буфере;
    IN  datatype  - тип элементов во входном буфере;
    IN  op        - операция, по которой выполняется редукция;
    IN  root      - номер процесса-получателя результата операции;
    IN  comm      - коммуникатор
  */
  MPI_Reduce(&t2, &mt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  //Выводим на экран информацию о кол-ве задач, потоков и времени работы програмы
  if(rank == 0){
    printf("Tasks = %d\nThreads = %d\nTime = %f\n", TASKS, size, mt);
  }

  MPI_Finalize();

  return 0;
}

int createThreads(){
  //Объект, задающий атрибуты потока имеет тип pthread_attr_t. Такой объект должен быть инициализирован с помощью функции pthread_attr_init
  pthread_attr_t attrs;
  int i;

  //Иницилизуем объект, задающий атрибуты потока
  /*
    int pthread_attr_init(pthread_attr_t *attr);

    В результате объект будет содержать набор свойств потока по умолчанию
  */
  if(pthread_attr_init(&attrs) != 0){
      perror("Cannot initialize attributes!");
      return 1;
  };

  //Устанавливаем в аттрибуте свойство "присоединяемости"(joinable) для потоков посредством вызова функции pthread_attr_setdetachstate
  /*
    int pthread_attr_setdetachstate(pthread_attr_t *attr, int detachstate);

    Поток может быть "присоединяемым" (joinable) или "оторванным" (detached), где
    etachstate можно установить в PTHREAD_CREATE_JOINABLE или в PTHREAD_CREATE_DETACHED соответственно.
  */
  if(pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE) != 0){
    perror("Error in setting attributes!");
    return 1;
  }

  //Создаем DO-поток
  /*
    int pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*start_routine)(void*), void *arg);
    
    thread        OUT     В результате успешного срабатывания функции по указанному адресу будет размещен описатель порожденного потока.
    attr          IN      Атрибуты потока. Задают свойства потока. Может быть NULL.
    start_routine IN      Указатель на функцию потока. Выполнение потока состоит в выполнении этой функции.
    arg           IN/OUT  Указатель, который будет передан функции потока в качестве параметра. 
                          OUT -- в том смысле, что функция потока может менять содержимое памяти с использованием этого указателя. 
                          pthread_create сожержимого не меняет, просто передает указатель функции потока. 
                          Функция потока сама интерпретирует содержание памяти по этому адресу.
  */
  if(pthread_create(&thrs[0], &attrs, doTasks, NULL) != 0){
    perror("Cannot create a DO thread!");
    return 1;
  }

  //Если установлен флаг балансировки, то создаем SEND-поток
  if(balance == 1){
    if(pthread_create(&thrs[1], &attrs, sendTasks, NULL) != 0){
      perror("Cannot create a SEND thread!");
      return 1;
    }
  }

  //Освобождаем ресурсы атрибута
  /*
    Ресурсы, которые могут использоваться в системе для хранения атрибутов освобождаются 
    вызовом функции (после того, как объект был использован в вызове pthread_create и больше не нужен)

    int pthread_attr_destroy(pthread_attr_t *attr);
  */
  pthread_attr_destroy(&attrs);

  /*
    int pthread_join(pthread_t thread, void **value_ptr);

    value_ptr (OUT) - это указатель на указатель, возвращенный функцией завершившегося потока.

    Поток, вызвавший эту функцию, останавливается, пока не окончится выполнение потока thread. 
    Если никто не вызывает pthread_join для присоединяемого потока, то завершившись, он не освобождает свои ресурсы, 
    а это может служить причиной утечки памяти в программе
  */
  if(pthread_join(thrs[0], NULL) != 0){
    perror("Cannot join a thread");
    return 1;
  }

  if(balance == 1){
    if(pthread_join(thrs[1], NULL) != 0){
      perror("Cannot join a thread");
      return 1;
    }
  }

  return 0;
}

void* doTasks(void* args) {
  int time;
  int i;
  int request;

  for(nextTask = 0; nextTask < processTaskNum;) {
    /*Блокировка: теперь к ресурсам имеет доступ только один поток, который владеет мьютексом. Он же единственный, 
    кто может его разблокировать. Поток блокируется, если он пытается заблокировать мьютекс, который уже заблокирован. 
    Он пробуждается, когда мьютекс разблокируется. 
    */
    pthread_mutex_lock(&mutex);
    //Получаем задачу
    time = tasks[nextTask];
    tasks[nextTask] = -1;
    nextTask++;
    pthread_mutex_unlock(&mutex);
    //Выполняем задачу
    sleep(time);
    printf("Rank#%d did task.\n", rank);
  }

  //Принимаем задачи 
  if (balance == 1) {
    //Проходимся по всем процессам  
    for (i = 0; i < size; i++) {
      request = 1;

      while (1) {
        /*
          Функция передачи сообщения MPI_Send
         
          int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
         
          IN buf - адрес начала расположения пересылаемых данных;
          IN count - число пересылаемых элементов;
          IN datatype - тип посылаемых элементов;
          IN dest - номер процесса-получателя в группе, связанной с коммуникатором comm;
          IN tag - идентификатор сообщения (аналог типа сообщения функций nread и nwrite PSE nCUBE2);
          IN comm - коммуникатор области связи.
        */
        MPI_Send(&request, 1, MPI_INT, (rank + i) % size, 0, MPI_COMM_WORLD);
        int answer;
        //Принимаем задачу
        /*
          int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)

          OUT buf - адрес начала расположения принимаемого сообщения;
          IN  count - максимальное число принимаемых элементов;
          IN  datatype  - тип элементов принимаемого сообщения;
          IN  source  - номер процесса-отправителя;
          IN  tag - идентификатор сообщения;
          IN  comm  - коммуникатор области связи;
          OUT status  - атрибуты принятого сообщения.
        */
        MPI_Recv(&answer, 1, MPI_INT, (rank + i) % size, 1, MPI_COMM_WORLD, 0);
        //Если задач на исполнение нет, то завершаем работу цикла
        if (answer == -1) break;
        //Выполняем задачу
        sleep(answer);

        printf("Rank#%d did alien task from Rank#%d.\n", rank, (rank + i) % size);
      }
    }

    //Синхронизируем все процессы
    MPI_Barrier(MPI_COMM_WORLD);

    //Таким образом информируем поток-отправки об окончании работы 
    request = 0;
    MPI_Send(&request, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
  }
  printf("Rank#%d finished tasks.\n", rank);
}

void* sendTasks(void* args) {
  MPI_Status st;
  int request;
  int answer;

  while(1) {
    //Ждем информации от других процессов о том, что они готовы принимать задачи на исполнение 
    MPI_Recv(&request, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &st);
    if(request == 0) break;
    //Блокируем очередь и забираем оттуда задачу
    pthread_mutex_lock(&mutex);

    //Если задач на отсылку больше нет, 
    if(nextTask == processTaskNum) {
      answer = -1;
    } else {
      answer = tasks[nextTask];
      nextTask++;
    }

    pthread_mutex_unlock(&mutex);

    MPI_Send(&answer, 1, MPI_INT, st.MPI_SOURCE, 1, MPI_COMM_WORLD);
  }
  printf("Rank#%d finished sending tasks.\n", rank);
}