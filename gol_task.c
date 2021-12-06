/**
 * @file gol_task.c
 * @author Abdullah Khan (mkhan94@uoguelph.ca) - 1101209
 * @brief 
 * @version 0.1
 * @date 2021-10-04
 * 
 * @copyright Copyright (c) 2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

/*Number of threads*/
#define NUMTHREADS 3

/*Point struct is needed in order to tell which queue is to be located next*/
typedef struct Point {
	int x;
	int y;
	struct Point * next;
} Point;

/*Queue struct needed for liveQueue and emptyQueue*/
typedef struct {
	Point * head;
	Point * tail;
    pthread_mutex_t mutex;
    //Add a condition variable
    pthread_cond_t cond;
} Queue;

/*Global variables*/
int **readFrom;
int **writeTo;

Queue *liveQueue;
Queue *emptyQueue;

int gridSize;
int nIterations;
int finished;

//Mutex used to protect from race conditions
pthread_mutex_t g_Mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t dataNotProduced = PTHREAD_COND_INITIALIZER;
pthread_cond_t dataNotConsumed = PTHREAD_COND_INITIALIZER;

/**
 * @brief Calculates the elapsed time from when it was called to when it ends
 * (Got function code from parallelSort.c provided by Prof. Denis)
 * @param start 
 * @return double time
 */
double calcTime(struct timeval start) {
    long long startusec, endusec;
    struct timeval	end;
    
    gettimeofday(&end, NULL);
    startusec = start.tv_sec * 1000000 + start.tv_usec;
    endusec = end.tv_sec * 1000000 + end.tv_usec;
    return (double)(endusec - startusec) / 1000000.0;
}

/**
 * @brief Displays the grid with the given current grid
 * @param currentGrid 
 * @param gridSize 
 */
void displayGrid(int **currentGrid, int gridSize) {
    /*Using unbuffered print statements*/
    for (int i = 0; i < gridSize; i++) {
       fprintf(stderr, "—+");
    }
    fprintf(stderr, "\n");
    
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            if (currentGrid[i][j] == 1) {
                fprintf(stderr, "#");
            } else {
                fprintf(stderr, " ");
            }
            fprintf(stderr, "|");
        }
        fprintf(stderr, "\n");
        for (int k = 0; k < gridSize; k++) {
            fprintf(stderr, "—+");
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n\n");
}

/**
 * @brief Sets the first generation grid with certain number of random live cells 
 * @param gridSize 
 */
void initializeArrays(int gridSize) {
    /*Random num generator*/
    srand(time(0));
    int randCells;

    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            /*Goes for potentionally only 40% of the grid*/
            randCells = (rand() % (5 - 1 + 1)) + 1;
            if (randCells == 1 || randCells == 2) {
                /*Alive cells*/
                readFrom[i][j] = 1;
            } else {
                /*Dead cells*/
                readFrom[i][j] = 0;
            }
        }
    }
}

/**
 * @brief Allocates the two 2D arrays with proper memory size given grid-size
 * @param gridSize 
 */
void allocateArrays(int gridSize) {
    /*Mallocs the readFrom array*/
	readFrom = (int**)malloc(gridSize * sizeof(int*));
    for (int i = 0; i < gridSize; i++) {
        readFrom[i] = (int*)malloc(gridSize * sizeof(int));
    }
    /*Mallocs the writeTo array*/
	writeTo = (int**)malloc(gridSize * sizeof(int*));
    for (int i = 0; i < gridSize; i++) {
        writeTo[i] = (int*)malloc(gridSize * sizeof(int));
    }
}

/**
 * @brief Copys the readFrom array grid to the next generation writeTo array grid
 * @param gridSize 
 */
void copyGrid(int gridSize) {
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            readFrom[i][j] = writeTo[i][j];
        }
    }
}

/******************************************
 * New functions needed for gol_task.c
 *****************************************/
/**
 * @brief Create a Queue object and initialize its mutex and condition variable
 * @return Queue* q
 */
Queue *createQueue() {
    Queue* q = (Queue*)malloc(sizeof(Queue));
    q->head = q->tail = NULL;
    pthread_mutex_init(&q->mutex, NULL);
    //Initialize the condition variable
    pthread_cond_init(&q->cond, NULL);
    return q;
}

/**
 * @brief Creates a Point object and initialize the variables inside
 * @return Point* p
 */
Point *createPoint(int x, int y) {
    Point *p = (Point *)malloc(sizeof(Point));
    p->next = NULL;
    p->x = x;
    p->y = y;
    return p;
}

/**
 * @brief Check if given queue is empty
 * @return int 
 */
int isEmpty(Queue *queue) {
    if (queue->head == NULL) {
        return 1;
    } else {
        return 0;
    }
}

/**
 * @brief Calculates the next generation grid based on Conway's game of life
 * @param arg 
 * @return ((void *)0)
 */
void *iterateNextGrid() {
    //long index = (long)arg;
    /*Param *param = (Param *)arg;*/
    //int threadNum = gridSize / nThreads;
    /*int startOffset = index * threadNum;
    int endOffset = (index + 1) * threadNum - 1;*/

    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int neighbors = 0;
            if (i > 0) {
                if (readFrom[i - 1][j] == 1) {
                    neighbors++;
                }
                if (j > 0) {
                    if (readFrom[i - 1][j - 1] == 1) {
                        neighbors++;
                    }
                }
                if (j < gridSize - 1) {
                    if (readFrom[i - 1][j + 1] == 1) {
                        neighbors++;
                    }
                }
            }
            if (i < (gridSize - 1)) {
                if (readFrom[i+1][j] == 1) {
					neighbors++;
                }
				if (j > 0) {
					if (readFrom[i+1][j-1] == 1) {
						neighbors++;
                    }
				}
				if (j < (gridSize-1)) {
					if (readFrom[i+1][j+1] == 1) {
						neighbors++;
                    }
				}
            }
            if (j > 0) {
                if (readFrom[i][j - 1] == 1) {
                    neighbors++;
                }
            }
            if (j < (gridSize - 1)) {
                if (readFrom[i][j + 1] == 1) {
                    neighbors++;
                }
            }
            /*Queue has a point like (10,10) which means in the next interation
            the location (10,10) should be occupied*/
            Point *point = createPoint(i, j);
            /*Lets know if queue is occupied or unoccupied*/
            if (readFrom[i][j] == 1) {
                if (neighbors == 2 || neighbors == 3) {
                    //This position is occupied
                    if (liveQueue->tail != NULL) {
                        liveQueue->tail->next = point;
                        liveQueue->tail = point;
                    } else {
                        liveQueue->head = point;
		                liveQueue->tail = point;
                    }
                    //writeTo[i][j] = 1;
                } else {
                    //This position is unoccupied
                    if (emptyQueue->tail != NULL) {
                        emptyQueue->tail->next = point;
                        emptyQueue->tail = point;
                    } else {
                        emptyQueue->head = point;
		                emptyQueue->tail = point;
                    }
                    //writeTo[i][j] = 0;
                }
            } else {
                if (neighbors == 3) {
                    //This position is occupied
                    if (liveQueue->tail != NULL) {
                        liveQueue->tail->next = point;
                        liveQueue->tail = point;
                    } else {
                        liveQueue->head = point;
		                liveQueue->tail = point;
                    }
                    //writeTo[i][j] = 1;
                } else {
                    //This position is unoccupied
                    if (emptyQueue->tail != NULL) {
                        emptyQueue->tail->next = point;
                        emptyQueue->tail = point;
                    } else {
                        emptyQueue->head = point;
		                emptyQueue->tail = point;
                    }
                    //writeTo[i][j] = 0;
                }
            }
        }
    }
    finished = 1;
    return((void *)0);
}

/**
 * @brief Reads through the live queue and updates the next iteration of the
 * game board. Marks the location of the board as occupied
 * @return void* NULL
 */
void *nextLiveQueue() {
    while (finished == 0 || isEmpty(liveQueue) == 0 || isEmpty(emptyQueue) == 0) {
        Point *p;
        if (liveQueue->head == NULL) {
            p = NULL;
        } else if (liveQueue->head->next == NULL) {
            Point * tempPoint = liveQueue->head;
            liveQueue->head = NULL;
            liveQueue->tail = NULL;
            p = tempPoint;
        } else {
            Point * tempPoint = liveQueue->head;
            liveQueue->head = liveQueue->head->next;
            p = tempPoint;
        }

		if (p != NULL) {
            //Marks as 1 for occupied
			writeTo[p->x][p->y] = 1;
            free(p);
		}
    }
    return NULL;
}

/**
 * @brief Reads through the empty queue and updates the next iteration
 * of the game board. Marks the location on the board as unoccupied
 * @return void* NULL
 */
void *nextEmptyQueue() {
    while (finished == 0 || isEmpty(liveQueue) == 0 || isEmpty(emptyQueue) == 0) {
        Point *p;
        if (emptyQueue->head == NULL) {
            p = NULL;
        } else if (emptyQueue->head->next == NULL) {
            Point * tempPoint = emptyQueue->head;
            emptyQueue->head = NULL;
            emptyQueue->tail = NULL;
            p = tempPoint;
        } else {
            Point * tempPoint = emptyQueue->head;
            emptyQueue->head = emptyQueue->head->next;
            p = tempPoint;
        }

		if (p != NULL) {
            //Marks as 0 for unoccupied
			writeTo[p->x][p->y] = 0;
            free(p);
		}
    }
    return NULL;
}

/**
 * @brief Frees the queue of any allocated memory
 * @param removedQueue 
 */
void freeQueue(Queue *removedQueue) {
    Point *p = removedQueue->head;
    while (p != NULL) {
        Point *tempPoint = p;
        p = p->next;
        free(tempPoint);
    }
    free(removedQueue);
}

int main(int argc, char const *argv[]) {
    int             display = 0;
    int             err;
    double          elapsed;
	struct timeval  start;
    pthread_t       tid[NUMTHREADS];

    /*Basic error-checking/initializing for command-line arguments*/
    if(argc == 4 && strcmp(argv[3], "-d") == 0) {
        display = 1;
    } else if(argc == 4 && strcmp(argv[3], "-d") != 0) {
        printf("~ Enter \"./gol_data [grid-size] [# of iterations] [-d]\" ~\n");
        exit(0);
    } else if(argc != 3) {
        printf("~ Enter \"./gol_data [grid-size] [# of iterations] [-d]\" ~\n");
        exit(0);
    }

    /*Initializing argument values to variables*/
	gridSize = atoi(argv[1]);
	nIterations = atoi(argv[2]);

    /*readFrom = malloc(gridSize * sizeof(*readFrom));
    writeTo = malloc(gridSize * sizeof(*writeTo));*/
    
    /*Initialize the grid-size of the board*/
    allocateArrays(gridSize);
    initializeArrays(gridSize);

    liveQueue = createQueue();
    emptyQueue = createQueue();

    /*If -d is in the arguments parameter, then display the grid*/
    if (display == 1) {
        displayGrid(readFrom, gridSize);
    }

    /*Get the timing data*/
	gettimeofday(&start, NULL);

    for (int i = 0; i < nIterations; i++) {
        finished = 0;
        err = pthread_create(&tid[0], NULL, iterateNextGrid, NULL);
        if (err != 0) {
            fprintf(stderr, "~ Cannot create thread, error: %d ~", err);
            exit(0);
        }
        err = pthread_create(&tid[1], NULL, nextLiveQueue, NULL);
        if (err != 0) {
            fprintf(stderr, "~ Cannot create thread, error: %d ~", err);
            exit(0);
        }
        err = pthread_create(&tid[2], NULL, nextEmptyQueue, NULL);
        if (err != 0) {
            fprintf(stderr, "~ Cannot create thread, error: %d ~", err);
            exit(0);
        }
        /*Wait for worker threads to finish and then copy the grid of arrays*/
        for (int i = 0; i < NUMTHREADS; i++) {
            pthread_join(tid[i], NULL);
        }
        /*Copys the readFrom array to the (next) writeTo array*/
        copyGrid(gridSize);
        /*If -d is in the arguments parameter, then display the grid*/
        if (display == 1) {
            displayGrid(readFrom, gridSize);   
        }
    }

    /*Get the elapsed time is took for the threads to do the work*/
    elapsed = calcTime(start);
	printf("This game took %.4f seconds\n\n", elapsed);

    /*Free the allocated memory*/
    for (int i = 0; i < gridSize; i++) {
        free(readFrom[i]);
        free(writeTo[i]);
    }
    free(readFrom);
	free(writeTo);

    freeQueue(liveQueue);
    freeQueue(emptyQueue);

    /*for (int i = 0; i < NUMTHREADS; i++) {
        args[i].workerId = i + 1;
        args[i].q = q;
        pthread_create(&tid[i], NULL, workerFunc, &args[i]);
    }*/
    /*Create worker threads
    for (i = 0; i < nThreads; i++) {
        args[i].workerId = i + 1;
        args[i].q = q;
        pthread_create(&tid[i], NULL, workerFunc, &args[i]);
    }*/
    /*for (int i = 0; i < gridSize; i++) {
        readGrid[i] = malloc(gridSize * sizeof(*readGrid[i]));
        writeGrid[i] = malloc(gridSize * sizeof(*writeGrid[i]));
    }*/
    return 0;
}