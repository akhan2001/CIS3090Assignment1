/**
 * @file gol_data.c
 * @author Abdullah Khan (mkhan94@uoguelph.ca) - 1101209
 * @brief This file does the Data Parallelism algorithm for the Conway's Game of Life
 * Uses two arrays to act as the grid of the game. Where one array is read-from and 
 * the other is written-to.
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

/*Global variables*/
int **readFrom;
int **writeTo;

int nThreads;
int gridSize;
int nIterations;

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
 * @brief Calculates the next generation grid based on Conway's game of life
 * @param arg 
 * @return ((void *)0)
 */
void *iterateNextGrid(void *arg) {
    long index = (long)arg;
    /*Param *param = (Param *)arg;*/
    int threadNum = gridSize / nThreads;
    int startOffset = index * threadNum;
    int endOffset = (index + 1) * threadNum - 1;

    for (int i = startOffset; i <= endOffset; i++) {
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

            if (readFrom[i][j] == 1) {
                if (neighbors == 2 || neighbors == 3) {
                    writeTo[i][j] = 1;
                } else {
                    writeTo[i][j] = 0;
                }
            } else {
                if (neighbors == 3) {
                    writeTo[i][j] = 1;
                } else {
                    writeTo[i][j] = 0;
                }
            }
        }
    }

    return((void *)0);
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

/**
 * @brief Main function
 * @param argc 
 * @param argv 
 * @return int num
 */
int main(int argc, char const *argv[]) {
    int             display = 0;
    int             err;
    double          elapsed;
	struct timeval  start;
    pthread_t       tid[nThreads];
    
    /*if (argc < 4) {
        printf("~ Enter \"./gol_data [# of threads] [grid-size] [# of iterations] [-d]\" ~\n");
        exit(0);
    }
    if (argc == 5) {
		if (strcmp(argv[4], "-d") == 0){
			display = 1;
		}
	}*/

    /*Basic error-checking/initializing for command-line arguments*/
    if(argc == 5 && strcmp(argv[4], "-d") == 0) {
        display = 1;
    } else if(argc == 5 && strcmp(argv[4], "-d") != 0) {
        printf("~ Enter \"./gol_data [# of threads] [grid-size] [# of iterations] [-d]\" ~\n");
        exit(0);
    } else if(argc != 4) {
        printf("~ Enter \"./gol_data [# of threads] [grid-size] [# of iterations] [-d]\" ~\n");
        exit(0);
    }

    nThreads = atoi(argv[1]);
	gridSize = atoi(argv[2]);
	nIterations = atoi(argv[3]);

    if (nThreads < 1) {
        printf("~ Threads has to be greater than zero (0) ~\n");
        exit(0);
    }
    if (nThreads > gridSize) {
        nThreads = gridSize;
    }

    /*Initialize the grid-size of the board*/
    allocateArrays(gridSize);
    initializeArrays(gridSize);

    /*If -d is in the arguments parameter, then display the grid*/
    if (display == 1) {
        displayGrid(readFrom, gridSize);
    }

    /*Get the timing data*/
	gettimeofday(&start, NULL);

    for (int i = 1; i < nIterations; i++) {
        /*Create number of threads for threads to calculate the work*/
        for (long j = 0; j < nThreads; j++) {
            //long startOff = j * threadNum;
            err = pthread_create(&tid[j], NULL, iterateNextGrid, (void*)j);
            if (err != 0) {
                fprintf(stderr, "~ Cannot create thread, error: %d ~", err);
                exit(0);
            }
        }
        /*Wait for worker threads to finish and then copy the grid of arrays*/
        for (long j = 0; j < nThreads; j++) {
            pthread_join(tid[j], NULL);
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

    return 0;
}