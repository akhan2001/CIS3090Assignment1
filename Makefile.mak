all: gol_data.o gol_task.o

gol_data.o: gol_data.c
	gcc gol_data.c -o gol_data -lpthread

gol_task.o: gol_task.c
	gcc gol_task.c -o gol_task -lpthread

clean:
	rm gol_data gol_task