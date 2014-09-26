CC = gcc
CFLAGS = -g -lm -Wall -std=c99
OBJECTS = scatter.o 
all : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o scatter


scatter: scatter.c
	$(CC) $(CFLAGS) -c scatter.c
