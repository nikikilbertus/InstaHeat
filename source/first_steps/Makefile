TARGET = first_steps
LIBS = -L /usr/local/lib -lfftw3_omp -lfftw3 -lm -fopenmp -lprofiler
CC = gcc-5
CFLAGS = -std=c99 -fopenmp -O3 -g -Wall -I /usr/local/include

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
	-rm -f file*.txt