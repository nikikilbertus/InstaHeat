TARGET = run
LIBS = -L /usr/local/lib -g -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm -fopenmp
CC = h5cc
CFLAGS = -std=c99 -g -Wall -fopenmp -O2 -march=native -I /usr/local/include
export HDF5_CC = gcc-6
export HDF5_CLINKER = gcc-6

.PHONY: default all clean

default: config $(TARGET)

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

config:
	./configure

clean:
	-rm -f *.o
	-rm -r run.dSYM
	-rm -f $(TARGET)
