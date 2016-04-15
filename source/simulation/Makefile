TARGET = run
LIBS = -L /usr/local/lib -lfftw3_omp -lfftw3 -lm -fopenmp
CC = h5cc
CFLAGS = -std=c99 -fopenmp -O2 -g -Wall -march=native -I /usr/local/include
export HDF5_CC = gcc-5
export HDF5_CLINKER = gcc-5

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
	-rm -f $(TARGET)
