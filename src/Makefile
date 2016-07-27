TARGET = run
CC = h5cc
CFLAGS = -std=c99 -g -Wall -fopenmp -O3

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBS = -L/lib -L/share/sw/free/gsl/2.1/lib -g -shlib -lfftw3_threads -lfftw3 -lgsl -lgslcblas -lm -fopenmp
    CFLAGS += -I/share -I/share/sw/free/gsl/2.1/include
endif
ifeq ($(UNAME_S),Darwin)
    LIBS = -L /usr/local/lib -g -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm -fopenmp
    CFLAGS += -march=native -I /usr/local/include
    export HDF5_CC = gcc-6
    export HDF5_CLINKER = gcc-6
endif

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
	./configure $(parameters)

clean:
	-rm -f *.o
	-rm -f $(TARGET)
