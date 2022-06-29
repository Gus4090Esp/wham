CC = g++
MPICC = mpic++
CFLAGS = -Wall -Wextra -O3 -fno-strict-aliasing -std=gnu++11
LFLAGS = -lm

all: wham_colvar wham2d_colvar_vark

wham_colvar: wham_colvar.cpp histogram.cpp
	$(CC) $(CFLAGS) -c histogram.cpp
	$(CC) $(CFLAGS) -c wham_colvar.cpp
	$(CC) $(CFLAGS) histogram.o wham_colvar.o -o wham_colvar $(LFLAGS)

wham2d_colvar_vark: wham2d_colvar_vark.cpp histogram2d.cpp
	$(CC) $(CFLAGS) -c histogram2d.cpp
	$(CC) $(CFLAGS) -c wham2d_colvar_vark.cpp
	$(CC) $(CFLAGS) histogram2d.o wham2d_colvar_vark.o -o wham2d_colvar_vark $(LFLAGS)

clean:
	rm wham_colvar wham_colvar.o histogram.o wham2d_colvar_vark wham2d_colvar_vark.o histogram2d.o

