CC = mpicxx
CFLAGS = -O3 --std=c++11 -I/usr/lib/x86_64-linux-gnu/openmpi/include
AFLAGS = -larmadillo -L/usr/lib/x86_64-linux-gnu
NP = 4  # Numero di processi MPI

all: main.exe

main.exe : main.o random.o Genetics.o Path.o City.o
	$(CC) Genetics.o Path.o City.o random.o main.o -o main.exe $(AFLAGS)

main.o : main.cpp Genetics.h Path.h City.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

Genetics.o : Genetics.cpp Genetics.h Path.h random.h
	$(CC) -c Genetics.cpp -o Genetics.o $(CFLAGS)

City.o : City.cpp City.h
	$(CC) -c City.cpp -o City.o $(CFLAGS)

Path.o : Path.cpp Path.h City.h
	$(CC) -c Path.cpp -o Path.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

run: main.exe
	mpiexec --use-hwthread-cpus -np $(NP) ./main.exe

clean :
	rm *.o main.exe seed.out