CC = g++
CFLAGS = -O3 --std=c++11
all: prova.exe
prova.exe : prova.o random.o
	$(CC) random.o prova.o -o prova.exe $(AFLAGS)
prova.o : prova.cpp random.h
	$(CC) -c prova.cpp -o prova.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o prova.exe seed.out
