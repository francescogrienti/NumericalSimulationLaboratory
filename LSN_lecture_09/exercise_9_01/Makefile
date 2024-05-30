CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: main.exe
main.exe : main.o random.o functions.o Path.o City.o
	$(CC) Path.o City.o random.o functions.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp Path.h City.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
City.o : City.cpp City.h
	$(CC) -c City.cpp -o City.o $(CFLAGS)
Path.o : Path.cpp Path.h functions.h City.h
	$(CC) -c Path.cpp -o Path.o $(CFLAGS)
functions.o : functions.cpp functions.h random.h City.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
