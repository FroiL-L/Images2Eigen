all: main

main: ReadImage.o image.o jacobi.o
	g++ -o main.exe main.cpp ReadImage.o image.o jacobi.o rgb.o

ReadImage.o: image.o rgb.o
	g++ -c ReadImage.cpp -o ReadImage.o

image.o: rgb.o
	g++ -c image.cpp -o image.o

rgb.o:
	g++ -c rgb.cpp -o rgb.o

jacobi.o:
	gcc -c jacobi.c -o jacobi.o -lm

clean:
	rm main.exe main.o

cleanall:
	rm main.exe *.o
