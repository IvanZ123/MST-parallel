all: mst.o
	g++ -o  mst mst.o
mst.o: MST.cpp 
	g++ -c -o mst.o MST.cpp
test: MST.cpp
	g++ -g -o mst -fopenmp MST.cpp
clean:
	rm -f *.o mst
