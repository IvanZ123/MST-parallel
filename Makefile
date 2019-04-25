test: MST.cpp
	g++ -g -o mst -fopenmp MST.cpp
stam: MST.cpp
	icc -o mst -qopenmp MST.cpp
clean:
	rm -f *.o mst
