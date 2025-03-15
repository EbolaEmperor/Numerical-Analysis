all: test.o TriangularInterpolation.o
	g++ test.o TriangularInterpolation.o -o test -O2 -O3 -Ofast -std=c++20

test.o: test.cpp include/TriangularInterpolation.h
	g++ -c test.cpp -Iinclude -O2 -O3 -Ofast -std=c++20

TriangularInterpolation.o: src/TriangularInterpolation.cpp include/TriangularInterpolation.h
	g++ -c src/TriangularInterpolation.cpp -Iinclude -O2 -O3 -Ofast -std=c++20

clean:
	rm *.o test