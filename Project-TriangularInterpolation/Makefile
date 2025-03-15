all: test.o testFilter.o TriangularInterpolation.o
	g++ test.o TriangularInterpolation.o -o test -O2 -O3 -Ofast -std=c++20
	g++ testFilter.o TriangularInterpolation.o -o testFilter -O2 -O3 -Ofast -std=c++20

test.o: test.cpp include/TriangularInterpolation.h
	g++ -c test.cpp -Iinclude -O2 -O3 -Ofast -std=c++20

testFilter.o: testFilter.cpp include/TriangularInterpolation.h
	g++ -c testFilter.cpp -Iinclude -O2 -O3 -Ofast -std=c++20

TriangularInterpolation.o: src/TriangularInterpolation.cpp include/TriangularInterpolation.h
	g++ -c src/TriangularInterpolation.cpp -Iinclude -O2 -O3 -Ofast -std=c++20

clean:
	rm *.o test testFilter