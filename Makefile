LIBS="-lm"

GCCF="-g0 -O -Wall" 
CXXF="-w0 -g0 -O -std strict_ansi"
KCCF="-g -O --strict"

cc:
	make CXX="c++ -Itoolbox/src -std=c++11" CXXFLAGS=${GCCF} all

gcc:
	make CXX="g++ -Itoolbox/src -std=c++11" CXXFLAGS=${GCCF} all

all:
	$(CXX) tests/opoly_test1.cc -o t1 $(LIBS)

clean:
	rm -f t1 t2
