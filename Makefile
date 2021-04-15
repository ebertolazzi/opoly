LIBS="-lm"

GCCF="-g0 -O -Wall" 
CXXF="-w0 -g0 -O -std strict_ansi"
KCCF="-g -O --strict"

cc:
	make CXX="c++" CXXFLAGS=${GCCF} all

gcc:
	make CXX="g++" CXXFLAGS=${GCCF} all

cxx:
	make CXX="cxx" CXXFLAGS=${CXXF} all

kcc:
	make CXX="KCC" CXXFLAGS=${KCCF} all

OBJS = opoly.test.o

all:	$(OBJS)
	$(CXX) $(OBJS) -o t $(LIBS)

opoly.test.o: opoly.test.cc opoly.hh
	$(CXX) $(CXXFLAGS) -c opoly.test.cc

clean:
	rm -f t $(OBJS) *~ pch/opolyPPC++ opoly.µ
	rm -rf "opoly Data"


tgz:
	rm -rf opoly
	mkdir              opoly
	cp Makefile        opoly
	cp opoly.hh        opoly
	cp opoly.test.cc   opoly
	cp license         opoly
	tar -cvf opoly.tar opoly
	gzip -9 opoly.tar
	mv opoly.tar.gz opoly.tgz
