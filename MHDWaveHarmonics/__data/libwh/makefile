# Compiler
CC=g++

# CFLAGS for CC
CFLAGS=-lm -std=c++11

All: o libwh.so

o:
	$(CC) -c -fpic diff.cc findharmonics.cc misfitobject.cc rkg.cc solvewave.cc spline.cc libwh.cc ApproxFundamental.cc neldermead.cc ArgSort.cc Linterp.cc -std=c++11

libwh.so:
	$(CC) -shared libwh.o diff.o findharmonics.o misfitobject.o rkg.o solvewave.o spline.o ApproxFundamental.o neldermead.o ArgSort.o Linterp.o -o libwh.so $(CFLAGS)
	
clean:
	rm -v *.o libwh.so
