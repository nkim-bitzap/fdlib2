CPP = g++
LD = g++
AR = ar

CPPFLAGS = `pkg-config --cflags eigen3` -DGNU_GCC -c -std=c++11
ARFLAGS = rcs

all: example

example: libfdlib.a example.o
	$(LD) example.o -o example -Wl,-L,.,-lfdlib,-lpthread

example.o: example.cpp
	$(CPP) $(CPPFLAGS) example.cpp

libfdlib.a: algorithm.o forcedensity.o
	$(AR) $(ARFLAGS) libfdlib.a algorithm.o forcedensity.o

forcedensity.o: forcedensity.cpp forcedensity.h
	$(CPP) $(CPPFLAGS) forcedensity.cpp

algorithm.o: algorithm.cpp algorithm.h
	$(CPP) $(CPPFLAGS) algorithm.cpp

.PHONY: clean

clean:
	rm -f libfdlib.a algorithm.o forcedensity.o \
	example.o example
