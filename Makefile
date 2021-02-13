OBJECTS = apdbifurcationmain4.cpp
EXECUTABLE = bif.out
CC = g++
#FLAGS = -gencode arch=compute_30,code=sm_30
FLAGS = -std=c++11

VARIABLE1 = ibarcafac
VARIABLE2 = nacafac
VARIABLE3 = ikrfac
VARIABLE4 = iksfac

numpcl = 1
numvar1 = 1
numvar2 = 1
BEATS = 2000
REMOVEBEATS = 1500
PRECTYPE = double

BIFMACROS = -D numpcl=${numpcl} -D numvar1=${numvar1} -D numvar2=${numvar2} -D BEATS=${BEATS} -D REMOVEBEATS=${REMOVEBEATS} -D PRECTYPE=${PRECTYPE}

DV_MAX = 0.1

default: all

all: TTbif 

TTbif:
	$(CC) $(FLAGS) $(OBJECTS) -o TT$(EXECUTABLE) $(BIFMACROS) -D TT -D VARIABLE1=${VARIABLE1} -D VARIABLE2=${VARIABLE2} -D VARIABLE3=${VARIABLE3} -D VARIABLE4=${VARIABLE4} -D stimulus=-26.0 -D stimduration=2.0 -D EPI -D xiaodong

clean:
	-rm -f *$(BIFEXECUTABLE) *$(RESTEXECUTABLE)

