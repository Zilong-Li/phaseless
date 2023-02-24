# detect OS architecture and add flags

HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
CXX      = g++
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
INC      = -I. -I../src -I../external -I$(HTSINC)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread

OBJS = test-main.o test-io.o test-phaseless.o test-fastphasek2.o test-fastphasek4.o

BINS = test-io.bin test-phaseless.bin test-fastphasek2.bin test-fastphasek4.bin

.PHONY: all test clean

all: test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

test: $(BINS) $(OBJS)
	./test-io.bin --success
	./test-phaseless.bin --success
	./test-fastphasek2.bin --success
	./test-fastphasek4.bin --success

clean:
	rm -f *.o *.bin
