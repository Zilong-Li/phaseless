# detect OS architecture and add flags

HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
CXX      = g++
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
INC      = -I. -I../src -I../external -I$(HTSINC)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread

OBJS = test-main.o test-joint.o test-phaseless.o test-fastphase.o test-forward-backward.o

BINS = test-joint.bin test-phaseless.bin test-fastphase.bin test-forward-backward.bin

DEPS = ../src/phaseless.o ../src/fastphase.o ../src/admixture.o

UTILS = test-cli test-io

.PHONY: all test clean

all: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o $(DEPS)
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(DEPS) $(LDFLAGS) $(LIBS)

test: $(BINS)
	./test-joint.bin --success
	./test-phaseless.bin --success
	./test-fastphase.bin --success
	./test-forward-backward.bin --success

clean:
	rm -f *.o *.bin
