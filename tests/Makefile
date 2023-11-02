HTSDIR   = ../inst/include/htslib-1.18
CXX      = g++
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
INC      = -I. -I../src -I../inst/include -I$(HTSDIR)
LDFLAGS  =  -L$(HTSDIR) -Wl,-rpath,$(HTSDIR)
LIBS     = -lhts -llzma -lbz2 -lm -lz -lpthread

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
