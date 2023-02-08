# detect OS architecture and add flags

HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
EIGEN    = /usr/local/include
CXX      = g++
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
INC      = -I. -I../src -I$(HTSINC) -I$(EIGEN)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread

OBJS = test-main.o test-io.o test-fastphasek4.o

BINS = test-io.bin test-fastphasek4.bin

.PHONY: all test clean

all: $(BINS) $(OBJS)

test: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.bin
