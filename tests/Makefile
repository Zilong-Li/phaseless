HTSDIR   = ../inst/include/htslib-1.18
CXX      = g++
CXXFLAGS = -std=c++17 -Wall -O3 -march=native -MMD -MP
INC      = -I../src -I../inst/include -I$(HTSDIR)
# LDFLAGS  =  -L$(HTSDIR) -Wl,-rpath,$(HTSDIR)
LIBS     = $(HTSDIR)/libhts.a -llzma -lbz2 -lm -lz -lpthread
CUDA     ?= 0

ifeq ($(strip $(CUDA)),1)
  LIBS += -lcudart
endif

OBJS = test-main.o test-joint.o test-admixture.o test-fastphase.o test-utils.o test-io.o test-cli.o test-cuda.o

BINS = test-joint.bin test-admixture.bin test-fastphase.bin test-utils.bin test-io.bin test-cli.bin test-cuda.bin

DEPS = ../src/libsrc.a

UTILS = test-cli test-io

.PHONY: all test clean

all: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o $(DEPS)
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(DEPS) $(LDFLAGS) $(LIBS) $(INC)

test: $(BINS)
	# ./test-admixture.bin --success
	# ./test-fastphase.bin --success
	./test-utils.bin --success
	./test-joint.bin --success
	./test-io.bin --success
	./test-cli.bin --success
	./test-cuda.bin --success

clean:
	rm -f *.o *.bin *.d

-include $(OBJS:.o=.d)
