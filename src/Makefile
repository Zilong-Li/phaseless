HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
EIGEN    = /usr/local/include
CXX      = g++
# CXXFLAGS = -std=c++17 -Wall -O3 -g -fsanitize=address
CXXFLAGS = -std=c++17 -Wall -O3 -g -march=native
INC      =  -I$(HTSINC) -I$(EIGEN)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread
# OBJS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
BINS     = $(patsubst %.cpp, %, $(wildcard *.cpp))

.PHONY: all clean

all: $(BINS)

%: %.cpp
	${CXX} ${CXXFLAGS} -o $@ $< ${INC} $(LIBS) $(LDFLAGS)

clean:
	rm -f $(BINS)
