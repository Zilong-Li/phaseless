HTSINC   = ../external
HTSLIB   = /usr/local/lib
CXX      = g++
# CXXFLAGS = -std=c++17 -Wall -O3 -g -fsanitize=address
CXXFLAGS = -std=c++17 -Wall -O3 -march=native -DNDEBUG
INC      = -I../external -I$(HTSINC) -I$(HTSLIB)
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
