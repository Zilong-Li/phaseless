HTSINC   = ./external
HTSLIB   = /usr/local/lib
CXX      = g++
# CXXFLAGS = -std=c++17 -Wall -O3 -g -fsanitize=address
# CXXFLAGS = -std=c++17 -Wall -O3 -march=native -DNDEBUG
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
INC      = -I../external -I$(HTSINC) -I$(HTSLIB)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread
# OBJS     = $(patsubst src/%.cpp, src/%.o, $(wildcard *.cpp))
OBJS     = src/admix.o src/impute.o src/parse.o src/convert.o
BINS     = phaseless

.PHONY: all clean

all: $(BINS)

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

$(BINS): $(OBJS) src/phaseless.o
	${CXX} ${CXXFLAGS} -o $@ src/phaseless.o $< ${INC} $(LIBS) $(LDFLAGS)

clean:
	rm -f $(BINS) src/*.o
