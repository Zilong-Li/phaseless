HTSINC   = /usr/local/include
HTSLIB   = /usr/local/lib
CXX      = g++
# CXXFLAGS = -std=c++17 -Wall -O3 -g -fsanitize=address
# CXXFLAGS = -std=c++17 -Wall -O3 -march=native -DNDEBUG
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
INC      = -I./src -I./external -I$(HTSINC) -I$(HTSLIB)
LDFLAGS  =  -L$(HTSLIB) -Wl,-rpath,$(HTSLIB)
LIBS     = -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread
# OBJS     = $(patsubst src/%.cpp, src/%.o, $(wildcard *.cpp))
OBJS     = src/phaseless.o src/fastphase.o src/admixture.o src/utils.o
BINS     = phaseless

.PHONY: all clean

all: $(BINS)

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

$(BINS): $(OBJS)
	${CXX} ${CXXFLAGS} -o $@ $(OBJS) ${INC} $(LIBS) $(LDFLAGS)

clean:
	rm -f $(BINS) $(OBJS)
