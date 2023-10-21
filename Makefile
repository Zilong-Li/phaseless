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
OBJS     = src/main.o src/phaseless.o src/fastphase.o src/admixture.o src/utils.o
BINS     = phaseless
FLOAT    = 0

ifeq ($(strip $(FLOAT)),1)
  $(info "use float in phaseless!")
  CXXFLAGS += -DUSE_FLOAT
endif

.PHONY: all clean

all: $(BINS)

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

$(BINS): $(OBJS)
	${CXX} ${CXXFLAGS} -o $@ $(OBJS) ${INC} $(LIBS) $(LDFLAGS)

clean:
	rm -f $(BINS) $(OBJS)

impute:
	./phaseless impute -g data/bgl.gz -c 10 -n 4 -S -i 100

impute2:
	./phaseless impute -g data/all.bgl.gz -c 10 -n 4 -S -i 100

joint:
	./phaseless -pr --pfile impute.P --rfile impute.recomb joint -g data/bgl.gz -c 10 -k 3 -n 4 -S

joint2:
	./phaseless -pr --pfile impute.P --rfile impute.recomb joint -g data/all.bgl.gz -c 10 -k 3 -n 4 -S
