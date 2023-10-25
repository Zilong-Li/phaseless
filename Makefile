HTSDIR   = ./inst/include/htslib-1.18
CXX      = g++

# CXXFLAGS = -std=c++17 -Wall -O3 -g -fsanitize=address
# CXXFLAGS = -std=c++17 -Wall -O3 -march=native -DNDEBUG
CXXFLAGS = -std=c++17 -Wall -O3 -march=native -fPIC -DNDEBUG
INC      = -I./src -I./inst/include -I$(HTSDIR)
LDFLAGS  =  -L$(HTSDIR) $(HTSDIR)/libhts.a
LIBS     = -llzma -lbz2 -lm -lz -lpthread
OBJS     = src/main.o src/phaseless.o src/fastphase.o src/admixture.o src/utils.o
BINS     = phaseless
FLOAT    = 0

ifeq ($(strip $(FLOAT)),1)
  $(info "use float in phaseless!")
  CXXFLAGS += -DUSE_FLOAT
endif

.PHONY: all clean htslib

all: $(BINS)

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

$(BINS): $(OBJS) htslib
	${CXX} ${CXXFLAGS} -o $@ $(OBJS) ${INC} $(LIBS) $(LDFLAGS)

htslib:
	cd $(HTSDIR) && ./configure --disable-libcurl && make -j10

clean:
	rm -f $(BINS) $(OBJS)
	cd $(HTSDIR) && make clean

impute:
	./phaseless impute -g data/bgl.gz -c 10 -n 4 -S -i 100

joint:
	./phaseless -Dpr --pfile impute.P --rfile impute.recomb joint -g data/bgl.gz -c 10 -k 3 -n 4 -S -i 100

parse:
	./phaseless -Dpr parse -j joint.pars.bin -n 4 -i 100

impute2:
	./phaseless impute -g data/all.bgl.gz -c 10 -n 4 -S -i 100

joint2:
	./phaseless -Dpr --pfile impute.P --rfile impute.recomb joint -g data/all.bgl.gz -c 10 -k 3 -n 4 -S
