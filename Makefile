HTSDIR   = ./inst/include/htslib-1.18
CXX      = g++

# CXXFLAGS = -std=c++17 -Wall -O3 -g -fsanitize=address
# CXXFLAGS = -std=c++17 -Wall -O3 -march=native -DNDEBUG
CXXFLAGS = -std=c++17 -Wall -O3 -DNDEBUG
INC      = -I./src -I./inst/include -I$(HTSDIR)
LDFLAGS  =  -L$(HTSDIR) -Wl,-rpath,$(HTSDIR)
LIBS     =  -llzma -lbz2 -lm -lz -lpthread
OBJS     = src/main.o src/phaseless.o src/fastphase.o src/admixture.o src/utils.o
BINS     = phaseless
libhts   = $(HTSDIR)/libhts.a
FLOAT    = 0

ifeq ($(strip $(FLOAT)),1)
  $(info "use float in phaseless!")
  CXXFLAGS += -DUSE_FLOAT
endif

.PHONY: all clean 

all: $(BINS) $(libhts)

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

$(BINS): $(OBJS) $(libhts)
	${CXX} ${CXXFLAGS} -o $@ $(OBJS) ${INC} $(libhts) $(LIBS) $(LDFLAGS)

$(libhts):
	cd $(HTSDIR) && ./configure --disable-libcurl --without-libdeflate && make -j6

clean:
	rm -f $(BINS) $(OBJS)
	cd $(HTSDIR) && make clean

impute:
	./phaseless -Dr impute -g data/bgl.gz -c 10 -n 4 -S -i 100

joint:
	./phaseless -Dpr --pfile impute.P --rfile impute.recomb joint -g data/bgl.gz -c 10 -k 3 -n 4 -S -i 100

parse:
	./phaseless -Dpr parse -j joint.pars.bin -n 4 -i 100

impute2:
	./phaseless impute -g data/all.bgl.gz -c 10 -n 4 -S -i 100

joint2:
	./phaseless -Dpr --pfile impute.P --rfile impute.recomb joint -g data/all.bgl.gz -c 10 -k 3 -n 4 -S
