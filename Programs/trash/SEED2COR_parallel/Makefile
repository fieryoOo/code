INST_DIR = $(HOME)/bin
BIN  = Seed2Cor

fflags = -O2 -Wall
cflags = -O3 -Wall

LDLIBS = -lfftw3 -lgfortran -lstdc++ -fopenmp

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(BIN)

FOBJS = Filter.o driver.o SacOps.o MakeSDB.o Pthread.o ProduceSac.o RemoveResp.o TempSpecNorm.o Whiten.o CrossCorr.o IFFT.o

$(BIN) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS) -o $(BIN)

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
