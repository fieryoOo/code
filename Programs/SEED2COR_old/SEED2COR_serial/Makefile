INST_DIR = $(HOME)/bin
BIN  = Seed2Cor

fflags = -O2 -Wall
cflags = -O3 -Wall

LDLIBS = -lfftw3 -lgfortran -lstdc++

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(BIN)

FOBJS = filter4.o driver.o SacOps.o DisAzi.o MakeSDB.o ProduceSac.o RemoveResp.o TempSpecNorm.o whiten_phamp.o CrossCorr.o ifft.o

$(BIN) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS) -o $(BIN)

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
