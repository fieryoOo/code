INST_DIR = $(HOME)/src/local/bin
BIN  = Seed2Cor

fflags = -O2 -Wall
cflags = -O3 -Wall

LDLIBS = -lm -lfftw3 -lstdc++ -lpthread -fopenmp
TMP = -lgfortran

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(BIN)

FOBJS = driver.o SacOps.o DisAzi.o SysTools.o MakeSDB.o Pthread.o ProduceSac.o RemoveResp.o TempSpecNorm.o Whiten.o CrossCorr.o FreqOps.o

$(BIN) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS) -o $(BIN)

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
