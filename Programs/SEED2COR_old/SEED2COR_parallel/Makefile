INST_DIR = $(HOME)/bin
EXECUTABLE  = Seed2Cor

fflags = -O2 -Wall
cflags = -O3 -Wall -I${HOME}/usr/include 

LDLIBS = -lm -L${HOME}/usr/lib -lfftw3 -lstdc++ -lpthread -fopenmp 
TMP = -lgfortran

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(EXECUTABLE)

FOBJS = driver.o SacOps.o DisAzi.o SysTools.o MakeSDB.o Pthread.o ProduceSac.o RemoveResp.o TempSpecNorm.o Whiten.o CrossCorr.o FreqOps.o

$(EXECUTABLE) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS) -o $(EXECUTABLE)

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
