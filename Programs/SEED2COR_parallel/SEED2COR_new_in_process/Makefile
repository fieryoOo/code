INST_DIR = $(HOME)/bin
EXECUTABLE  = Seed2Cor

fflags = -O2 -Wall
cflags = -O3 -Wall
cflag2 = -O3 -Wall -I${HOME}/usr/include

LDLIBS = -L${HOME}/usr/lib -lfftw3 -lstdc++ -lpthread -fopenmp 
TMP = -lgfortran

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++
ICC = icpc

DBG = -g

all : $(EXECUTABLE)

FOBJS = driver.o SacOps.o DisAzi.o SysTools.o MakeSDB.o Pthread.o ProduceSac.o RemoveResp.o TempSpecNorm.o Whiten.o CrossCorr.o FreqOps.o

$(EXECUTABLE) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $^ -o $@

Whiten.o : Whiten.c
	$(CC) $(cflag2) -c $<

FreqOps.o : FreqOps.c
	$(CC) $(cflag2) -c $<

%.o : %.c
	$(CC) $(cflags) -c $<

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

install_ariadne : $(EXECUTABLE)
	cp $(EXECUTABLE) $(INST_DIR)/$(EXECUTABLE)_ariadne

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
