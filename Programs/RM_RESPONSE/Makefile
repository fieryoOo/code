INST_DIR = $(HOME)/bin
BIN  = RemoveRESP

fflags = -O2 -Wall
cflags = -O3 -Wall

LDLIBS = -lfftw3 -lgfortran -lstdc++

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

FOBJS =  driver_c.o SysTools.o FreqOps.o


$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -static -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
