INST_DIR = $(HOME)/bin
BIN  = RemoveRESP

fflags = -O2 -Wall
cflags = -O3

LDLIBS = -lfftw3 -lstdc++

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(cflags)

FC = gfortran

CC = g++

DBG = -g

FOBJS =  driver_c.o SysTools.o FreqOps.o


$(BIN) : $(FOBJS)
	$(CC) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
