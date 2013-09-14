#
INST_DIR = $(HOME)/bin
BIN  = whiten_rej_phamp_no_fil

fflags = -O2 -Wall -ffixed-line-length-none -m64
cflags = -O2 -Wall -m64

LDLIBS =  -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = gcc

DBG = -g

all : $(BIN)

FOBJS = driver_c.o whiten_phamp.o swapn.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
