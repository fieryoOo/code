#
INST_DIR = $(HOME)/bin
BIN  =  SURF_DISP

fflags = -O2 -Wall -ffixed-line-length-none -m64
cflags = -O2 -Wall -m64

#LDLIBS = -L/home/barmin/FFTW_64/lib -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = g77

CC = gcc

DBG =

FOBJS = SURF_DISP.o flat1.o surfa.o calcul_deep.o init_deep.o

$(BIN) : $(FOBJS)
	$(FC) $(CFLAGS) $(FOBJS) -static -o $(BIN)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
