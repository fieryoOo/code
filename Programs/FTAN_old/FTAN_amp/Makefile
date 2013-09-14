#
INST_DIR = $(HOME)/bin
BIN  = aftani_c_pgl_amp

fflags = -O2 -Wall -ffixed-line-length-none 
cflags = -O2 -Wall 

LDLIBS =   -L/home/barmin/FFTW_64/lib -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = g77

CC = gcc

DBG =

FOBJS = driver_c.o swapn.o aftanpg.o aftanipg.o ftfilt.o fmax.o taper.o \
	trigger.o tapers.o tgauss.o mspline.o pred_cur.o misc.o phtovel.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -static -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
