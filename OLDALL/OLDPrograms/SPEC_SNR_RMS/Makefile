#
INST_DIR = $(HOME)/bin
BIN  = lf_spec_snr_rms_fast

fflags = -O2 -Wall -ffixed-line-length-none
cflags = -O2 -Wall

LDLIBS =   -L/home/barmin/FFTW/lib -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = g77

CC = gcc

DBG = -g

FOBJS = lf_spec_snr_rms_fast.o filter4.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -static -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
