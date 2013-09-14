#
INST_DIR = $(HOME)/bin
BIN  = normalize_sac_by_noise

fflags = -O2 -Wall -ffixed-line-length-none
cflags = -O2 -Wall

LDLIBS = -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = g77

CC = gcc

DBG = -g

FOBJS = normalize_sac_by_noise.o filter4.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
