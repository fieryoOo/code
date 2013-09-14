#
INST_DIR = $(HOME)/bin
BIN  = amp_length_normalize

fflags = -O2 -Wall -ffixed-line-length-none
cflags = -O2 -Wall

LDLIBS = -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = g77

CC = gcc

DBG = -g

FOBJS = amp_length_normalize.o filter4.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
