#
INST_DIR = $(HOME)/bin
BIN  = pred_from_transfer_function

fflags = -O2 -Wall -ffixed-line-length-none -m32
cflags = -O2 -Wall -m32

LDLIBS =   -L/home/linf/lib -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = g77

CC = gcc

DBG =

FOBJS =  rmresponse.o driver_c.o


$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -static -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
