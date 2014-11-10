INST_DIR = $(HOME)/usr/bin
BIN  = justCOR

fflags = -O2 -ffixed-line-length-none -m32
cflags = -O2 -m32
#fflags = -fast -Mextend -Msave
#fflags =  -O -e

# if you don't use shared libraries comment the next line
#shlib = -R/opt/SUNWspro/lib

FFLAGS = $(DBG)  $(fflags) $(shlib)
CFLAGS = $(DBG)  $(cflags)

LDLIBS = -L$(HOME)/usr/lib -lfftw3
FC = g77
CC = gcc
DBG = 
FOBJS = dcommon.o dmultifft.o justCOR.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -static -o $(BIN) $(LDLIBS)

install : $(BIN)
	install  -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
