#
INST_DIR = $(HOME)/bin
#BIN  = filter4_test
BIN  = filter4

fflags = -O2 -Wall -ffixed-line-length-none
cflags = -O2 -Wall

#LDLIBS =   -L/home/zheng/FFTW/lib -lfftw3

LDLIBS = -lfftw3 -lgfortran

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = gcc

DBG = -g

all : $(BIN)

FOBJS = driver_c.o filter4.o swapn.o

driver_c.o: driver_c.c
	$(CC) -c driver_c.c -o driver_c.o -m64

filter4.o: filter4.f
	$(FC) -c -o filter4.o filter4.f -m64

swapn.o: swapn.c
	$(CC) -c -o swapn.o swapn.c -m64


$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) -m64

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
