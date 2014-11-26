#
INST_DIR = $(HOME)/bin
BIN  = aftani_c_pgl_amp

fflags = -O2 -Wall -ffixed-line-length-132
#cflags = -O3 -Wall  -I"/home/tianye/code/Programs/head/"
#cflags = -O3 -Wall  -I/usr/include -I/home/tianye/code/Programs/head/
cflags = -O3 -Wall  -I/usr/include -I/projects/yeti4009/code/Programs/head/

#LDLIBS =   -L/home/barmin/FFTW_64/lib -lfftw3
LDLIBS = -L/usr/lib64 -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = gcc

DBG =

FOBJS = driver_c.o swapn.o aftanpg.o aftanipg.o ftfilt.o Filter.o fmax.o \
        taper.o trigger.o tapers.o tgauss.o mspline.o pred_cur.o misc.o \
        phtovel.o whiten_phamp.o get_snr.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -static -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
