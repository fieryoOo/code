#
INST_DIR = ../../bin
BIN  = surfsyn

#fflags =  -Wall -ffixed-line-length-none
fflags =  -Wall -O2 -ffixed-line-length-none
cflags =  -Wall -O2
#fflags =  -Wall  -O2 -m32 -ffixed-line-length-none

LDLIBS = -L/usr/lib64/libg2c.so.0

FFLAGS = $(DBG) $(fflags)

FC = gfortran
CC = gcc

DBG = 

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FOBJS = surfsyn.o source.o surfread.o intpol.o pha.o  \
	syn.o FFT.o calspecr.o calspecl.o new_tapwin1.o namer.o \
	force.o angles2tensor.o azd.o wsac_.o rw_sac.o \
	atracer.o rbimod.o read_rect_model.o tracer.o vect.o

$(BIN) : $(FOBJS)
	$(FC)  $(FFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
