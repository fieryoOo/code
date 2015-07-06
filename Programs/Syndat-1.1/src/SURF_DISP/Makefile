#
INST_DIR = ../../bin
BIN  = SURF_DISP

fflags =  -Wall -O2 -ffixed-line-length-none
#fflags =  -Wall  -O2 -m32 -ffixed-line-length-none

FFLAGS = $(DBG) $(fflags)

FC = gfortran

DBG = 

FOBJS =  SURF_DISP.o flat1.o surfa.o calcul_deep.o init_deep.o

$(BIN) : $(FOBJS)
	$(FC)  $(FFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
