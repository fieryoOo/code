#
INST_DIR = ../../bin
BIN  = mod2tabs

fflags = -O -Wall -ffixed-line-length-none
#fflags =  -Wall  -O2 -m32 -ffixed-line-length-none

FFLAGS = $(DBG) $(fflags)

FC = gfortran

DBG =

FOBJS =  mod2tabs.o

$(BIN) : $(FOBJS)
	$(FC)  $(FFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
