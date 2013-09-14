INST_DIR = $(HOME)/bin
BIN  = csta_atn_aplf_term

fflags = -O2 -Wall
cflags = -O2 -Wall

LDLIBS =

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(BIN)

#FOBJS = driver_c.o whiten_phamp.o fft_phamp.o swapn.o
FOBJS =  driver.o grdt_lplc.o STACK_dayflag_csta.o bin_hex.o sac_io.o dist_azimuth.o local_amp_grdt.o ph_amp_map.o Combiner.o

$(BIN) : $(FOBJS)
	$(CC) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)

