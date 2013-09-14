INST_DIR = $(HOME)/bin
BIN  = Fortune

fflags = -O2 -Wall
cflags = -O3 -Wall

LDLIBS = -lfftw3 -lgfortran -lstdc++

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(BIN)

FOBJS = edgelist.o  geometry.o  heap.o  main.o  memory.o  output.o  voronoi.o

$(BIN) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS) -o $(BIN)

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
