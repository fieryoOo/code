INST_DIR = $(HOME)/bin
BIN  = LocatePrecursor

fflags = -O2 -Wall
cflags = -O3 -Wall

LDLIBS = -lstdc++ -fopenmp

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = g++

DBG = -g

all : $(BIN)

FOBJS = DisAzi.o SacOps.o PathAverage.o driver.o 

$(BIN) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS) -o $(BIN)

%.o: %.cpp
	$(CC) $(LDLIBS) $(CFLAGS) -c $^

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
