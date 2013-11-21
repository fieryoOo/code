INST_DIR = $(HOME)/bin
BIN1 = LocatePrecursor
BIN2 = PredictPrecursor

cflags = -O3 -Wall

LDLIBS = -fopenmp

CFLAGS = $(DBG) $(cflags)

CC = g++

DBG = -g

all : $(BIN1) $(BIN2)

FOBJS1 = DisAzi.o SacOps.o PathAverage.o $(BIN1).o 

FOBJS2 = DisAzi.o SacOps.o PathAverage.o $(BIN2).o

$(BIN1) : $(FOBJS1)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS1) -o $(BIN1)

$(BIN2) : $(FOBJS2)
	$(CC) $(LDLIBS) $(CFLAGS) $(FOBJS2) -o $(BIN2)

%.o: %.cpp
	$(CC) $(LDLIBS) $(CFLAGS) -c $^

install : $(BIN1) $(BIN2)
	install -s $(BIN1) $(INST_DIR)
	install -s $(BIN2) $(INST_DIR)

clean :
	rm -f $(BIN1) $(BIN2) $(FOBJS1) $(FOPBS2)
