INST_DIR = $(HOME)/bin
BIN1 = LocatePrecursor
BIN2 = PredictPrecursor_PointSource
BIN3 = PredictPrecursor_SourceGroup
BIN4 = PrecursorSNR_SourceGroup

cflags = -O3 -std=c++0x

LDLIBS = -fopenmp -lSacRec -lMap -lfftw3

CFLAGS = $(DBG) $(cflags) #-gdwarf-2 -gstrict-dwarf

CC = g++

DBG =

all : $(BIN1) $(BIN2) $(BIN3) $(BIN4)

FOBJS1 = DisAzi.o SacOps.o PathAverage.o $(BIN1).o 

FOBJS2 = DisAzi.o SacOps.o PathAverage.o $(BIN2).o

FOBJS3 = DisAzi.o SacOps.o PathAverage.o $(BIN3).o

FOBJS4 = DisAzi.o SacOps.o PathAverage.o $(BIN4).o

$(BIN1) : $(FOBJS1)
	$(CC) $(LDLIBS) $(CFLAGS) $^ -o $@

$(BIN2) : $(FOBJS2)
	$(CC) $(LDLIBS) $(CFLAGS) $^ -o $@

$(BIN3) : $(FOBJS3)
	$(CC) $(LDLIBS) $(CFLAGS) $^ -o $@

$(BIN4) : $(FOBJS4)
	$(CC) $(LDLIBS) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(LDLIBS) $(CFLAGS) -c $^

install : $(BIN1) $(BIN2) $(BIN3)
	install -s $(BIN1) $(INST_DIR)
	install -s $(BIN2) $(INST_DIR)
	install -s $(BIN3) $(INST_DIR)
	install -s $(BIN4) $(INST_DIR)

clean :
	rm -f $(BIN1) $(BIN2) $(BIN3) $(BIN4) $(FOBJS1) $(FOBJS2) $(FOBJS3) $(FOBJS4)
