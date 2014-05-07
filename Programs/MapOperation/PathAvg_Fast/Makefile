INST_DIR = $(HOME)/bin
BIN1 = PathAvg
BIN2 = PointAvg

cflags = -O3 -Wall -std=c++0x -fopenmp

LDLIBS = -lstdc++ -fopenmp

CC = g++

DBG = -g

all : $(BIN1) $(BIN2)

FOBJS = Map_Fast.o

$(BIN1) : $(FOBJS) $(BIN1).o
	$(CC) $(LDLIBS) $^ -o $@

$(BIN2) : $(FOBJS) $(BIN2).o
	$(CC) $(LDLIBS) $^ -o $@

%.o : %.cpp
	$(CC) $(cflags) -c $^ -o $@ 

install : $(BIN1) $(BIN2)
	install -s $(BIN1) $(INST_DIR)
	install -s $(BIN2) $(INST_DIR)

clean :
	rm -f $(BIN1) $(BIN1).o $(BIN2) $(BIN2).o core $(FOBJS)
