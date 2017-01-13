INST_DIR = $(HOME)/bin
BIN1 = PathAvg
BIN2 = PointAvg

cflags = -O3 -Wall -std=c++11 -lMap

LDLIBS = -lstdc++ -fopenmp

CC = g++

DBG = -g

all : $(BIN1) $(BIN2)

FOBJS =

$(BIN1) : $(FOBJS)
	$(CC) $(LDLIBS) $(cflags) $^ -o $@ $@.cpp

$(BIN2) : $(FOBJS)
	$(CC) $(LDLIBS) $(cflags) $^ -o $@ $@.cpp

%.o : %.cpp
	$(CC) $(cflags) -c $^ -o $@ 

install : $(BIN1) $(BIN2)
	install -s $(BIN1) $(INST_DIR)
	install -s $(BIN2) $(INST_DIR)

clean :
	rm -f $(BIN1) $(BIN2) core $(FOBJS)
