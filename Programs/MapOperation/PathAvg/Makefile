INST_DIR = $(HOME)/bin
BIN  = MapAvg

cflags = -O3 -Wall -std=c++0x -fopenmp

LDLIBS = -lstdc++ -fopenmp

CC = g++

DBG = -g

all : $(BIN)

FOBJS = driver.o DisAzi.o Map.o

$(BIN) : $(FOBJS)
	$(CC) $(LDLIBS) $(FOBJS) -o $(BIN)

%.o : %.cpp
	$(CC) $(cflags) -c $^ -o $@ 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
