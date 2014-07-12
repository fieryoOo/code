INST_DIR = $(HOME)/bin
EXECUTABLE  = Seed2Cor

cflags = -std=c++0x -O3 -Wall #-I${HOME}/usr/include

LDLIBS = -L${HOME}/usr/lib -lfftw3 -lpthread
LTMP = -lstdc++ -fopenmp

CFLAGS = $(DBG) $(cflags)

CC = g++

DBG = -g

all : $(EXECUTABLE)

FOBJS = CCDatabase.o SeedRec.o SacRec.o SysTools.o $(EXECUTABLE).o

$(EXECUTABLE) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) -o $@ $^

%.o : %.cpp
	$(CC) $(cflags) -c $<

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

install_ariadne : $(EXECUTABLE)
	cp $(EXECUTABLE) $(INST_DIR)/$(EXECUTABLE)_ariadne

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
