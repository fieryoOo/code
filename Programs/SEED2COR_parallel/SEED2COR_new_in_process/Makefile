INST_DIR = $(HOME)/bin
EXECUTABLE  = Seed2Cor

cflags = -O3 -Wall
cflag2 = -O3 -Wall -I${HOME}/usr/include

LDLIBS = -L${HOME}/usr/lib -lfftw3 -lpthread
LTMP = -lstdc++ -fopenmp

CFLAGS = $(DBG) $(cflags)

CC = g++

DBG = -g

all : $(EXECUTABLE)

FOBJS = driver.o CCDatabase.o

$(EXECUTABLE) : $(FOBJS)
	$(CC) $(LDLIBS) $(CFLAGS) $^ -o $@

Whiten.o : Whiten.cpp
	$(CC) $(cflag2) -c $<

FreqOps.o : FreqOps.cpp
	$(CC) $(cflag2) -c $<

%.o : %.cpp
	$(CC) $(cflags) -c $<

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

install_ariadne : $(EXECUTABLE)
	cp $(EXECUTABLE) $(INST_DIR)/$(EXECUTABLE)_ariadne

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
