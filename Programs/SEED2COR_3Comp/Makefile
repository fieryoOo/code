INST_DIR = $(HOME)/bin
EXECUTABLE  = Seed2Cor

cflags = -std=c++0x -g -fopenmp #-Wall -I${HOME}/usr/include
#cflags = -std=c++0x -g -fopenmp -I/curc/tools/x_86_64/rh6/fftw/3.3.3/openmpi/1.7.4/gcc/4.8.2/double/include

LDLIBS = -L${HOME}/usr/lib -fopenmp -lfftw3 
#LDLIBS = -L/curc/tools/x_86_64/rh6/fftw/3.3.3/openmpi/1.7.4/gcc/4.8.2/double/lib -fopenmp -lfftw3 -lm

CFLAGS = $(DBG) $(cflags)

CC = g++

DBG = -g

all : $(EXECUTABLE)

FOBJS = StaPair.o CCList.o FTNorm.o CCRec.o vincenty.o CCDatabase.o SeedRec.o SacRec.o SysTools.o $(EXECUTABLE).o

$(EXECUTABLE) : $(FOBJS)
	$(CC) -o $@ $^ $(LDLIBS) $(CFLAGS)

%.o : %.cpp
	$(CC) $(cflags) -c $<

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

install_ariadne : $(EXECUTABLE)
	cp $(EXECUTABLE) $(INST_DIR)/$(EXECUTABLE)_ariadne

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
