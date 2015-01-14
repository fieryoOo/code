

BIN = Rename

CC = g++

LIBs = -fopenmp -lSysTools # -L${HOME}/usr/lib

cflags = -std=c++11 -O3 # -fopenmp

all : $(BIN)

OBJs = Rename.o

$(BIN) : $(OBJs)
	$(CC) $(LIBs) $(cflags) -o $@ $^

%.o : %.cpp
	$(CC) $(cflags) -c $<

clean :
	rm -f $(BIN) $(OBJs) core.*
