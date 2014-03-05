BIN = Seed2amph

OBJS = $(BIN).o SysTools.o DailyRec.o SacRec.o

LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread -lstdc++ -fopenmp

cflags = -std=c++0x -fopenmp -O3

all : $(BIN)

$(BIN) : $(OBJS)
	g++ $(LIBS) -o $@ $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean :
	rm -f $(BIN) $(OBJS)
