BIN1 = Seed2amph

BIN2 = Seed2fSAC

OBJS = SysTools.o DailyRec.o SacRec.o

LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread -lstdc++ -fopenmp

cflags = -std=c++0x -O3 -fopenmp

all : $(BIN1) $(BIN2)

$(BIN1) : $(OBJS) $(BIN1).o
	g++ $(LIBS) -o $@ $^

$(BIN2) : $(OBJS) $(BIN2).o
	g++ $(LIBS) -o $@ $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean :
	rm -f $(BIN1) $(BIN1).o $(BIN2) $(BIN2).o $(OBJS)
