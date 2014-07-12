BIN1 = Seed2Sac

OBJS = SysTools.o SeedRec.o

LIBS = -L${HOME}/usr/lib -lfftw3 -lstdc++ -fopenmp -lSacRec

cflags = -std=c++0x -g -fopenmp

all : $(BIN1)

$(BIN1) : $(OBJS) $(BIN1).o
	g++ $(LIBS) -o $@ $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean :
	rm -f $(BIN1) $(BIN1).o $(OBJS)
