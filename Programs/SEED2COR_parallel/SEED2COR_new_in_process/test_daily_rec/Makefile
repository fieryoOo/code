BIN = Seed2amph

OBJS = driver.o SysTools.o DailyRec.o SacRec.o

LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread

cflags = -lstdc++ -std=c++0x -fopenmp -O3

all : $(BIN)

$(BIN) : $(OBJS)
	g++ $(cflags) $(LIBS) -o $@ $^

%.o : %.cpp
	g++ $(LIBS) $(cflags) -c $^

clean :
	rm -f $(BIN) $(OBJS)
