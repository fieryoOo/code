BIN = StalstFromSACs

#LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread -fopenmp -lstdc++
LIBS = -lfftw3 -lpthread -fopenmp -lstdc++ -lSacRec

cflags = -std=c++0x -O3

OBJS = 

all : $(BIN)

$(BIN) : $(OBJS)
	g++ $(LIBS) $(cflags) -o $(BIN) $(BIN).cpp $^

%.o : %.cpp
	g++ $(LIBS) $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN)
