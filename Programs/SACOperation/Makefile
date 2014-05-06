BIN1 = chhdr

LIBS = -L${HOME}/usr/lib -lfftw3 -lSacRec

cflags = -lstdc++ -std=c++0x -O3

#OBJS = SacRec.o
OBJS = 

all : $(BIN1)

$(BIN1) : $(OBJS)
	g++ $(LIBS) $(cflags) -o $(BIN1) $(BIN1).cpp $^

#$(BIN2) : $(OBJS) $(BIN2).o
#	g++ $(LIBS) $(cflags) -o $(BIN2) $^

#$(BIN3) : $(OBJS) $(BIN3).o
#	g++ $(LIBS) $(cflags) -o $(BIN3) $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN1)
