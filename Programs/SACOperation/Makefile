BIN1 = chhdr
BIN2 = filter
BIN3 = SAC_addall

LIBS = -L${HOME}/usr/lib -lfftw3 -lSacRec

cflags = -lstdc++ -std=c++0x -O3

#OBJS = SacRec.o
OBJS = 

all : $(BIN1) $(BIN2) ${BIN3}

$(BIN1) : $(OBJS) $(BIN1).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN2) : $(OBJS) $(BIN2).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN3) : $(OBJS) $(BIN3).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

#$(BIN2) : $(OBJS) $(BIN2).o
#	g++ $(LIBS) $(cflags) -o $(BIN2) $^

#$(BIN3) : $(OBJS) $(BIN3).o
#	g++ $(LIBS) $(cflags) -o $(BIN3) $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN1) $(BIN2) $(BIN3)
