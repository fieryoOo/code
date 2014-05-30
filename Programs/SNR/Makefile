BIN1 = ComputeSNR

BIN2 = TestFFT

OBJs = SacRec.o

cflags = -std=c++0x -O3

LIBs = -lfftw3

all : $(BIN1) $(BIN2)

$(BIN1) : $(BIN1).o $(OBJs)
	g++ $(cflags) $(LIBs) -o $@ $^

$(BIN2) : $(BIN2).o $(OBJs)
	g++ $(cflags) $(LIBs) -o $@ $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean :
	rm -f $(BIN) $(OBJs)
