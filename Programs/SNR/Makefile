BIN1 = ComputeSNR

BIN2 = TestFFT

BIN3 = TestSNRofNoise

BIN4 = ComputeSingleSNR

OBJs = SacRec.o

cflags = -std=c++0x -O3

LIBs = -lfftw3

all : $(BIN1) $(BIN2) $(BIN3) $(BIN4)

$(BIN1) : $(BIN1).o $(OBJs)
	g++ $(cflags) $(LIBs) -o $@ $^

$(BIN2) : $(BIN2).o $(OBJs)
	g++ $(cflags) $(LIBs) -o $@ $^

$(BIN3) : $(BIN3).o $(OBJs)
	g++ $(cflags) $(LIBs) -o $@ $^

$(BIN4) : $(BIN4).o $(OBJs)
	g++ $(cflags) $(LIBs) -o $@ $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean :
	rm -f $(BIN1) $(BIN1).o $(BIN2) $(BIN2).o $(BIN3) $(BIN3).o $(BIN4) $(BIN4).o $(OBJs)
