BIN = BSpline1D

cflags = -std=c++11 
libs = -lmultivariate-splines-1-2

all : $(BIN)

$(BIN) : $(BIN).o
	g++ $(cflags) $(libs) -O3 -o $@ $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean :
	rm -f $(BIN) $(BIN).o
