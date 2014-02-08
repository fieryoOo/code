BIN1 = TestBasic

BIN2 = TestMerge

LIBCC = -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7 -I/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7/include

LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread -fopenmp 

cflags = -lstdc++ -std=c++0x -O3

all : $(BIN1) $(BIN2)

$(BIN1) : SacRec.o $(BIN1).o
	g++ $(LIBS) $(cflags) -o $(BIN1) $^

$(BIN2) : SacRec.o $(BIN2).o
	g++ $(LIBS) $(cflags) -o $(BIN2) $^

%.o : %.cpp
	g++ $(LIBS) $(cflags) -c $^

clean:
	rm -f SacRec.o $(BIN1).o $(BIN2).o $(BIN1) $(BIN2)
