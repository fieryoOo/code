BIN = test

LIBCC = -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7 -I/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7/include

LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread -fopenmp 

OBJS = driver.o SacRec.o

cflags = -lstdc++ -std=c++0x -O3

all : $(BIN)

$(BIN) : $(OBJS)
	g++ $(LIBS) $(cflags) -o $(BIN) $^

%.o : %.cpp
	g++ $(LIBS) $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN)
