BIN = test

LIBS = -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7 -I/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7/include

OBJS = driver.o SacRec.o

cflags = -lstdc++ -std=c++0x -O3

all : $(BIN)

$(BIN) : $(OBJS)
	g++44 $(cflags) -o $(BIN) $^

%.o : %.cpp
	g++44 $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN)
