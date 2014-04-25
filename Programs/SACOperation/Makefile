BIN1 = chhdr

BIN2 = TestMerge

BIN3 = TestResampleRmRESP

LIBCC = -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7 -I/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7/include

LIBS = -L${HOME}/usr/lib -lfftw3

cflags = -lstdc++ -std=c++0x -O3

OBJS = SacRec.o SysTools.o

all : $(BIN1)

$(BIN1) : $(OBJS) $(BIN1).o
	g++ $(LIBS) $(cflags) -o $(BIN1) $^

#$(BIN2) : $(OBJS) $(BIN2).o
#	g++ $(LIBS) $(cflags) -o $(BIN2) $^

#$(BIN3) : $(OBJS) $(BIN3).o
#	g++ $(LIBS) $(cflags) -o $(BIN3) $^

%.o : %.cpp
	g++ $(LIBS) $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN1).o $(BIN2).o $(BIN3).o $(BIN1) $(BIN2) $(BIN3)
