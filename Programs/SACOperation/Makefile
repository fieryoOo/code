BIN1 = SAC_chhdr
BIN2 = SAC_filter
BIN3 = SAC_addall
BIN4 = SAC_dump
BIN5 = SAC_envelope
BIN6 = SAC_FTAN

LIBS = -L${HOME}/usr/lib -lfftw3 -lSacRec

cflags = -lstdc++ -std=c++11 -O3

#OBJS = SacRec.o
OBJS = 

all : $(BIN1) $(BIN2) ${BIN3} $(BIN4) $(BIN5) $(BIN6)

$(BIN1) : $(OBJS) $(BIN1).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN2) : $(OBJS) $(BIN2).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN3) : $(OBJS) $(BIN3).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN4) : $(OBJS) $(BIN4).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN5) : $(OBJS) $(BIN5).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

$(BIN6) : $(OBJS) $(BIN6).cpp
	g++ $(LIBS) $(cflags) -o $@ $^

#$(BIN2) : $(OBJS) $(BIN2).o
#	g++ $(LIBS) $(cflags) -o $(BIN2) $^

#$(BIN3) : $(OBJS) $(BIN3).o
#	g++ $(LIBS) $(cflags) -o $(BIN3) $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN1) $(BIN2) $(BIN3) $(BIN4) $(BIN5) $(BIN6)

install:
	cp $(BIN1) /home/yeti4009/bin
