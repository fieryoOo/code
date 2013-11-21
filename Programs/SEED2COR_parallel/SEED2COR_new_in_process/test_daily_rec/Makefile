BIN = Seed2amph

OBJS = driver.o SysTools.o DailyRec.o

LIBS = -L${HOME}/usr/lib -lfftw3 -lpthread

all : $(BIN)

$(BIN) : $(OBJS)
	g++ -O3 $(LIBS) -o $@ $^


clean :
	rm -f $(BIN) $(OBJS)
