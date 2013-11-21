BIN = test

OBJS = driver.o SacRec.o

cflags = -lstdc++ -std=c++11 -O3

all : $(BIN)

$(BIN) : $(OBJS)
	gcc $(cflags) -o $(BIN) $^

%.o : %.cpp
	gcc $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN)
