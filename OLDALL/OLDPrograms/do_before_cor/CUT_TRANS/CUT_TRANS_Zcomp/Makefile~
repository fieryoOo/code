INST_DIR = $(HOME)/bin
BIN  = cut_trans_z

# need to include the math library when the math.h header file is used
cflags = -O -lm
CFLAGS = $(DBG) $(cflags)

FC = g77
CC = gcc

DBG = 
FOBJS = cut_trans_z.c

$(BIN) : $(FOBJS)
	$(CC) $(CFLAGS) $(FOBJS) -o $(BIN)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core
