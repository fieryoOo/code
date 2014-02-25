#
BIN=rad_pattern4_Love

fflags = -e -g -ffixed-line-length-132
# -Nl30
FC=f77
cflags = -I$(X11_INCLUDE)

ldlibs = -lmath -lX11 -lm -lev -lio -lmap 

OBJS = $(BIN).o source.o surfread.o unwrap_contin.o intpol.o pha.o  angles2tensor.o

include $(DSAPMAKE)
DBG=-g

$(BIN) : $(OBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $^ -o $(BIN) $(LDFLAGS) $(LDLIBS)

%.o : %.f
	$(FC) $(fflags) -c $^

clean :
	rm -f $(OBJS) $(BIN)
