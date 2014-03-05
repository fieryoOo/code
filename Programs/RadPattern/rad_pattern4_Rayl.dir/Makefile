BIN = rad_pattern4_Rayl

fflags = -e -g -ffixed-line-length-132
# -Nl30
FC = gfortran
cflags = -I$(X11_INCLUDE)

ldlibs = -lmath -lX11 -lm -lev -lio -lmap

OBJS = $(BIN).o source.o surfread.o angles2tensor.o pha.o new_unwrap.o 
#intpol.o
#unwrap_contin.o

include $(DSAPMAKE)
DBG=-g

$(BIN) : $(OBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $^ -o $(BIN) $(LDFLAGS) $(LDLIBS)

%.o : %.f
	$(FC) $(fflags) -c $^

clean :
	rm -f $(OBJS) $(BIN)
