#
BIN=rad_pattern4_Rayl        

fflags = -e -g -Nl30 
FC=f77
cflags = -I$(X11_INCLUDE)

ldlibs =  -lmath   -lX11 -lm  -lev    -lio -lmap 

include $(DSAPMAKE)
DBG=-g

$(BIN) : $(BIN).o source.o surfread.o unwrap_contin.o intpol.o angles2tensor.o pha.o
	$(FC) $(FFLAGS) $(CFLAGS)  $(BIN).o source.o surfread.o unwrap_contin.o intpol.o pha.o  angles2tensor.o -o $(BIN) $(LDFLAGS) $(LDLIBS)
