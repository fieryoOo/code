BIN = RadPattern

VPATH = RadPattern_src

fflags = -e $(dbg) -ffixed-line-length-132
# -Nl30
cflags = $(dbg) -std=c++0x #-I$(X11_INCLUDE)

dbg = -g

FC = gfortran

CC = g++

ldlibs = -lstdc++ -lX11 -lm
# -lmath -lev -lio -lmap
OBJS = driver.o RadPattern.o $(VPATH)/rad_pattern4_Love.o $(VPATH)/rad_pattern4_Rayl.o $(VPATH)/source.o $(VPATH)/surfread.o \
       $(VPATH)/angles2tensor.o $(VPATH)/pha.o $(VPATH)/unwrapR.o $(VPATH)/unwrapL.o
#intpol.o unwrap_contin.o
#include $(DSAPMAKE)

$(BIN) : $(OBJS)
	$(FC) $^ -o $(BIN) $(ldlibs)

%.o : %.f
	$(FC) $(fflags) -c $< -o $@

%.o : %.cpp
	$(CC) $(cflags) -c $< -o $@

clean :
	rm -f $(OBJS) $(BIN)
