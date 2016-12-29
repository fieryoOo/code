BIN1 = WaveformFromPhaseDispersion
BIN2 = ComputePhaseShift
BIN3 = CurveDiff
BIN4 = DispersionFromSAC
BIN5 = PhaseToGroup
#BIN5 = SAC_envelope

LIBS = -L${HOME}/usr/lib -lfftw3 -lSacRec -lmultivariate-splines-1-2 -fopenmp

cflags = -lstdc++ -std=c++11 -g

#OBJS = SacRec.o
OBJS = 

all : $(BIN1) $(BIN2) $(BIN3) $(BIN4) $(BIN5)
#$(BIN2) ${BIN3} $(BIN4) $(BIN5)

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

#$(BIN2) : $(OBJS) $(BIN2).o
#	g++ $(LIBS) $(cflags) -o $(BIN2) $^

#$(BIN3) : $(OBJS) $(BIN3).o
#	g++ $(LIBS) $(cflags) -o $(BIN3) $^

%.o : %.cpp
	g++ $(cflags) -c $^

clean:
	rm -f $(OBJS) $(BIN1) $(BIN2) $(BIN3) $(BIN4) $(BIN5)

