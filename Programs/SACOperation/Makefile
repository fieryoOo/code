BIN1 = SAC_chhdr
BIN2 = SAC_filter
BIN3 = SAC_addall
BIN4 = SAC_dump
BIN5 = SAC_printHD
BIN6 = SAC_envelope
BIN7 = SAC_FTAN
BIN8 = SAC_mul
BIN9 = SAC_norm
BIN10 = SAC_Correlation
BIN11 = SAC_CrossCorr
BIN12 = SAC_Integrate
BIN13 = SAC_Amplitudes
BIN14 = SAC_ampratio
BIN15 = SAC_PHSM
BIN16 = SAC_taper
BIN17 = SAC_merge
BIN18 = SAC_cut
BIN18 = SAC_event

LIBS = -L${HOME}/usr/lib -lfftw3 -lSacRec

cflags = -lstdc++ -std=c++11 -O3

#OBJS = SacRec.o
OBJS = 

BINALL1 := $(BIN1) $(BIN2) ${BIN3} $(BIN4) $(BIN5) $(BIN6) $(BIN7) $(BIN8) $(BIN9) $(BIN10) 
BINALL1 := $(BINALL1) $(BIN11) $(BIN12) $(BIN13) $(BIN14) $(BIN16) $(BIN17) $(BIN18) $(BIN19)
BINALL = $(BINALL1) $(BIN15)
all : $(BINALL)

$(BIN15) : $(OBJS) $(BIN15).cpp
	g++ $(LIBS) -lmultivariate-splines-1-2 $(cflags) -o $@ $^

define makebin
$(1) : $(OBJS) $(1).cpp
	g++ $(LIBS) $(cflags) -o $$@ $$^
endef

$(foreach bin,$(BINALL1),$(eval $(call makebin,$(bin))))

%.o : %.cpp
	g++ $(cflags) -c $^

.PHONY : clean
clean:
	rm -f $(OBJS) $(BINALL)

.PHONY : install
install:
	cp $(BIN1) /home/yeti4009/bin
